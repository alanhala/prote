use std::{
    collections::{HashMap, HashSet},
    f64::MIN_POSITIVE,
    hash::Hash,
    usize, vec,
};

use nalgebra::constraint;

use crate::{
    atom_properties::{AtomProperties, Hybridization},
    atom_types::AtomTypes,
    bonding::{covalent_partners, Bonding},
    element::Element,
    geometry::{Plane3, Point3},
    mask::Mask,
    ring::Ring,
    selection::Selection,
    spatial_index::SpatialIndex,
    structure::Structure,
};

#[derive(Debug)]
pub enum Contact {
    VanDerWaals { atom_1: usize, atom_2: usize, distance: f64 },
    Hydrophobic { atom_1: usize, atom_2: usize, distance: f64 },
    MetalCoordination { atom_1: usize, atom_2: usize, distance: f64 },
    FaceToFace { ring_1: usize, ring_2: usize, distance: f64, plane_angle: f64, normal_to_center_angle: f64 },
    EdgeToFace { ring_1: usize, ring_2: usize, distance: f64, plane_angle: f64, normal_to_center_angle: f64, intersect_distance: f64 },
    HydrogenBond { donor: usize, acceptor: usize, distance: f64 },
    Ionic { cation: usize, anion: usize, distance: f64 },
    SaltBridge { cation: usize, anion: usize, distance: f64 },
    WaterBridge { atom_1: usize, atom_2: usize, water: usize, distance_1_to_water: f64, distance_2_to_water: f64 },
    CationPi { cation: usize, ring: usize, distance: f64, normal_to_cation_angle: f64 },
    Halogen { halogen: usize, acceptor: usize, distance: f64 },
}

struct AtomPairCheck {
    distance: f64,
    evaluate: fn(&Structure, &Bonding, &SpatialIndex, &AtomTypes, &AtomProperties, usize, usize, f64) -> Option<Contact>,
}

const ATOM_PAIR_CHECKS: [AtomPairCheck; 4] = [
    AtomPairCheck { distance: 4.5, evaluate: hydrophobic },
    AtomPairCheck { distance: 2.8, evaluate: metal_coordination },
    AtomPairCheck { distance: 4.5, evaluate: electrostatic },
    AtomPairCheck { distance: 3.5, evaluate: halogen },
];

fn hydrophobic(
    _: &Structure,
    _: &Bonding,
    _: &SpatialIndex,
    types: &AtomTypes,
    _: &AtomProperties,
    atom_1: usize,
    atom_2: usize,
    distance: f64,
) -> Option<Contact> {
    (types.hydrophobic.contains(atom_1) && types.hydrophobic.contains(atom_2)).then_some(Contact::Hydrophobic { atom_1, atom_2, distance })
}

fn halogen(
    structure: &Structure,
    bonding: &Bonding,
    index: &SpatialIndex,
    types: &AtomTypes,
    _: &AtomProperties,
    atom_1: usize,
    atom_2: usize,
    distance: f64,
) -> Option<Contact> {
    let (halogen, acceptor) = if types.halogen_bond_donor.contains(atom_1) && types.halogen_bond_acceptor.contains(atom_2) {
        Some((atom_1, atom_2))
    } else if types.halogen_bond_donor.contains(atom_2) && types.halogen_bond_acceptor.contains(atom_1) {
        Some((atom_2, atom_1))
    } else {
        None
    }?;
    let carbon_within_angle = covalent_partners(structure, bonding, halogen)
        .filter(|partner| partner.element == Element::C)
        .filter_map(|partner| partner.atom)
        .any(|carbon| {
            let carbon_halogen = index.positions[carbon] - index.positions[halogen];
            let acceptor_halogen = index.positions[acceptor] - index.positions[halogen];
            (130.0..=180.0).contains(&carbon_halogen.angle(&acceptor_halogen).to_degrees())
        });
    let neighbour_within_angle = covalent_partners(structure, bonding, acceptor).filter_map(|partner| partner.atom).any(|neighbour| {
        let halogen_acceptor = index.positions[halogen] - index.positions[acceptor];
        let neighbour_acceptor = index.positions[neighbour] - index.positions[acceptor];
        (80.0..=140.0).contains(&halogen_acceptor.angle(&neighbour_acceptor).to_degrees())
    });
    (neighbour_within_angle && carbon_within_angle).then(|| Contact::Halogen { halogen, acceptor, distance })
}

fn metal_coordination(
    _: &Structure,
    _: &Bonding,
    _: &SpatialIndex,
    types: &AtomTypes,
    _: &AtomProperties,
    atom_1: usize,
    atom_2: usize,
    distance: f64,
) -> Option<Contact> {
    let binds = |metal: usize, binder: usize| types.metal.contains(metal) && types.metal_binder.contains(binder);
    (binds(atom_1, atom_2) || binds(atom_2, atom_1)).then_some(Contact::MetalCoordination { atom_1, atom_2, distance })
}

fn electrostatic(
    structure: &Structure,
    bonding: &Bonding,
    index: &SpatialIndex,
    types: &AtomTypes,
    properties: &AtomProperties,
    atom_1: usize,
    atom_2: usize,
    distance: f64,
) -> Option<Contact> {
    let charged = opposite_charges(types, atom_1, atom_2);
    let bonded = if distance <= 3.5 { h_bond(structure, bonding, index, types, properties, atom_1, atom_2, distance) } else { None };
    match (charged, bonded) {
        (Some((cation, anion)), Some(_)) => Some(Contact::SaltBridge { cation, anion, distance }),
        (Some((cation, anion)), None) => Some(Contact::Ionic { cation, anion, distance }),
        (None, bonded) => bonded,
    }
}

fn opposite_charges(types: &AtomTypes, atom_1: usize, atom_2: usize) -> Option<(usize, usize)> {
    if types.cation.contains(atom_1) && types.anion.contains(atom_2) {
        Some((atom_1, atom_2))
    } else if types.anion.contains(atom_1) && types.cation.contains(atom_2) {
        Some((atom_2, atom_1))
    } else {
        None
    }
}

fn h_bond(
    structure: &Structure,
    bonding: &Bonding,
    index: &SpatialIndex,
    types: &AtomTypes,
    properties: &AtomProperties,
    atom_1: usize,
    atom_2: usize,
    distance: f64,
) -> Option<Contact> {
    let atom_1_type = (types.hydrogen_bond_acceptor.contains(atom_1), types.hydrogen_bond_donor.contains(atom_1));
    let atom_2_type = (types.hydrogen_bond_acceptor.contains(atom_2), types.hydrogen_bond_donor.contains(atom_2));
    let h_bond_roles: Vec<(usize, usize)> = match (atom_1_type, atom_2_type) {
        ((true, true), (true, true)) => vec![(atom_1, atom_2), (atom_2, atom_1)],
        ((true, _), (_, true)) => vec![(atom_1, atom_2)],
        ((_, true), (true, _)) => vec![(atom_2, atom_1)],
        _ => vec![],
    };
    h_bond_roles.into_iter().find_map(|(acceptor, donor)| {
        let (hydrogens, heavy): (Vec<usize>, Vec<usize>) = covalent_partners(structure, bonding, donor)
            .filter_map(|partner| partner.atom)
            .partition(|&atom| structure.atoms[atom].element == Element::H);
        let points_at_acceptor = |hydrogen: usize| {
            (index.positions[donor] - index.positions[hydrogen]).angle(&(index.positions[acceptor] - index.positions[hydrogen])).to_degrees() >= 130.0
        };
        if hydrogens.into_iter().any(points_at_acceptor) {
            return Some(Contact::HydrogenBond { donor, acceptor, distance });
        }
        // We checked all the explicit hydrogens and none passed the check
        if bonding.implicit_hydrogens[donor] == 0 {
            return None;
        }
        match (properties.hybridization[donor], heavy.len()) {
            (Some(Hybridization::Sp2), 2) => {
                let midpoint = Point3::centroid(&[index.positions[heavy[0]], index.positions[heavy[1]]]);
                let hydrogen_direction = index.positions[donor] - midpoint;
                ((index.positions[acceptor] - index.positions[donor]).angle(&hydrogen_direction).to_degrees() <= 35.0)
                    .then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            (Some(Hybridization::Sp2), 1) => {
                let neighbor_donor = index.positions[heavy[0]] - index.positions[donor];
                let acceptor_donor = index.positions[acceptor] - index.positions[donor];
                let angle = neighbor_donor.angle(&acceptor_donor).to_degrees();
                (angle >= 85.0 && angle <= 155.0).then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            (Some(Hybridization::Sp3), 1) => {
                let neighbor_donor = index.positions[heavy[0]] - index.positions[donor];
                let acceptor_donor = index.positions[acceptor] - index.positions[donor];
                let angle = neighbor_donor.angle(&acceptor_donor).to_degrees();
                (angle >= 75.0 && angle <= 145.0).then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            (Some(Hybridization::Sp3), 2) => {
                let midpoint = Point3::centroid(&[index.positions[heavy[0]], index.positions[heavy[1]]]);
                let away_from_neighbors = index.positions[donor] - midpoint;
                ((index.positions[acceptor] - index.positions[donor]).angle(&away_from_neighbors).to_degrees() <= 90.0)
                    .then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            (Some(Hybridization::Sp3), 3) => {
                let midpoint = Point3::centroid(&[index.positions[heavy[0]], index.positions[heavy[1]], index.positions[heavy[2]]]);
                let hydrogen_direction = index.positions[donor] - midpoint;
                ((index.positions[acceptor] - index.positions[donor]).angle(&hydrogen_direction).to_degrees() <= 35.0)
                    .then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            // Nothing to measure an angle against (water has no neighbor in the file).
            _ => Some(Contact::HydrogenBond { donor, acceptor, distance }),
        }
    })
}

fn ring_pair_contact(bonding: &Bonding, index: &SpatialIndex, ring_1: usize, ring_2: usize) -> Option<Contact> {
    let plane_1 = &bonding.rings[ring_1].plane(index.positions);
    let plane_2 = &bonding.rings[ring_2].plane(index.positions);
    let centre = plane_2.centroid - plane_1.centroid;
    let distance = plane_1.centroid.distance(&plane_2.centroid);
    let plane_angle = plane_1.angle(plane_2);
    let normal_to_center_angle = plane_1.normal_angle(&centre).min(plane_2.normal_angle(&centre)).to_degrees();
    let plane_angle = plane_angle.to_degrees();

    if distance <= 5.5 && plane_angle <= 35.0 && normal_to_center_angle <= 33.0 {
        return Some(Contact::FaceToFace { ring_1, ring_2, distance, plane_angle, normal_to_center_angle });
    }
    if distance <= 6.5 && plane_angle >= 50.0 && normal_to_center_angle <= 30.0 {
        // In a T, one ring's plane cuts through the middle of the other ring.
        let line = plane_1.intersection(plane_2)?;
        let intersect_distance = line.distance(&plane_1.centroid).min(line.distance(&plane_2.centroid));
        return (intersect_distance <= 1.5).then_some(Contact::EdgeToFace {
            ring_1,
            ring_2,
            distance,
            plane_angle,
            normal_to_center_angle,
            intersect_distance,
        });
    }
    None
}

pub fn find(
    structure: &Structure,
    bonding: &Bonding,
    types: &AtomTypes,
    properties: &AtomProperties,
    index: &SpatialIndex,
    query: &Selection,
    search: &Selection,
) -> Vec<Contact> {
    let (mut contacts, unexplained) = atom_pair_contacts(structure, bonding, types, properties, index, query, search);
    let rings = ring_contacts(bonding, types, index, query, search);
    let van_der_waals = van_der_waals(structure, bonding, &unexplained, &rings);
    contacts.extend(rings);
    contacts.extend(water_bridges(structure, bonding, types, properties, index, query, search));
    contacts.extend(van_der_waals);
    contacts
}

fn atom_pair_contacts(
    structure: &Structure,
    bonding: &Bonding,
    types: &AtomTypes,
    properties: &AtomProperties,
    index: &SpatialIndex,
    query: &Selection,
    search: &Selection,
) -> (Vec<Contact>, Vec<(usize, usize, f64)>) {
    let max_distance = ATOM_PAIR_CHECKS.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap().distance;
    let candidates = index
        .pairs_within(query, search, |atom_id| structure.atoms[atom_id].element.van_der_waals_radius().unwrap_or(0.0).max(max_distance / 2.0))
        .filter(|&(atom_1, atom_2, _)| !bonding.bonds.within_steps(atom_1, atom_2, 3));
    let mut contacts = vec![];
    let mut unexplained = vec![];
    for (atom_1, atom_2, distance) in candidates {
        let mut explained = false;
        for check in &ATOM_PAIR_CHECKS {
            if distance <= check.distance {
                if let Some(contact) = (check.evaluate)(structure, bonding, index, types, properties, atom_1, atom_2, distance) {
                    explained = true;
                    contacts.push(contact);
                }
            }
        }
        if !explained {
            unexplained.push((atom_1, atom_2, distance));
        }
    }
    (contacts, unexplained)
}

fn ring_contacts(bonding: &Bonding, types: &AtomTypes, index: &SpatialIndex, query: &Selection, search: &Selection) -> Vec<Contact> {
    let mut contacts = vec![];
    let aromatic: Mask<Ring> = Mask::from_bools(bonding.rings.iter().map(|ring| ring.aromatic).collect());
    let query_rings: Mask<Ring> =
        aromatic.and(&Mask::from_bools(bonding.rings.iter().map(|ring| ring.atoms.iter().all(|&atom| query.contains(atom))).collect()));
    let search_rings: Mask<Ring> =
        aromatic.and(&Mask::from_bools(bonding.rings.iter().map(|ring| ring.atoms.iter().all(|&atom| search.contains(atom))).collect()));
    for ring_1 in query_rings.indices() {
        for ring_2 in search_rings.indices() {
            let found_from_ring_2 = ring_2 < ring_1 && query_rings.contains(ring_2) && search_rings.contains(ring_1);
            let shares_atoms = bonding.rings[ring_1].atoms.iter().any(|atom| bonding.rings[ring_2].atoms.contains(atom));
            if ring_1 == ring_2 || found_from_ring_2 || shares_atoms {
                continue;
            }
            contacts.extend(ring_pair_contact(bonding, index, ring_1, ring_2));
        }
    }

    for ring in aromatic.indices() {
        for cation in types.cation.indices() {
            let wanted = (query_rings.contains(ring) && search.contains(cation)) || (query.contains(cation) && search_rings.contains(ring));
            if !wanted {
                continue;
            }

            if bonding.rings[ring].atoms.contains(&cation) {
                continue;
            }

            let plane = bonding.rings[ring].plane(index.positions);
            let to_cation = index.positions[cation] - plane.centroid;
            let distance = to_cation.magnitude();
            let normal_to_cation_angle = plane.normal_angle(&to_cation).to_degrees();
            if distance <= 4.5 && normal_to_cation_angle <= 30.0 {
                contacts.push(Contact::CationPi { cation, ring, distance, normal_to_cation_angle });
            }
        }
    }
    contacts
}

fn water_bridges(
    structure: &Structure,
    bonding: &Bonding,
    types: &AtomTypes,
    properties: &AtomProperties,
    index: &SpatialIndex,
    query: &Selection,
    search: &Selection,
) -> Vec<Contact> {
    let mut contacts = vec![];
    let bonded_to_water = |selection: &Selection| -> Vec<(usize, usize, f64)> {
        index
            .pairs_closer_than(selection, &types.water, 3.5)
            .filter(|&(atom, water, distance)| {
                !types.water.contains(atom) && h_bond(structure, bonding, index, types, properties, atom, water, distance).is_some()
            })
            .collect()
    };
    let mut search_by_water: HashMap<usize, Vec<(usize, f64)>> = HashMap::new();
    for (atom, water, distance) in bonded_to_water(search) {
        search_by_water.entry(water).or_default().push((atom, distance));
    }
    for (atom_1, water, distance_1) in bonded_to_water(query) {
        for &(atom_2, distance_2) in search_by_water.get(&water).into_iter().flatten() {
            if atom_1 == atom_2 {
                continue;
            }
            if atom_2 < atom_1 && query.contains(atom_2) && search.contains(atom_1) {
                continue;
            }
            if bonding.bonds.within_steps(atom_1, atom_2, 3) {
                continue;
            }
            let angle = (index.positions[atom_1] - index.positions[water]).angle(&(index.positions[atom_2] - index.positions[water])).to_degrees();
            if (71.0..=140.0).contains(&angle) {
                contacts.push(Contact::WaterBridge { atom_1, atom_2, water, distance_1_to_water: distance_1, distance_2_to_water: distance_2 });
            }
        }
    }
    contacts
}

/// The leftover kind: touching pairs that no atom pair check explained and that are not part of a ring contact.
fn van_der_waals(structure: &Structure, bonding: &Bonding, unexplained: &[(usize, usize, f64)], ring_contacts: &[Contact]) -> Vec<Contact> {
    let mut ring_pairs: HashSet<(usize, usize)> = HashSet::new();
    for contact in ring_contacts {
        match contact {
            Contact::FaceToFace { ring_1, ring_2, .. } | Contact::EdgeToFace { ring_1, ring_2, .. } => {
                for &atom_1 in &bonding.rings[*ring_1].atoms {
                    for &atom_2 in &bonding.rings[*ring_2].atoms {
                        ring_pairs.insert(unordered(atom_1, atom_2));
                    }
                }
            }
            Contact::CationPi { cation, ring, .. } => {
                for &ring_atom in &bonding.rings[*ring].atoms {
                    ring_pairs.insert(unordered(ring_atom, *cation));
                }
            }
            _ => {}
        }
    }
    unexplained
        .iter()
        .filter(|&&(atom_1, atom_2, distance)| !ring_pairs.contains(&unordered(atom_1, atom_2)) && touching(structure, atom_1, atom_2, distance))
        .map(|&(atom_1, atom_2, distance)| Contact::VanDerWaals { atom_1, atom_2, distance })
        .collect()
}

fn unordered(atom_1: usize, atom_2: usize) -> (usize, usize) {
    if atom_1 < atom_2 {
        (atom_1, atom_2)
    } else {
        (atom_2, atom_1)
    }
}

fn touching(structure: &Structure, atom_1: usize, atom_2: usize, distance: f64) -> bool {
    let radius = |atom: usize| structure.atoms[atom].element.van_der_waals_radius();
    radius(atom_1).zip(radius(atom_2)).is_some_and(|(radius_1, radius_2)| distance <= radius_1 + radius_2)
}
