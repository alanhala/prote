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
    bond::BondKind,
    bond_graph,
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
    // later: CationPi { cation, ring, .. }
}

struct AtomPairCheck {
    distance: f64,
    evaluate: fn(&Structure, &SpatialIndex, &AtomTypes, &AtomProperties, usize, usize, f64) -> Option<Contact>,
}

const ATOM_PAIR_CHECKS: [AtomPairCheck; 3] = [
    AtomPairCheck { distance: 4.5, evaluate: hydrophobic },
    AtomPairCheck { distance: 2.8, evaluate: metal_coordination },
    AtomPairCheck { distance: 3.5, evaluate: h_bond },
];

fn hydrophobic(
    _: &Structure,
    _: &SpatialIndex,
    types: &AtomTypes,
    _: &AtomProperties,
    atom_1: usize,
    atom_2: usize,
    distance: f64,
) -> Option<Contact> {
    (types.hydrophobic.contains(atom_1) && types.hydrophobic.contains(atom_2)).then_some(Contact::Hydrophobic { atom_1, atom_2, distance })
}

fn metal_coordination(
    _: &Structure,
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

fn h_bond(
    structure: &Structure,
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
        let present: Vec<usize> =
            structure.covalent_partners(donor).filter(|partner| partner.element != Element::H).filter_map(|partner| partner.atom).collect();
        match (properties.hybridization[donor], present.len()) {
            (Some(Hybridization::Sp2), 2) => {
                let midpoint = Point3::centroid(&[index.positions[present[0]], index.positions[present[1]]]);
                let hydrogen_direction = index.positions[donor] - midpoint;
                ((index.positions[acceptor] - index.positions[donor]).angle(&hydrogen_direction).to_degrees() <= 35.0)
                    .then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            (Some(Hybridization::Sp2), 1) => {
                let neighbor_donor = index.positions[present[0]] - index.positions[donor];
                let acceptor_donor = index.positions[acceptor] - index.positions[donor];
                let angle = neighbor_donor.angle(&acceptor_donor).to_degrees();
                (angle >= 85.0 && angle <= 155.0).then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            (Some(Hybridization::Sp3), 1) => {
                let neighbor_donor = index.positions[present[0]] - index.positions[donor];
                let acceptor_donor = index.positions[acceptor] - index.positions[donor];
                let angle = neighbor_donor.angle(&acceptor_donor).to_degrees();
                (angle >= 75.0 && angle <= 145.0).then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            (Some(Hybridization::Sp3), 2) => {
                let midpoint = Point3::centroid(&[index.positions[present[0]], index.positions[present[1]]]);
                let away_from_neighbors = index.positions[donor] - midpoint;
                ((index.positions[acceptor] - index.positions[donor]).angle(&away_from_neighbors).to_degrees() <= 90.0)
                    .then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            (Some(Hybridization::Sp3), 3) => {
                let midpoint = Point3::centroid(&[index.positions[present[0]], index.positions[present[1]], index.positions[present[2]]]);
                let hydrogen_direction = index.positions[donor] - midpoint;
                ((index.positions[acceptor] - index.positions[donor]).angle(&hydrogen_direction).to_degrees() <= 35.0)
                    .then(|| Contact::HydrogenBond { donor, acceptor, distance })
            }
            _ => Some(Contact::HydrogenBond { donor, acceptor, distance }),
        }
    })
}

fn ring_pair_contact(structure: &Structure, index: &SpatialIndex, ring_1: usize, ring_2: usize) -> Option<Contact> {
    let plane_1 = &structure.bonding.rings[ring_1].plane(index.positions);
    let plane_2 = &structure.bonding.rings[ring_2].plane(index.positions);
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
    types: &AtomTypes,
    properties: &AtomProperties,
    index: &SpatialIndex,
    query: &Selection,
    search: &Selection,
) -> Vec<Contact> {
    let max_distance = ATOM_PAIR_CHECKS.iter().max_by(|a, b| a.distance.partial_cmp(&b.distance).unwrap()).unwrap().distance;
    let mut contacts = vec![];
    let candidates: Vec<(usize, usize, f64)> = index
        .pairs_within(query, search, |atom_id| structure.atoms[atom_id].element.van_der_waals_radius().unwrap_or(0.0).max(max_distance / 2.0))
        .filter(|&(atom_1, atom_2, _)| !structure.bonding.bonds.within_steps(atom_1, atom_2, 3))
        .collect();
    let mut explained = vec![false; candidates.len()];
    let mut ring_pairs: HashSet<(usize, usize)> = HashSet::new();
    for (i, &(atom_1, atom_2, distance)) in candidates.iter().enumerate() {
        for check in &ATOM_PAIR_CHECKS {
            if distance <= check.distance {
                if let Some(contact) = (check.evaluate)(structure, index, types, properties, atom_1, atom_2, distance) {
                    explained[i] = true;
                    contacts.push(contact);
                }
            }
        }
    }

    let aromatic: Mask<Ring> = Mask::from_bools(structure.bonding.rings.iter().map(|ring| ring.aromatic).collect());
    let query_rings: Mask<Ring> =
        aromatic.and(&Mask::from_bools(structure.bonding.rings.iter().map(|ring| ring.atoms.iter().all(|&atom| query.contains(atom))).collect()));
    let search_rings: Mask<Ring> =
        aromatic.and(&Mask::from_bools(structure.bonding.rings.iter().map(|ring| ring.atoms.iter().all(|&atom| search.contains(atom))).collect()));
    for ring_1 in query_rings.indices() {
        for ring_2 in search_rings.indices() {
            let found_from_ring_2 = ring_2 < ring_1 && query_rings.contains(ring_2) && search_rings.contains(ring_1);
            let shares_atoms = structure.bonding.rings[ring_1].atoms.iter().any(|atom| structure.bonding.rings[ring_2].atoms.contains(atom));
            if ring_1 == ring_2 || found_from_ring_2 || shares_atoms {
                continue;
            }
            if let Some(contact) = ring_pair_contact(structure, index, ring_1, ring_2) {
                contacts.push(contact);
                for &atom_1 in &structure.bonding.rings[ring_1].atoms {
                    for &atom_2 in &structure.bonding.rings[ring_2].atoms {
                        ring_pairs.insert(unordered(atom_1, atom_2));
                    }
                }
            }
        }
    }

    for (i, &(atom_1, atom_2, distance)) in candidates.iter().enumerate() {
        if !explained[i] && !ring_pairs.contains(&unordered(atom_1, atom_2)) && touching(structure, atom_1, atom_2, distance) {
            contacts.push(Contact::VanDerWaals { atom_1, atom_2, distance });
        }
    }
    contacts
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
