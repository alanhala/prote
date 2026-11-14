use std::{
    collections::{HashMap, HashSet},
    f64::MIN_POSITIVE,
    hash::Hash,
};

use nalgebra::constraint;

use crate::{
    atom_types::AtomTypes, bond::BondKind, bond_graph, geometry::Plane3, mask::Mask, ring::Ring, selection::Selection, spatial_index::SpatialIndex,
    structure::Structure,
};

#[derive(Debug)]
pub enum Contact {
    VanDerWaals { atom_1: usize, atom_2: usize, distance: f64 },
    Hydrophobic { atom_1: usize, atom_2: usize, distance: f64 },
    MetalCoordination { atom_1: usize, atom_2: usize, distance: f64 },
    FaceToFace { ring_1: usize, ring_2: usize, distance: f64, plane_angle: f64, normal_to_center_angle: f64 },
    EdgeToFace { ring_1: usize, ring_2: usize, distance: f64, plane_angle: f64, normal_to_center_angle: f64, intersect_distance: f64 },
    // later: HydrogenBond { donor, acceptor, .. }, CationPi { cation, ring, .. }
}

struct AtomPairCheck {
    distance: f64,
    evaluate: fn(&AtomTypes, usize, usize, f64) -> Option<Contact>,
}

const ATOM_PAIR_CHECKS: [AtomPairCheck; 2] =
    [AtomPairCheck { distance: 4.5, evaluate: hydrophobic }, AtomPairCheck { distance: 2.8, evaluate: metal_coordination }];

fn hydrophobic(types: &AtomTypes, atom_1: usize, atom_2: usize, distance: f64) -> Option<Contact> {
    (types.hydrophobic.contains(atom_1) && types.hydrophobic.contains(atom_2)).then_some(Contact::Hydrophobic { atom_1, atom_2, distance })
}

fn metal_coordination(types: &AtomTypes, atom_1: usize, atom_2: usize, distance: f64) -> Option<Contact> {
    let binds = |metal: usize, binder: usize| types.metal.contains(metal) && types.metal_binder.contains(binder);
    (binds(atom_1, atom_2) || binds(atom_2, atom_1)).then_some(Contact::MetalCoordination { atom_1, atom_2, distance })
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

pub fn find(structure: &Structure, types: &AtomTypes, index: &SpatialIndex, query: &Selection, search: &Selection) -> Vec<Contact> {
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
                if let Some(contact) = (check.evaluate)(types, atom_1, atom_2, distance) {
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
