use std::vec;

use crate::{
    atom_pointer::AtomPointer,
    bond::Bond,
    bond_angle::BondAngle,
    bond_order::{BondOrder, BondOrderSource},
    conformer::Conformer,
    proximity::Proximity,
    spatial_index::SpatialIndex,
};

const ANGLE_TOLERANCE: f64 = 5.0;

#[derive(Debug, PartialEq)]
pub enum Hybridization {
    SP,
    SP2,
    SP3,
}

enum AngleShape {
    Linear,
    Trigonal,
    Tetrahedral,
    Other,
}

fn categorize_angle(angle: f64) -> AngleShape {
    if (angle - 180.0).abs() <= ANGLE_TOLERANCE {
        AngleShape::Linear
    } else if (angle - 120.0).abs() <= ANGLE_TOLERANCE {
        AngleShape::Trigonal
    } else if (angle - 109.5).abs() <= ANGLE_TOLERANCE {
        AngleShape::Tetrahedral
    } else {
        AngleShape::Other
    }
}

#[derive(Debug)]
pub struct BondGraph<'a> {
    pub conformers: Vec<&'a Conformer<'a>>,
    pub bonds: Vec<Vec<Bond>>,
    pub angles: Vec<Vec<BondAngle>>,
}

impl<'a> BondGraph<'a> {
    pub fn new(search: &'a Conformer, spatial_index: &SpatialIndex, query: &'a Conformer) -> Self {
        let same_conformer = search == query;
        let query_conformer_index = if same_conformer { 0 } else { 1 };
        let conformers = if same_conformer {
            vec![search]
        } else {
            vec![search, query]
        };
        let query_offset = if same_conformer { 0 } else { search.topology.atoms.len() };
        let atoms = if same_conformer {
            search.topology.atoms.len()
        } else {
            search.topology.atoms.len() + query.topology.atoms.len()
        };
        let mut bonds = (0..atoms).map(|_| Vec::new()).collect::<Vec<Vec<Bond>>>();
        let mut angles = (0..atoms).map(|_| Vec::new()).collect::<Vec<Vec<BondAngle>>>();
        for (query_atom_index, query_atom) in query.topology.atoms.iter().enumerate() {
            let candidates = spatial_index.candidates_within(
                &search.frame.positions,
                &query.frame.positions[query_atom_index],
                Proximity::search_radius(Proximity::Bonded, &query_atom),
            );
            for (search_atom_index, distance) in candidates {
                let query_atom_pointer = AtomPointer::new(query_conformer_index, query_offset, query_atom_index);
                let search_atom = &search.topology.atoms[search_atom_index];
                let search_atom_pointer = AtomPointer::new(0, 0, search_atom_index);
                let different_atoms = search_atom_pointer != query_atom_pointer;
                let bonded = Proximity::classify(&search_atom, &query_atom, distance) == Some(Proximity::Bonded);
                if different_atoms && bonded {
                    if !bonds[search_atom_pointer.index()]
                        .iter()
                        .any(|b| b.atom_2 == query_atom_pointer)
                    {
                        bonds[search_atom_pointer.index()].push(Bond::new(
                            search_atom_pointer,
                            query_atom_pointer,
                            distance,
                        ));
                    }
                    if !bonds[query_atom_pointer.index()]
                        .iter()
                        .any(|b| b.atom_2 == search_atom_pointer)
                    {
                        bonds[query_atom_pointer.index()].push(Bond::new(
                            query_atom_pointer,
                            search_atom_pointer,
                            distance,
                        ));
                    }
                }
            }
        }
        for (atom_index, atom_bonds) in bonds.iter().enumerate() {
            match atom_bonds.len() {
                0 | 1 => {}
                bonds_len => {
                    let vertex_atom_pointer = AtomPointer::new(
                        atom_bonds[0].atom_1.conformer_index,
                        atom_bonds[0].atom_1.conformer_offset,
                        atom_bonds[0].atom_1.atom_index,
                    );
                    let vertex_position = &conformers[vertex_atom_pointer.conformer_index].frame.positions
                        [vertex_atom_pointer.atom_index];
                    for i in 0..(bonds_len - 1) {
                        for j in (i + 1)..bonds_len {
                            let atom_pointer_1 = &atom_bonds[i].atom_2;
                            let atom_pointer_2 = &atom_bonds[j].atom_2;
                            let atom_1_position =
                                &conformers[atom_pointer_1.conformer_index].frame.positions[atom_pointer_1.atom_index];
                            let atom_2_position =
                                &conformers[atom_pointer_2.conformer_index].frame.positions[atom_pointer_2.atom_index];
                            angles[atom_index].push(BondAngle::new(
                                (*atom_pointer_1, atom_1_position),
                                (vertex_atom_pointer, vertex_position),
                                (*atom_pointer_2, atom_2_position),
                            ))
                        }
                    }
                }
            }
        }
        Self {
            conformers,
            bonds,
            angles,
        }
    }

    pub fn hybridization_at(&self, atom_index: usize) -> Option<Hybridization> {
        let atom_angles = self.angles.get(atom_index)?;
        match atom_angles.as_slice() {
            [a] => match categorize_angle(a.value) {
                AngleShape::Linear => Some(Hybridization::SP),
                AngleShape::Trigonal => Some(Hybridization::SP2),
                AngleShape::Tetrahedral => Some(Hybridization::SP3),
                AngleShape::Other => None,
            },
            [a, b, c] => match (
                categorize_angle(a.value),
                categorize_angle(b.value),
                categorize_angle(c.value),
            ) {
                (AngleShape::Tetrahedral, AngleShape::Tetrahedral, AngleShape::Tetrahedral) => Some(Hybridization::SP3),
                (AngleShape::Trigonal, AngleShape::Trigonal, AngleShape::Trigonal)
                    if (a.value + b.value + c.value - 360.0).abs() <= ANGLE_TOLERANCE =>
                {
                    Some(Hybridization::SP2)
                }
                _ => None,
            },
            [a, b, c, d, e, f] => {
                let all_tetrahedral = [a, b, c, d, e, f]
                    .iter()
                    .all(|angle| matches!(categorize_angle(angle.value), AngleShape::Tetrahedral));
                if all_tetrahedral {
                    Some(Hybridization::SP3)
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    // TODO: We are computing bound order in both directions. Unnecessary.
    pub fn assign_orders(&mut self, sources: &[&dyn BondOrderSource]) {
        for bonds in self.bonds.iter_mut() {
            for bond in bonds.iter_mut() {
                bond.order = sources
                    .iter()
                    .find_map(|source| source.bond_order(bond, &self.conformers))
            }
        }
    }

    // fn leftover_capacity(atom: &Atom, bonds: &[Bond]) -> Option<u8> {
    //     let used: u8 = bonds.iter().filter_map(|bond| bond.order).map(BondOrder::as_u8).sum();
    //     atom.typical_valence().map(|typical| typical.saturating_sub(used))
    // }

    // pub fn assign_remaining_orders(&mut self, ensembles: &[Ensemble]) {
    //     let mut updates: Vec<(usize, usize, BondOrder)> = Vec::new();
    //     for (search_index, bonds) in self.bonds.iter().enumerate() {
    //         for (bond_index, bond) in bonds.iter().enumerate() {
    //             if bond.order.is_some() {
    //                 continue;
    //             }
    //             let (atom_1, _) = resolve(ensembles, &bond.atom_1);
    //             let (atom_2, _) = resolve(ensembles, &bond.atom_2);
    //             let leftover_1 = Self::leftover_capacity(atom_1, &self.bonds[bond.atom_1.index]);
    //             let leftover_2 = Self::leftover_capacity(atom_2, &self.bonds[bond.atom_2.index]);
    //             let order = leftover_1
    //                 .zip(leftover_2)
    //                 .and_then(|(a, b)| BondOrder::from_u8(a.min(b)));
    //             if let Some(order) = order {
    //                 updates.push((search_index, bond_index, order));
    //             }
    //         }
    //     }
    //     for (search_index, bond_index, order) in updates {
    //         self.bonds[search_index][bond_index].order = Some(order);
    //     }
    // }
}
