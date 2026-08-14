use std::{iter::Cycle, vec};

use crate::{
    bond::Bond, bond_angle::BondAngle, bond_order::BondOrderSource, conformer::Conformer, proximity::Proximity,
    proximity_graph::ProximityGraph, spatial_index::SpatialIndex,
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
        let proximity_graph = ProximityGraph::new(search, spatial_index, query, Proximity::Bonded);
        let conformers = proximity_graph.conformers;
        let bonds: Vec<Vec<Bond>> = proximity_graph
            .edges
            .into_iter()
            .map(|edges| {
                edges
                    .into_iter()
                    .map(|(atom_1, atom_2, distance)| Bond::new(atom_1, atom_2, distance))
                    .collect()
            })
            .collect();
        let mut angles = (0..bonds.len()).map(|_| Vec::new()).collect::<Vec<Vec<BondAngle>>>();
        for (atom_index, atom_bonds) in bonds.iter().enumerate() {
            match atom_bonds.len() {
                0 | 1 => {}
                bonds_len => {
                    let vertex_atom_pointer = atom_bonds[0].atom_1;
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

    pub fn dfs(&self) {
        let mut visited = vec![false; self.bonds.len()];
        let mut on_path = vec![false; self.bonds.len()];
        let mut cycles: Vec<Vec<usize>> = vec![];
        let mut path: Vec<usize> = vec![];
        for i in 0..self.bonds.len() {
            if !visited[i] && self.hybridization_at(i) == Some(Hybridization::SP2) {
                self.visit(i, None, &mut visited, &mut path, &mut on_path, &mut cycles);
            }
        }
        println!("{:?}", cycles)
    }

    fn visit(
        &self,
        node: usize,
        parent: Option<usize>,
        visited: &mut Vec<bool>,
        path: &mut Vec<usize>,
        on_path: &mut Vec<bool>,
        cycles: &mut Vec<Vec<usize>>,
    ) {
        visited[node] = true;
        for bond in &self.bonds[node] {
            let neighbor = bond.atom_2.index();
            if !visited[neighbor] && self.hybridization_at(neighbor) == Some(Hybridization::SP2) {
                path.push(node);
                on_path[node] = true;
                self.visit(neighbor, Some(node), visited, path, on_path, cycles);
                path.pop();
                on_path[node] = false;
            } else if parent == Some(neighbor) {
                continue;
            } else if on_path[neighbor] {
                let idx = path.iter().position(|&x| x == neighbor).unwrap();
                let mut cycle = path[idx..].to_vec();
                cycle.push(node);
                cycles.push(cycle);
            }
        }
    }

    pub fn to_dot(&self) -> String {
        let mut out = String::from("graph G {\n");
        for (i, bonds) in self.bonds.iter().enumerate() {
            for bond in bonds {
                let j = bond.atom_2.index();
                if i < j {
                    out += &format!("  {i} -- {j};\n");
                }
            }
        }
        out += "}\n";
        out
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
