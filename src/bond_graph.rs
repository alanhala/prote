// use std::{collections::VecDeque, vec};

// use crate::{
//     atom::Atom,
//     atom_pointer::{self, AtomPointer},
//     bond::Bond,
//     bond_angle::BondAngle,
//     bond_order::{BondOrder, BondOrderSource},
//     conformer::{self, Conformer},
//     proximity::Proximity,
//     proximity_graph_old::ProximityGraph,
//     ring::Ring,
//     spatial_index::SpatialIndex,
// };

// const ANGLE_TOLERANCE: f64 = 5.0;

// #[derive(Debug, PartialEq)]
// pub enum Hybridization {
//     SP,
//     SP2,
//     SP3,
// }

// enum AngleShape {
//     Linear,
//     Trigonal,
//     Tetrahedral,
//     Other,
// }

// fn categorize_angle(angle: f64) -> AngleShape {
//     if (angle - 180.0).abs() <= ANGLE_TOLERANCE {
//         AngleShape::Linear
//     } else if (angle - 120.0).abs() <= ANGLE_TOLERANCE {
//         AngleShape::Trigonal
//     } else if (angle - 109.5).abs() <= ANGLE_TOLERANCE {
//         AngleShape::Tetrahedral
//     } else {
//         AngleShape::Other
//     }
// }

// #[derive(Debug)]
// pub struct BondGraph<'a> {
//     pub conformers: Vec<&'a Conformer<'a>>,
//     pub bonds: Vec<Vec<Bond>>,
//     pub angles: Vec<Vec<BondAngle>>,
//     // TODO: `rings` can contain the same ring more than once — a ring found directly by DFS
//     // and the same ring produced again as a split_cycle() piece of a different oversized
//     // cycle aren't deduplicated against each other. Harmless for hybridization_at (it only
//     // checks whether *any* referenced ring is coplanar), but will double-count real rings
//     // once something iterates `rings` directly (aromatic stacking, cation-pi). Dedup by
//     // sorted atom-index set before that happens.
//     pub rings: Vec<Ring>,
//     pub atom_rings: Vec<Vec<usize>>,
// }

// impl<'a> BondGraph<'a> {
//     pub fn new(search: &'a Conformer, spatial_index: &SpatialIndex, query: &'a Conformer) -> Self {
//         let proximity_graph = ProximityGraph::new(search, spatial_index, query, Proximity::Bonded);
//         let conformers = proximity_graph.conformers;
//         let bonds: Vec<Vec<Bond>> = proximity_graph
//             .edges
//             .into_iter()
//             .map(|edges| {
//                 edges
//                     .into_iter()
//                     .map(|(atom_1, atom_2, distance)| Bond::new(atom_1, atom_2, distance))
//                     .collect()
//             })
//             .collect();
//         let mut angles = (0..bonds.len()).map(|_| Vec::new()).collect::<Vec<Vec<BondAngle>>>();
//         for (atom_index, atom_bonds) in bonds.iter().enumerate() {
//             match atom_bonds.len() {
//                 0 | 1 => {}
//                 bonds_len => {
//                     let vertex_atom_pointer = atom_bonds[0].atom_1;
//                     let vertex_position = &conformers[vertex_atom_pointer.conformer_index].frame.positions
//                         [vertex_atom_pointer.atom_index];
//                     for i in 0..(bonds_len - 1) {
//                         for j in (i + 1)..bonds_len {
//                             let atom_pointer_1 = &atom_bonds[i].atom_2;
//                             let atom_pointer_2 = &atom_bonds[j].atom_2;
//                             let atom_1_position =
//                                 &conformers[atom_pointer_1.conformer_index].frame.positions[atom_pointer_1.atom_index];
//                             let atom_2_position =
//                                 &conformers[atom_pointer_2.conformer_index].frame.positions[atom_pointer_2.atom_index];
//                             angles[atom_index].push(BondAngle::new(
//                                 (*atom_pointer_1, atom_1_position),
//                                 (vertex_atom_pointer, vertex_position),
//                                 (*atom_pointer_2, atom_2_position),
//                             ))
//                         }
//                     }
//                 }
//             }
//         }
//         let (rings, atom_rings) = Self::find_rings(&bonds, &conformers);
//         Self {
//             conformers,
//             bonds,
//             angles,
//             rings,
//             atom_rings,
//         }
//     }

//     pub fn hybridization_at(&self, atom_index: usize) -> Option<Hybridization> {
//         if self.atom_rings[atom_index]
//             .iter()
//             .any(|&ring_index| self.rings[ring_index].coplanar)
//         {
//             return Some(Hybridization::SP2);
//         }
//         let atom_angles = self.angles.get(atom_index)?;

//         match atom_angles.as_slice() {
//             [a] => match categorize_angle(a.value) {
//                 AngleShape::Linear => Some(Hybridization::SP),
//                 AngleShape::Trigonal => Some(Hybridization::SP2),
//                 AngleShape::Tetrahedral => Some(Hybridization::SP3),
//                 AngleShape::Other => None,
//             },
//             [a, b, c] => match (
//                 categorize_angle(a.value),
//                 categorize_angle(b.value),
//                 categorize_angle(c.value),
//             ) {
//                 (AngleShape::Tetrahedral, AngleShape::Tetrahedral, AngleShape::Tetrahedral) => Some(Hybridization::SP3),
//                 (AngleShape::Trigonal, AngleShape::Trigonal, AngleShape::Trigonal)
//                     if (a.value + b.value + c.value - 360.0).abs() <= ANGLE_TOLERANCE =>
//                 {
//                     Some(Hybridization::SP2)
//                 }
//                 _ => None,
//             },
//             [a, b, c, d, e, f] => {
//                 let all_tetrahedral = [a, b, c, d, e, f]
//                     .iter()
//                     .all(|angle| matches!(categorize_angle(angle.value), AngleShape::Tetrahedral));
//                 if all_tetrahedral {
//                     Some(Hybridization::SP3)
//                 } else {
//                     None
//                 }
//             }
//             _ => None,
//         }
//     }

//     /// The atom behind a graph index (a conformer offset plus that conformer's atom index).
//     pub fn atom_at(&self, index: usize) -> &Atom {
//         let mut offset = 0;
//         for conformer in &self.conformers {
//             let atoms = &conformer.topology.atoms;
//             if index < offset + atoms.len() {
//                 return &atoms[index - offset];
//             }
//             offset += atoms.len();
//         }
//         panic!("graph index {index} is past the last atom")
//     }

//     /// The bond from `atom_1` to `atom_2`, if they are bonded.
//     pub fn bond_between(&self, atom_1: usize, atom_2: usize) -> Option<&Bond> {
//         self.bonds[atom_1].iter().find(|bond| bond.atom_2.index() == atom_2)
//     }

//     /// Sets the order on both stored directions of the bond between `atom_1` and `atom_2`.
//     pub fn set_order(&mut self, atom_1: usize, atom_2: usize, order: BondOrder) {
//         for (from, to) in [(atom_1, atom_2), (atom_2, atom_1)] {
//             if let Some(bond) = self.bonds[from].iter_mut().find(|bond| bond.atom_2.index() == to) {
//                 bond.order = Some(order);
//             }
//         }
//     }

//     pub fn to_dot(&self) -> String {
//         let mut dot = String::from("graph BondGraph {\n");
//         for atom_index in 0..self.bonds.len() {
//             let atom = self.atom_at(atom_index);
//             let (style, fillcolor) = if self.atom_rings[atom_index].is_empty() {
//                 ("filled", "white")
//             } else {
//                 ("filled", "lightblue")
//             };
//             dot.push_str(&format!(
//                 "  {atom_index} [label=\"{} {} {}\", style={style}, fillcolor={fillcolor}];\n",
//                 atom_index,
//                 atom.symbol(),
//                 atom.name
//             ));
//         }
//         for atom_bonds in &self.bonds {
//             for bond in atom_bonds {
//                 let (i, j) = (bond.atom_1.index(), bond.atom_2.index());
//                 if i < j {
//                     let attributes = match bond.order {
//                         Some(BondOrder::Double) => " [penwidth=3]",
//                         Some(BondOrder::Triple) => " [penwidth=5]",
//                         Some(BondOrder::Aromatic) => " [penwidth=2, style=dashed]",
//                         _ => "",
//                     };
//                     dot.push_str(&format!("  {i} -- {j}{attributes};\n"));
//                 }
//             }
//         }
//         dot.push_str("}\n");
//         dot
//     }

//     // TODO: We are computing bound order in both directions. Unnecessary.
//     pub fn assign_orders(&mut self, sources: &[&dyn BondOrderSource]) {
//         for bonds in self.bonds.iter_mut() {
//             for bond in bonds.iter_mut() {
//                 bond.order = sources
//                     .iter()
//                     .find_map(|source| source.bond_order(bond, &self.conformers))
//             }
//         }
//     }

//     // TODO: This DFS-based cycle basis is *a* valid cycle basis, but for graphs where three or
//     // more rings are all fused together sharing multiple edges (e.g. steroid ring systems,
//     // porphyrins/heme), it isn't guaranteed to match the SSSR (smallest set of smallest rings)
//     // a chemist would draw. Fine for standard amino acid side chains (at most 2 fused rings),
//     // revisit if this ever runs against ligands/cofactors with more complex fused ring systems.
//     // To research if/when this needs fixing: "minimum cycle basis" (MCB) is the graph-theory
//     // name for the general problem; Horton's algorithm (1987) is the classic polynomial-time
//     // solution, later improved on by de Pina's algorithm. Figueras' 1996 ring-perception paper
//     // (breadth-first search based) is the more chemistry-specific version that OpenBabel/RDKit's
//     // SSSR implementations descend from. RING 4.0 covers this in the specific context of
//     // protein/ligand ring perception.
//     //
//     // Atoms with no typical valence (metals, in practice) are never part of a ring: a chelate
//     // loop such as N-Fe-N through a heme iron is not a ring in the chemical sense, and letting
//     // the search run through the metal both fabricates coplanar "rings" around it and hides the
//     // real macrocycle behind them.
//     fn find_rings(bonds: &[Vec<Bond>], conformers: &[&Conformer]) -> (Vec<Ring>, Vec<Vec<usize>>) {
//         let mut visited: Vec<bool> = (0..bonds.len())
//             .map(|atom_index| Self::excluded_from_rings(bonds, conformers, atom_index))
//             .collect();
//         let mut on_path = vec![false; bonds.len()];
//         let mut rings: Vec<Ring> = vec![];
//         let mut atom_rings = (0..bonds.len()).map(|_| Vec::new()).collect::<Vec<Vec<usize>>>();
//         let mut path: Vec<usize> = vec![];
//         for i in 0..bonds.len() {
//             if !visited[i] {
//                 Self::visit(
//                     bonds,
//                     conformers,
//                     i,
//                     None,
//                     &mut visited,
//                     &mut path,
//                     &mut on_path,
//                     &mut rings,
//                     &mut atom_rings,
//                 );
//             }
//         }
//         (rings, atom_rings)
//     }

//     fn excluded_from_rings(bonds: &[Vec<Bond>], conformers: &[&Conformer], atom_index: usize) -> bool {
//         match bonds[atom_index].first() {
//             None => false,
//             Some(bond) => {
//                 let pointer = bond.atom_1;
//                 conformers[pointer.conformer_index].topology.atoms[pointer.atom_index]
//                     .typical_valence()
//                     .is_none()
//             }
//         }
//     }

//     fn visit(
//         bonds: &[Vec<Bond>],
//         conformers: &[&Conformer],
//         node: usize,
//         parent: Option<usize>,
//         visited: &mut Vec<bool>,
//         path: &mut Vec<usize>,
//         on_path: &mut Vec<bool>,
//         rings: &mut Vec<Ring>,
//         atom_rings: &mut Vec<Vec<usize>>,
//     ) {
//         visited[node] = true;
//         for bond in &bonds[node] {
//             let neighbor = bond.atom_2.index();
//             if !visited[neighbor] {
//                 path.push(node);
//                 on_path[node] = true;
//                 Self::visit(
//                     bonds,
//                     conformers,
//                     neighbor,
//                     Some(node),
//                     visited,
//                     path,
//                     on_path,
//                     rings,
//                     atom_rings,
//                 );
//                 path.pop();
//                 on_path[node] = false;
//             } else if parent == Some(neighbor) {
//                 continue;
//             } else if on_path[neighbor] {
//                 let idx = path.iter().position(|&x| x == neighbor).unwrap();
//                 let cycle: Vec<usize> = path[idx..].iter().copied().chain(std::iter::once(node)).collect();
//                 for piece in Self::split_cycle(cycle, bonds) {
//                     let atoms: Vec<AtomPointer> = piece.iter().map(|&n| bonds[n][0].atom_1).collect();
//                     if let Some(ring) = Ring::new(atoms, conformers) {
//                         let ring_index = rings.len();
//                         for atom_pointer in &ring.atom_pointers {
//                             atom_rings[atom_pointer.index()].push(ring_index);
//                         }
//                         rings.push(ring);
//                     };
//                 }
//             }
//         }
//     }

//     fn split_cycle(cycle: Vec<usize>, bonds: &[Vec<Bond>]) -> Vec<Vec<usize>> {
//         let n = cycle.len();
//         let mut chord = None;
//         'search: for i in 0..n {
//             let prev = cycle[(i + n - 1) % n];
//             let next = cycle[(i + 1) % n];
//             for bond in &bonds[cycle[i]] {
//                 let neighbor = bond.atom_2.index();
//                 if neighbor == prev || neighbor == next {
//                     continue; // one of the cycle's own edges, not a chord
//                 }
//                 if let Some(j) = cycle.iter().position(|&c| c == neighbor) {
//                     chord = Some((i, j));
//                     break 'search;
//                 }
//             }
//         }

//         match chord {
//             None => vec![cycle],
//             Some((i, j)) => {
//                 let (lo, hi) = (i.min(j), i.max(j));
//                 let piece_a: Vec<usize> = cycle[lo..=hi].to_vec();
//                 let piece_b: Vec<usize> = cycle[hi..].iter().chain(cycle[..=lo].iter()).copied().collect();
//                 let mut pieces = Self::split_cycle(piece_a, bonds);
//                 pieces.extend(Self::split_cycle(piece_b, bonds));
//                 pieces
//             }
//         }
//     }
// }
