use crate::{
    atom_pointer::AtomPointer,
    bond::Bond,
    bond_angle::BondAngle,
    bond_order::{BondOrder, BondOrderSource},
    conformer::Conformer,
    proximity::Proximity,
    spatial_index::SpatialIndex,
};

#[derive(Debug)]
pub struct BondGraph<'a> {
    pub bonds: Vec<Vec<Bond<'a>>>,
    pub angles: Vec<Vec<BondAngle<'a>>>,
}

impl<'a> BondGraph<'a> {
    pub fn new(search: &'a Conformer, spatial_index: &SpatialIndex, query: &'a Conformer) -> Self {
        let mut bonds: Vec<Vec<Bond>> = (0..search.topology.atoms.len()).map(|_| Vec::new()).collect();
        for (query_atom_index, query_atom) in query.topology.atoms.iter().enumerate() {
            let position = &query.frame.positions[query_atom_index];
            let search_radius = query_atom.covalent_radius() + 3.0;
            let hits = spatial_index.candidates_within(
                &search.frame.positions,
                position,
                search_radius,
                |search_atom_index, distance| {
                    let same_atom = search == query && query_atom_index == search_atom_index;
                    let search_atom = &search.topology.atoms[search_atom_index];
                    !same_atom && Proximity::classify(query_atom, search_atom, distance) == Some(Proximity::Bonded)
                },
            );
            for (neighbor_atom_index, _) in hits {
                let search_pointer = AtomPointer::new(search, neighbor_atom_index);
                let query_pointer = AtomPointer::new(query, query_atom_index);
                bonds[neighbor_atom_index].push(Bond::new(search_pointer, query_pointer));
            }
        }

        let mut angles: Vec<Vec<BondAngle>> = (0..search.topology.atoms.len()).map(|_| Vec::new()).collect();
        for (atom_index, atom_bonds) in bonds.iter().enumerate() {
            match atom_bonds.len() {
                0 | 1 => {}
                bonds_len => {
                    let vertex_pointer = AtomPointer::new(search, atom_index);
                    let vertex_position = &search.frame.positions[atom_index];
                    for i in 0..(bonds_len - 1) {
                        for j in (i + 1)..bonds_len {
                            let atom_pointer_1 = &atom_bonds[i].atom_2;
                            let atom_pointer_2 = &atom_bonds[j].atom_2;
                            let atom_1_position = &query.frame.positions[atom_pointer_1.index];
                            let atom_2_position = &query.frame.positions[atom_pointer_2.index];
                            angles[atom_index].push(BondAngle::new(
                                (*atom_pointer_1, atom_1_position),
                                (vertex_pointer, vertex_position),
                                (*atom_pointer_2, atom_2_position),
                            ))
                        }
                    }
                }
            }
        }
        Self { bonds, angles }
    }

    pub fn assign_orders(&mut self, sources: &[&dyn BondOrderSource]) {
        for bonds in self.bonds.iter_mut() {
            for bond in bonds.iter_mut() {
                bond.order = Self::resolve_bond_order(bond, sources);
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

    fn resolve_bond_order(bond: &Bond, sources: &[&dyn BondOrderSource]) -> Option<BondOrder> {
        sources.iter().find_map(|source| source.bond_order(bond))
    }
}
