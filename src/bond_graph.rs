use crate::{
    atom::Atom,
    atom_pointer::AtomPointer,
    bond::Bond,
    bond_angle::BondAngle,
    bond_order::{BondOrder, BondOrderSource},
    ensemble::{resolve, Ensemble},
    geometry::Point3,
    molecule::Molecule,
    spatial_index::SpatialIndex,
    topology::Topology,
};

#[derive(Debug)]
pub struct BondGraph {
    pub bonds: Vec<Vec<Bond>>,
    pub angles: Vec<Vec<BondAngle>>,
}

impl BondGraph {
    pub fn new(search: &Molecule, spatial_index: &SpatialIndex, query: &Molecule) -> Self {
        let mut bonds: Vec<Vec<Bond>> = (0..search.topology.atoms.len()).map(|_| Vec::new()).collect();
        for (query_index, query_atom) in query.topology.atoms.iter().enumerate() {
            let position = &query.conformer.positions[query_index];
            let search_radius = query_atom.covalent_radius() + 3.0;
            let hits = spatial_index.candidates_within(
                &search.conformer.positions,
                position,
                search_radius,
                |neighbor, distance| {
                    let same_atom = search.ensemble_id == query.ensemble_id
                        && search.conformer_id == query.conformer_id
                        && neighbor == query_index;
                    !same_atom
                        && distance
                            <= query_atom.covalent_radius() + search.topology.atoms[neighbor].covalent_radius() * 1.3
                },
            );
            for (neighbor, _) in hits {
                let search_pointer = AtomPointer::new(search.ensemble_id, search.conformer_id, neighbor);
                let query_pointer = AtomPointer::new(query.ensemble_id, query.conformer_id, query_index);
                bonds[neighbor].push(Bond::new(search_pointer, query_pointer));
            }
        }

        let mut angles: Vec<Vec<BondAngle>> = (0..search.topology.atoms.len()).map(|_| Vec::new()).collect();
        for (atom_index, atom_bonds) in bonds.iter().enumerate() {
            match atom_bonds.len() {
                0 | 1 => {}
                bonds_len => {
                    let vertex_pointer = AtomPointer::new(search.ensemble_id, search.conformer_id, atom_index);
                    let vertex_position = &search.conformer.positions[atom_index];
                    for i in 0..(bonds_len - 1) {
                        for j in (i + 1)..bonds_len {
                            let atom_pointer_1 = atom_bonds[i].atom_2;
                            let atom_pointer_2 = atom_bonds[j].atom_2;
                            let atom_1_position = &query.conformer.positions[atom_pointer_1.index];
                            let atom_2_position = &query.conformer.positions[atom_pointer_2.index];
                            angles[atom_index].push(BondAngle::new(
                                (atom_pointer_1, atom_1_position),
                                (vertex_pointer, vertex_position),
                                (atom_pointer_2, atom_2_position),
                            ))
                        }
                    }
                }
            }
        }
        Self { bonds, angles }
    }

    pub fn assign_orders(&mut self, topology: &Topology, sources: &[&dyn BondOrderSource]) {
        for bonds in self.bonds.iter_mut() {
            for bond in bonds.iter_mut() {
                bond.order = Self::resolve_bond_order(topology, &bond.atom_1, &bond.atom_2, sources);
            }
        }
    }

    /// `bonds` is that atom's whole bond list, unresolved entries included —
    /// `filter_map` drops anything still `None` on its own, so there's no need
    /// to exclude the bond currently being resolved by hand.
    fn leftover_capacity(atom: &Atom, bonds: &[Bond]) -> Option<u8> {
        let used: u8 = bonds.iter().filter_map(|bond| bond.order).map(BondOrder::as_u8).sum();
        atom.typical_valence().map(|typical| typical.saturating_sub(used))
    }

    /// A second pass, run after `assign_orders`: for whatever's still `None`,
    /// take each side's leftover valence capacity (from whichever of its other
    /// bonds the dictionary/rules already resolved) and use the smaller of the
    /// two. The smaller side is trustworthy even when the other side's count is
    /// inflated by something missing from the data (e.g. an unresolved
    /// hydrogen) — it can't suggest a higher order than the atom actually has
    /// room for, so it can never be the wrong, inflated number.
    pub fn assign_remaining_orders(&mut self, ensembles: &[Ensemble]) {
        let mut updates: Vec<(usize, usize, BondOrder)> = Vec::new();
        for (search_index, bonds) in self.bonds.iter().enumerate() {
            for (bond_index, bond) in bonds.iter().enumerate() {
                if bond.order.is_some() {
                    continue;
                }
                let (atom_1, _) = resolve(ensembles, &bond.atom_1);
                let (atom_2, _) = resolve(ensembles, &bond.atom_2);
                let leftover_1 = Self::leftover_capacity(atom_1, &self.bonds[bond.atom_1.index]);
                let leftover_2 = Self::leftover_capacity(atom_2, &self.bonds[bond.atom_2.index]);
                let order = leftover_1
                    .zip(leftover_2)
                    .and_then(|(a, b)| BondOrder::from_u8(a.min(b)));
                if let Some(order) = order {
                    updates.push((search_index, bond_index, order));
                }
            }
        }
        for (search_index, bond_index, order) in updates {
            self.bonds[search_index][bond_index].order = Some(order);
        }
    }

    fn resolve_bond_order(
        topology: &Topology,
        atom_1: &AtomPointer,
        atom_2: &AtomPointer,
        sources: &[&dyn BondOrderSource],
    ) -> Option<BondOrder> {
        sources
            .iter()
            .find_map(|source| source.bond_order(topology, atom_1.index, atom_2.index))
    }
}
