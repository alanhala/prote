use crate::{
    atom_pointer::AtomPointer, bond::Bond, bond_angle::BondAngle, molecule::Molecule, spatial_index::SpatialIndex,
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
}
