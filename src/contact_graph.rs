use crate::{
    atom::Atom, atom_pointer::AtomPointer, bond_graph::BondGraph, contact::Contact, geometry::Point3,
    molecule::Molecule, proximity::Proximity, spatial_index::SpatialIndex,
};

#[derive(Debug)]
pub struct ContactGraph {
    contacts: Vec<Vec<Contact>>,
}

impl ContactGraph {
    pub fn new(search: &Molecule, spatial_index: &SpatialIndex, query: &Molecule) -> Self {
        let mut contacts: Vec<Vec<Contact>> = (0..search.topology.atoms.len()).map(|_| Vec::new()).collect();
        for (query_index, query_atom) in query.topology.atoms.iter().enumerate() {
            let position = &query.conformer.positions[query_index];
            let hits = spatial_index.candidates_within(
                &search.conformer.positions,
                position,
                Proximity::search_radius(&Proximity::InContact, query_atom),
                |neighbor, distance| {
                    Proximity::classify(&query_atom, &search.topology.atoms[neighbor], distance)
                        == Some(Proximity::InContact)
                },
            );
            for (neighbor, distance) in hits {
                let neighbor_pointer = AtomPointer::new(search.ensemble_id, search.conformer_id, neighbor);
                let query_pointer = AtomPointer::new(query.ensemble_id, query.conformer_id, query_index);
                contacts[neighbor].push(Contact::new(neighbor_pointer, query_pointer, distance));
            }
        }
        Self { contacts }
    }

    pub fn new_intramolecular(molecule: &Molecule, spatial_index: &SpatialIndex) -> Self {
        let topology = molecule.topology;
        let atoms = &topology.atoms;
        let positions = &molecule.conformer.positions;
        let mut contacts: Vec<Vec<Contact>> = (0..atoms.len()).map(|_| Vec::new()).collect();
        for (atom_index, atom) in atoms.iter().enumerate() {
            let residue_index = topology.residue_index_for(atom_index);
            let hits = spatial_index.candidates_within(
                positions,
                &positions[atom_index],
                Proximity::search_radius(&Proximity::InContact, atom),
                |neighbor, distance| {
                    neighbor > atom_index
                        && topology.residue_index_for(neighbor) != residue_index
                        && Proximity::classify(&atom, &molecule.topology.atoms[neighbor], distance)
                            == Some(Proximity::InContact)
                },
            );
            for (neighbor, distance) in hits {
                let atom_pointer = AtomPointer::new(molecule.ensemble_id, molecule.conformer_id, atom_index);
                let neighbor_pointer = AtomPointer::new(molecule.ensemble_id, molecule.conformer_id, neighbor);
                contacts[atom_index].push(Contact::new(atom_pointer, neighbor_pointer, distance));
            }
        }
        Self { contacts }
    }

    pub fn contact_count(&self) -> usize {
        self.contacts.iter().map(Vec::len).sum()
    }

    pub fn contacts(&self) -> impl Iterator<Item = &Contact> {
        self.contacts.iter().flatten()
    }
}
