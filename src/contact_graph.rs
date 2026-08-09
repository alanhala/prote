use crate::{
    atom_pointer::AtomPointer, conformer::Conformer, contact::Contact, proximity::Proximity,
    spatial_index::SpatialIndex,
};

#[derive(Debug)]
pub struct ContactGraph<'a> {
    contacts: Vec<Vec<Contact<'a>>>,
}

impl<'a> ContactGraph<'a> {
    pub fn new(search: &'a Conformer, spatial_index: &SpatialIndex, query: &'a Conformer) -> Self {
        let mut contacts: Vec<Vec<Contact>> = (0..search.topology.atoms.len()).map(|_| Vec::new()).collect();
        for (query_atom_index, query_atom) in query.topology.atoms.iter().enumerate() {
            let position = &query.frame.positions[query_atom_index];
            let hits = spatial_index.candidates_within(
                &search.frame.positions,
                position,
                Proximity::search_radius(&Proximity::InContact, query_atom),
                |search_atom_index, distance| {
                    let same_atom = search == query && query_atom_index == search_atom_index;
                    let search_atom = &search.topology.atoms[search_atom_index];
                    !same_atom && Proximity::classify(query_atom, search_atom, distance) == Some(Proximity::InContact)
                },
            );
            for (neighbor_atom_index, distance) in hits {
                let search_pointer = AtomPointer::new(search, neighbor_atom_index);
                let query_pointer = AtomPointer::new(query, query_atom_index);
                contacts[neighbor_atom_index].push(Contact::new(search_pointer, query_pointer, distance));
            }
        }
        Self { contacts }
    }

    pub fn new_intramolecular(conformer: &'a Conformer, spatial_index: &SpatialIndex) -> Self {
        let topology = conformer.topology;
        let atoms = &topology.atoms;
        let positions = &conformer.frame.positions;
        let mut contacts: Vec<Vec<Contact>> = (0..atoms.len()).map(|_| Vec::new()).collect();
        for (query_atom_index, query_atom) in atoms.iter().enumerate() {
            let hits = spatial_index.candidates_within(
                positions,
                &positions[query_atom_index],
                Proximity::search_radius(&Proximity::InContact, query_atom),
                |search_atom_index, distance| {
                    let same_atom = query_atom_index == search_atom_index;
                    let search_atom = &conformer.topology.atoms[search_atom_index];
                    !same_atom && Proximity::classify(query_atom, search_atom, distance) == Some(Proximity::InContact)
                },
            );
            for (neighbor_atom_index, distance) in hits {
                let search_pointer = AtomPointer::new(conformer, neighbor_atom_index);
                let query_pointer = AtomPointer::new(conformer, query_atom_index);
                contacts[neighbor_atom_index].push(Contact::new(search_pointer, query_pointer, distance));
            }
        }
        Self { contacts }
    }
}
