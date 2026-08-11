use crate::{
    atom_pointer::AtomPointer, conformer::Conformer, contact::Contact, proximity::Proximity,
    spatial_index::SpatialIndex,
};

#[derive(Debug)]
pub struct ContactGraph<'a> {
    pub conformers: Vec<&'a Conformer<'a>>,
    contacts: Vec<Vec<Contact>>,
}

impl<'a> ContactGraph<'a> {
    pub fn new(search: &'a Conformer, spatial_index: &SpatialIndex, query: &'a Conformer) -> Self {
        let search_offset = 0;
        let query_offset = if search == query {
            0
        } else {
            search.topology.atoms.len()
        };
        let conformers = vec![search, query];
        let atoms = query_offset + query.topology.atoms.len();
        let mut contacts = (0..atoms).map(|_| Vec::new()).collect::<Vec<Vec<Contact>>>();
        for (query_atom_index, query_atom) in query.topology.atoms.iter().enumerate() {
            let candidates = spatial_index.candidates_within(
                &search.frame.positions,
                &query.frame.positions[query_atom_index],
                Proximity::search_radius(Proximity::InContact, &query_atom),
            );
            for (search_atom_index, distance) in candidates {
                let query_atom_pointer = AtomPointer::new(1, query_offset, query_atom_index);
                let search_atom = &search.topology.atoms[search_atom_index];
                let search_atom_pointer = AtomPointer::new(0, search_offset, search_atom_index);
                let different_atoms = search_atom_pointer != query_atom_pointer;
                let bonded = Proximity::classify(&search_atom, &query_atom, distance) == Some(Proximity::InContact);
                if different_atoms && bonded {
                    if !contacts[search_atom_pointer.index()]
                        .iter()
                        .any(|b| b.atom_2 == query_atom_pointer)
                    {
                        contacts[search_atom_pointer.index()].push(Contact::new(
                            search_atom_pointer,
                            query_atom_pointer,
                            distance,
                        ));
                    }
                    if !contacts[query_atom_pointer.index()]
                        .iter()
                        .any(|b| b.atom_2 == search_atom_pointer)
                    {
                        contacts[query_atom_pointer.index()].push(Contact::new(
                            query_atom_pointer,
                            search_atom_pointer,
                            distance,
                        ));
                    }
                }
            }
        }
        Self { conformers, contacts }
    }
}
