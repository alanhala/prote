use crate::{atom_pointer::AtomPointer, conformer::Conformer, proximity::Proximity, spatial_index::SpatialIndex};

#[derive(Debug)]
pub struct ProximityGraph<'a> {
    pub conformers: Vec<&'a Conformer<'a>>,
    pub edges: Vec<Vec<(AtomPointer, AtomPointer, f64)>>,
}

impl<'a> ProximityGraph<'a> {
    pub fn new(
        search: &'a Conformer,
        spatial_index: &SpatialIndex,
        query: &'a Conformer,
        proximity: Proximity,
    ) -> Self {
        let same_conformer = search == query;
        let query_conformer_index = if same_conformer { 0 } else { 1 };
        let conformers = if same_conformer {
            vec![search]
        } else {
            vec![search, query]
        };
        let query_offset = if same_conformer { 0 } else { search.topology.atoms.len() };
        let atoms = query_offset + query.topology.atoms.len();
        let mut edges = (0..atoms)
            .map(|_| Vec::new())
            .collect::<Vec<Vec<(AtomPointer, AtomPointer, f64)>>>();
        for (query_atom_index, query_atom) in query.topology.atoms.iter().enumerate() {
            let candidates = spatial_index.candidates_within(
                &search.frame.positions,
                &query.frame.positions[query_atom_index],
                proximity.search_radius(query_atom),
            );
            for (search_atom_index, distance) in candidates {
                let query_atom_pointer = AtomPointer::new(query_conformer_index, query_offset, query_atom_index);
                let search_atom = &search.topology.atoms[search_atom_index];
                let search_atom_pointer = AtomPointer::new(0, 0, search_atom_index);
                let different_atoms = search_atom_pointer != query_atom_pointer;
                let matches = Proximity::classify(search_atom, query_atom, distance) == Some(proximity);
                if different_atoms && matches {
                    if !edges[search_atom_pointer.index()]
                        .iter()
                        .any(|(_, other, _)| *other == query_atom_pointer)
                    {
                        edges[search_atom_pointer.index()].push((search_atom_pointer, query_atom_pointer, distance));
                    }
                    if !edges[query_atom_pointer.index()]
                        .iter()
                        .any(|(_, other, _)| *other == search_atom_pointer)
                    {
                        edges[query_atom_pointer.index()].push((query_atom_pointer, search_atom_pointer, distance));
                    }
                }
            }
        }
        Self { conformers, edges }
    }
}
