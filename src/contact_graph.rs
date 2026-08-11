use crate::{conformer::Conformer, contact::Contact, proximity::Proximity, proximity_graph::ProximityGraph, spatial_index::SpatialIndex};

#[derive(Debug)]
pub struct ContactGraph<'a> {
    pub conformers: Vec<&'a Conformer<'a>>,
    contacts: Vec<Vec<Contact>>,
}

impl<'a> ContactGraph<'a> {
    pub fn new(search: &'a Conformer, spatial_index: &SpatialIndex, query: &'a Conformer) -> Self {
        let proximity_graph = ProximityGraph::new(search, spatial_index, query, Proximity::InContact);
        let contacts: Vec<Vec<Contact>> = proximity_graph
            .edges
            .into_iter()
            .map(|edges| {
                edges
                    .into_iter()
                    .map(|(atom_1, atom_2, distance)| Contact::new(atom_1, atom_2, distance))
                    .collect()
            })
            .collect();
        Self {
            conformers: proximity_graph.conformers,
            contacts,
        }
    }
}
