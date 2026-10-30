use crate::{selection::Selection, spatial_index::SpatialIndex, structure::Structure};

#[derive(Debug)]
pub struct Contact {
    pub atom_1: usize,
    pub atom_2: usize,
    pub distance: f64,
}

impl Contact {
    pub fn new(atom_1: usize, atom_2: usize, distance: f64) -> Self {
        Self { atom_1, atom_2, distance }
    }

    pub fn contacts(structure: &Structure, index: &SpatialIndex, query: &Selection, search: &Selection) -> Vec<Contact> {
        index
            .pairs_within(query, search, |atom_id| structure.atoms[atom_id].element.van_der_waals_radius())
            .filter(|&(atom_1, atom_2, _)| !structure.bonds.within_steps(atom_1, atom_2, 3))
            .map(|(atom_1, atom_2, distance)| Contact::new(atom_1, atom_2, distance))
            .collect()
    }
}
