use crate::{selection::Selection, spatial_index::IndexedStructure};

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

    pub fn contacts(indexed: &IndexedStructure, query: &Selection, search: &Selection) -> Vec<Contact> {
        indexed
            .pairs_within(query, search, |atom| atom.element.van_der_waals_radius())
            .filter(|&(atom_1, atom_2, distance)| {
                distance
                    > indexed.conformer().structure().atoms[atom_1].element.covalent_radius()
                        + indexed.conformer().structure().atoms[atom_2].element.covalent_radius()
            })
            .map(|(atom_1, atom_2, distance)| Contact::new(atom_1, atom_2, distance))
            .collect()
    }
}
