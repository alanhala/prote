use crate::{bond_order::BondOrder, selection::Selection, spatial_index::IndexedStructure};

pub struct TemplateBond {
    pub atom_id_1: String,
    pub atom_id_2: String,
    pub order: BondOrder,
}
#[derive(Debug)]
pub struct Bond {
    pub atom_1: usize,
    pub atom_2: usize,
}

impl Bond {
    pub fn new(atom_1: usize, atom_2: usize) -> Self {
        Self { atom_1, atom_2 }
    }

    // pub fn bonds(indexed: &IndexedStructure, atoms: &Selection) -> Vec<Bond> {
    //     indexed
    //         // PMC3557220, eq. 1: 0.8 < d <= r_1 + r_2 + 0.4
    //         .pairs_within(atoms, atoms, |atom| atom.element.covalent_radius() + 0.2)
    //         .filter(|&(_, _, distance)| distance > 0.8)
    //         .map(|(atom_1, atom_2, distance)| Bond::new(atom_1, atom_2, distance))
    //         .collect()
    // }
}
