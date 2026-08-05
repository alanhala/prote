use crate::{atom_pointer::AtomPointer, bond_order::BondOrder};

#[derive(Debug)]
pub struct Bond {
    pub atom_1: AtomPointer,
    pub atom_2: AtomPointer,
    pub order: Option<BondOrder>,
}

impl Bond {
    pub fn new(atom_1: AtomPointer, atom_2: AtomPointer) -> Self {
        Self {
            atom_1,
            atom_2,
            order: None,
        }
    }
}
