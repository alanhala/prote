use crate::{atom_pointer::AtomPointer, bond_order::BondOrder};

#[derive(Debug)]
pub struct Bond<'a> {
    pub atom_1: AtomPointer<'a>,
    pub atom_2: AtomPointer<'a>,
    pub order: Option<BondOrder>,
}

impl<'a> Bond<'a> {
    pub fn new(atom_1: AtomPointer<'a>, atom_2: AtomPointer<'a>) -> Self {
        Self {
            atom_1,
            atom_2,
            order: None,
        }
    }
}
