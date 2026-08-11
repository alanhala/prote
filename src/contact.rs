use crate::atom_pointer::AtomPointer;

#[derive(Debug)]
pub struct Contact {
    pub atom_1: AtomPointer,
    pub atom_2: AtomPointer,
    pub distance: f64,
}

impl Contact {
    pub fn new(atom_1: AtomPointer, atom_2: AtomPointer, distance: f64) -> Self {
        Self {
            atom_1,
            atom_2,
            distance,
        }
    }
}
