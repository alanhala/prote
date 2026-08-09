use crate::atom_pointer::AtomPointer;

#[derive(Debug)]
pub struct Contact<'a> {
    pub atom_1: AtomPointer<'a>,
    pub atom_2: AtomPointer<'a>,
    pub distance: f64,
}

impl<'a> Contact<'a> {
    pub fn new(atom_1: AtomPointer<'a>, atom_2: AtomPointer<'a>, distance: f64) -> Self {
        Self {
            atom_1,
            atom_2,
            distance,
        }
    }
}
