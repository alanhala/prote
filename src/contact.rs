#[derive(Debug)]
pub struct Contact {
    pub atom_1: usize,
    pub atom_2: usize,
    pub distance: f64,
}

impl Contact {
    pub fn new(atom_1: usize, atom_2: usize, distance: f64) -> Self {
        Self {
            atom_1,
            atom_2,
            distance,
        }
    }
}
