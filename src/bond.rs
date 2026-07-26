#[derive(Debug)]
pub struct Bond {
    pub atom_1: usize,
    pub atom_2: usize,
    // order: BondOrder
}

impl Bond {
    pub fn new(atom_1: usize, atom_2: usize) -> Self {
        Self { atom_1, atom_2 }
    }
}
