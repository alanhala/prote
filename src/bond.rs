use crate::bond_order::BondOrder;

#[derive(Debug)]
pub struct Bond {
    pub atom_1: usize,
    pub atom_2: usize,
    pub order: Option<BondOrder>,
    pub distance: f64,
}

impl Bond {
    pub fn new(atom_1: usize, atom_2: usize, distance: f64) -> Self {
        Self {
            atom_1,
            atom_2,
            order: None,
            distance,
        }
    }
}
