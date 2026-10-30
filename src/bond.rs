use crate::{bond_order::BondOrder, frame::Frame, graph::Edge};

pub struct ResidueTemplateBond {
    pub atom_name_1: String,
    pub atom_name_2: String,
    pub order: BondOrder,
    pub aromatic: bool,
}

#[derive(Debug)]
pub struct Bond {
    pub atom_1: usize,
    pub atom_2: usize,
    pub order: BondOrder,
    pub aromatic: bool,
}

impl Edge for Bond {
    fn nodes(&self) -> (usize, usize) {
        (self.atom_1, self.atom_2)
    }
}

impl Bond {
    pub fn length(&self, frame: &Frame) -> f64 {
        frame.positions[self.atom_1].distance(&frame.positions[self.atom_2])
    }
}
