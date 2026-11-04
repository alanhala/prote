use crate::{bond_order::BondOrder, element::Element, frame::Frame, graph::Edge};

#[derive(Debug)]
pub struct ResidueTemplateBond {
    pub atom_name_1: String,
    pub element_1: Element,
    pub atom_name_2: String,
    pub element_2: Element,
    pub order: BondOrder,
    pub aromatic: bool,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum LinkType {
    Covalent,
    Disulfide,
    MetalCoordination,
    HydrogenBond,
}

#[derive(Debug)]
pub struct Link {
    pub atom_1: usize,
    pub atom_2: usize,
    pub link_type: LinkType,
    pub order: Option<BondOrder>,
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
