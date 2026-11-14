use crate::{bond_order::BondOrder, element::Element, frame::Frame, geometry::Point3, graph::Edge};

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
    pub kind: BondKind,
}

#[derive(Debug)]
pub enum BondKind {
    Covalent { order: BondOrder, aromatic: bool, source: BondSource },
    Coordination,
}

#[derive(Debug, Clone, Copy)]
pub enum BondSource {
    Template,
    Backbone,
    Disulfide,
    CovalentLink,
}

impl Edge for Bond {
    fn nodes(&self) -> (usize, usize) {
        (self.atom_1, self.atom_2)
    }
}

impl Bond {
    pub fn length(&self, positions: &[Point3]) -> f64 {
        positions[self.atom_1].distance(&positions[self.atom_2])
    }
}
