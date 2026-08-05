use crate::topology::Topology;

#[derive(Debug, Clone, Copy)]
pub enum BondOrder {
    Single,
    Double,
    Triple,
    Quadruple,
    Pi,
    Aromatic,
}

pub trait BondOrderSource {
    fn bond_order(&self, topology: &Topology, atom_1: usize, atom_2: usize) -> Option<BondOrder>;
}

impl BondOrder {
    pub fn as_u8(self) -> u8 {
        match self {
            BondOrder::Single => 1,
            BondOrder::Double => 2,
            BondOrder::Triple => 3,
            BondOrder::Quadruple => 4,
            BondOrder::Pi | BondOrder::Aromatic => 1,
        }
    }

    pub fn from_u8(value: u8) -> Option<BondOrder> {
        match value {
            1 => Some(BondOrder::Single),
            2 => Some(BondOrder::Double),
            3 => Some(BondOrder::Triple),
            4 => Some(BondOrder::Quadruple),
            _ => None,
        }
    }
}
