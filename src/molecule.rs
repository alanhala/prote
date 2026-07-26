use crate::{conformer::Conformer, topology::Topology};

#[derive(Debug)]
pub struct Molecule {
    pub topology: Topology,
    pub conformer: Conformer,
}

impl Molecule {
    pub fn new(topology: Topology, conformer: Conformer) -> Self {
        Self { topology, conformer }
    }
}
