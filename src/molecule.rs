use crate::{conformer::Conformer, topology::Topology};

#[derive(Debug)]
pub struct Molecule<'a> {
    pub topology: &'a Topology,
    pub conformer: &'a Conformer,
}

impl<'a> Molecule<'a> {
    pub fn new(topology: &'a Topology, conformer: &'a Conformer) -> Self {
        Self { topology, conformer }
    }
}
