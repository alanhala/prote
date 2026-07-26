use crate::{conformer::Conformer, topology::Topology};

#[derive(Debug)]
pub struct Ensemble {
    pub topology: Topology,
    pub conformers: Vec<Conformer>,
}

impl Ensemble {
    pub fn new(topology: Topology, conformers: Vec<Conformer>) -> Self {
        Self { topology, conformers }
    }
}
