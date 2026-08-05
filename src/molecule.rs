use crate::{conformer::Conformer, topology::Topology};

#[derive(Debug)]
pub struct Molecule<'a> {
    pub ensemble_id: usize,
    pub conformer_id: usize,
    pub topology: &'a Topology,
    pub conformer: &'a Conformer,
}

impl<'a> Molecule<'a> {
    pub fn new(ensemble_id: usize, conformer_id: usize, topology: &'a Topology, conformer: &'a Conformer) -> Self {
        Self {
            ensemble_id,
            conformer_id,
            topology,
            conformer,
        }
    }
}
