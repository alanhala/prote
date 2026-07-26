use crate::{conformer::Conformer, molecule::Molecule, topology::Topology};

#[derive(Debug)]
pub struct Ensemble {
    pub topology: Topology,
    pub conformers: Vec<Conformer>,
}

impl Ensemble {
    pub fn new(topology: Topology, conformers: Vec<Conformer>) -> Self {
        Self { topology, conformers }
    }

    // TODO: Handle missing conformer
    pub fn molecule<'a>(&'a self, conformer: usize) -> Molecule<'a> {
        Molecule::new(&self.topology, &self.conformers[conformer])
    }
}
