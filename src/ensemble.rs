use crate::{
    atom::Atom, atom_pointer::AtomPointer, conformer::Conformer, geometry::Point3, molecule::Molecule,
    topology::Topology,
};

#[derive(Debug)]
pub struct Ensemble {
    pub id: usize,
    pub topology: Topology,
    pub conformers: Vec<Conformer>,
}

impl Ensemble {
    pub fn new(id: usize, topology: Topology, conformers: Vec<Conformer>) -> Self {
        Self {
            id,
            topology,
            conformers,
        }
    }

    // TODO: Handle missing conformer
    pub fn molecule(&self, conformer_index: usize) -> Molecule<'_> {
        Molecule {
            ensemble_id: self.id,
            conformer_id: conformer_index,
            topology: &self.topology,
            conformer: &self.conformers[conformer_index],
        }
    }
}

pub fn resolve<'a>(ensembles: &'a [Ensemble], pointer: &AtomPointer) -> (&'a Atom, &'a Point3) {
    let ensemble = ensembles
        .iter()
        .find(|ensemble| ensemble.id == pointer.ensemble_id)
        .expect("AtomPointer references an ensemble that isn't in this list");
    (
        &ensemble.topology.atoms[pointer.index],
        &ensemble.conformers[pointer.conformer_id].positions[pointer.index],
    )
}
