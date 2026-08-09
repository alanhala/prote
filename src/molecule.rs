use crate::{conformer::Conformer, frame::Frame, topology::Topology};

#[derive(Debug)]
pub struct Molecule {
    pub id: usize,
    pub name: String,
    pub topology: Topology,
    pub frames: Vec<Frame>,
}

impl Molecule {
    pub fn new(id: usize, name: String, topology: Topology, frames: Vec<Frame>) -> Self {
        Self {
            id,
            name,
            topology,
            frames,
        }
    }

    pub fn conformer(&self, index: usize) -> Option<Conformer<'_>> {
        let frame = self.frames.get(index)?;
        Some(Conformer::new(&self.topology, frame))
    }
}
