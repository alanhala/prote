use crate::{conformer::Conformer, frame::Frame, topology::Topology};

#[derive(Debug)]
pub struct Protein {
    pub name: String,
    pub topology: Topology,
    frames: Vec<Frame>,
    residue_numbers: Vec<i64>,
}

impl Protein {
    pub fn new(name: String, topology: Topology, frames: Vec<Frame>, residue_numbers: Vec<i64>) -> Self {
        Self {
            name,
            topology,
            frames,
            residue_numbers,
        }
    }

    pub fn conformer(&self, index: usize) -> Option<Conformer<'_>> {
        let frame = self.frames.get(index)?;
        Some(Conformer::new(&self.topology, frame))
    }
}
