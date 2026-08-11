use crate::{frame::Frame, topology::Topology};

#[derive(Debug)]
pub struct Conformer<'a> {
    pub topology: &'a Topology,
    pub frame: &'a Frame,
}

impl<'a> PartialEq for Conformer<'a> {
    fn eq(&self, other: &Self) -> bool {
        self.topology == other.topology && self.frame == other.frame
    }
}

impl<'a> Eq for Conformer<'a> {}

impl<'a> Conformer<'a> {
    pub fn new(topology: &'a Topology, frame: &'a Frame) -> Self {
        Conformer { topology, frame }
    }
}
