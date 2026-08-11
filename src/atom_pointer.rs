// TODO: This is not longer an AtomPointer and instead a graph pointer
#[derive(Debug, PartialEq, Eq, Copy, Clone)]
pub struct AtomPointer {
    pub conformer_index: usize,
    pub conformer_offset: usize,
    pub atom_index: usize,
}

impl AtomPointer {
    pub fn new(conformer_index: usize, conformer_offset: usize, atom_index: usize) -> Self {
        Self {
            conformer_index,
            conformer_offset,
            atom_index,
        }
    }

    pub fn index(&self) -> usize {
        self.conformer_offset + self.atom_index
    }
}
