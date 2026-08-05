#[derive(Debug, Clone, Copy)]
pub struct AtomPointer {
    pub ensemble_id: usize,
    pub conformer_id: usize,
    pub index: usize,
}

impl AtomPointer {
    pub fn new(ensemble_id: usize, conformer_id: usize, index: usize) -> Self {
        Self {
            ensemble_id,
            conformer_id,
            index,
        }
    }
}
