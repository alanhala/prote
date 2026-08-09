use crate::conformer::Conformer;

#[derive(Debug, Copy, Clone)]
pub struct AtomPointer<'a> {
    pub conformer: &'a Conformer<'a>,
    pub index: usize,
}

impl<'a> AtomPointer<'a> {
    pub fn new(conformer: &'a Conformer, index: usize) -> Self {
        Self { conformer, index }
    }
}
