use std::fmt;

use crate::conformer::Conformer;

#[derive(Copy, Clone)]
pub struct AtomPointer<'a> {
    pub conformer: &'a Conformer<'a>,
    pub index: usize,
}

impl<'a> AtomPointer<'a> {
    pub fn new(conformer: &'a Conformer, index: usize) -> Self {
        Self { conformer, index }
    }
}

impl fmt::Debug for AtomPointer<'_> {
    // Not derived: `conformer` is a full `&Conformer` (whole topology + frame), and deriving
    // would recursively print every atom in the molecule for a single pointer. Print just the
    // one atom this pointer actually refers to.
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let atom = &self.conformer.topology.atoms[self.index];
        f.debug_struct("AtomPointer")
            .field("index", &self.index)
            .field("atom", &format!("{}{}", atom.symbol(), atom.name))
            .finish()
    }
}
