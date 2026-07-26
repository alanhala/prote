use std::ops::Range;

#[derive(Debug)]
pub struct Residue {
    pub name: String,
    pub atom_range: Range<usize>,
    pub is_hetero: bool,
}

impl Residue {
    pub fn new(name: String, atom_range: Range<usize>, is_hetero: bool) -> Self {
        Self {
            name,
            atom_range,
            is_hetero,
        }
    }
}
