use std::ops::Range;

#[derive(Debug, Clone)]
pub struct Component {
    pub name: String,
    pub atom_range: Range<usize>,
}

impl Component {
    pub fn new(name: String, atom_range: Range<usize>) -> Self {
        Self { name, atom_range }
    }
}
