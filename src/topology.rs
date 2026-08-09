use crate::{atom::Atom, component::Component};

#[derive(Debug)]
pub struct Topology {
    pub id: usize,
    pub atoms: Vec<Atom>,
    pub components: Vec<Component>,
}

impl PartialEq for Topology {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}

impl Topology {
    pub fn new(id: usize, atoms: Vec<Atom>, components: Vec<Component>) -> Self {
        Self { id, atoms, components }
    }

    pub fn component_index_for(&self, atom_index: usize) -> usize {
        self.components
            .iter()
            .position(|component| component.atom_range.contains(&atom_index))
            .expect("every atom belongs to a component")
    }
}
