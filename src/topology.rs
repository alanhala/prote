use crate::atom::Atom;
use crate::residue::Residue;

#[derive(Debug)]
pub struct Topology {
    pub name: String,
    pub atoms: Vec<Atom>,
    pub residues: Vec<Residue>,
}

impl Topology {
    pub fn new(name: String, atoms: Vec<Atom>, residues: Vec<Residue>) -> Self {
        Self { name, atoms, residues }
    }

    pub fn residue_index_for(&self, atom_index: usize) -> usize {
        self.residues
            .iter()
            .position(|residue| residue.atom_range.contains(&atom_index))
            .expect("every atom belongs to a residue")
    }
}
