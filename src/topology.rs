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
}
