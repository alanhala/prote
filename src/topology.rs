use crate::atom::Atom;
use crate::residue::Residue;

#[derive(Debug)]
pub struct Topology {
    name: String,
    atoms: Vec<Atom>,
    residues: Vec<Residue>,
}

impl Topology {
    pub fn new(name: String, atoms: Vec<Atom>, residues: Vec<Residue>) -> Self {
        Self { name, atoms, residues }
    }

    pub fn atoms(&self) -> &[Atom] {
        &self.atoms
    }

    pub fn residues(&self) -> &[Residue] {
        &self.residues
    }
}
