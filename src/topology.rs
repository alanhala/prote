use crate::atom::Atom;
use crate::bond::Bond;
use crate::residue::Residue;

#[derive(Debug)]
pub struct Topology {
    pub name: String,
    pub atoms: Vec<Atom>,
    pub residues: Vec<Residue>,
    pub bonds: Vec<Bond>,
}

impl Topology {
    pub fn new(name: String, atoms: Vec<Atom>, residues: Vec<Residue>, bonds: Vec<Bond>) -> Self {
        Self {
            name,
            atoms,
            residues,
            bonds,
        }
    }
}
