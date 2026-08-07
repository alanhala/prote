use crate::atom::Atom;

#[derive(PartialEq)]
pub enum Proximity {
    Bonded,
    InContact,
}

impl Proximity {
    pub fn classify(atom_1: &Atom, atom_2: &Atom, distance: f64) -> Option<Self> {
        if distance <= atom_1.covalent_radius() + atom_2.covalent_radius() * 1.3 {
            Some(Self::Bonded)
        } else if distance <= atom_1.van_der_waals_radius() + atom_2.van_der_waals_radius() {
            Some(Self::InContact)
        } else {
            None
        }
    }

    pub fn search_radius(&self, atom: &Atom) -> f64 {
        match self {
            Self::Bonded => atom.covalent_radius() + 3.0,
            Self::InContact => atom.van_der_waals_radius() + 3.0,
        }
    }
}
