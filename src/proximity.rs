use crate::atom::Atom;

#[derive(Debug, PartialEq, Clone, Copy)]
pub enum Proximity {
    Bonded,
    InContact,
}

impl Proximity {
    pub fn classify(atom_1: &Atom, atom_2: &Atom, distance: f64) -> Option<Self> {
        // Bond-detection bounds from the rule-based bond-order paper (PMC3557220,
        // eq. 1): 0.8 < d_ij < r_i + r_j + 0.4. The lower bound rejects
        // near-coincident atoms that would otherwise register as bonded under a
        // purely upper-bounded check.
        //
        // TODO: distance alone can't distinguish "never bonded" from "should be
        // bonded but this instance's coordinates are disordered" (e.g. a proline
        // ring pucker collapsed by poor local density - confirmed via B-factor on
        // 4HHB Pro58 CG-CD, 0.226 A apart). For components with a CCD entry,
        // cross-check candidate pairs against `chem_comp_bonds` (known chemical
        // connectivity, independent of this structure's geometry) before falling
        // back to pure distance - and for anything geometry-only rejects, consider
        // surfacing it as a distinct "Suspect" case rather than silently dropping
        // it, since a hard `None` here loses the information that something was
        // close enough to be worth a second look.
        if distance > 0.8 && distance <= atom_1.covalent_radius() + atom_2.covalent_radius() + 0.4 {
            Some(Self::Bonded)
        } else if distance <= atom_1.van_der_waals_radius() + atom_2.van_der_waals_radius() {
            Some(Self::InContact)
        } else {
            None
        }
    }

    pub fn search_radius(self, atom: &Atom) -> f64 {
        match self {
            Self::Bonded => atom.covalent_radius() + 3.0,
            Self::InContact => atom.van_der_waals_radius() + 3.0,
        }
    }
}
