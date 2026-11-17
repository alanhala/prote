use crate::{
    atom_properties::{self, AtomProperties, Hybridization},
    bond::BondKind,
    element::{self, Element},
    mask::Mask,
    selection::Selection,
    structure::Structure,
};

pub struct AtomTypes {
    pub hydrophobic: Selection,
    pub hydrogen_bond_donor: Selection,
    pub hydrogen_bond_acceptor: Selection,
    pub metal: Selection,
    pub metal_binder: Selection,
}

impl AtomTypes {
    pub fn new(structure: &Structure, atom_properties: &AtomProperties) -> Self {
        AtomTypes {
            hydrophobic: Mask::from_bools((0..structure.atoms.len()).map(|atom| is_hydrophobic(structure, atom)).collect()),
            hydrogen_bond_donor: Mask::from_bools((0..structure.atoms.len()).map(|atom| is_hydrogen_bond_donor(structure, atom)).collect()),
            hydrogen_bond_acceptor: Mask::from_bools(
                (0..structure.atoms.len()).map(|atom| is_hydrogen_bond_acceptor(structure, atom_properties, atom)).collect(),
            ),
            metal: Selection::by_atom(structure, |atom| atom.element.is_metal()),
            metal_binder: Selection::by_atom(structure, |atom| matches!(atom.element, Element::N | Element::O | Element::S)),
        }
    }
}

fn is_hydrophobic(structure: &Structure, atom: usize) -> bool {
    let mut partners = structure.covalent_partners(atom).map(|partner| partner.element);
    match structure.atoms[atom].element {
        Element::C => partners.all(|element| !matches!(element, Element::N | Element::O | Element::F)),
        Element::S => {
            // A sulfur holding a metal is not hydrophobic, and the metal is not among its covalent partners.
            let holds_metal = structure
                .bonding
                .bonds
                .neighbors(atom)
                .iter()
                .any(|&(_, edge)| matches!(structure.bonding.bonds.edges()[edge].kind, BondKind::Coordination));
            partners.all(|element| matches!(element, Element::C | Element::S)) && structure.bonding.implicit_hydrogens[atom] == 0 && !holds_metal
        }
        // A free halide ion has no neighbors and is a charged anion, not a hydrophobe.
        Element::Br | Element::I => {
            let elements: Vec<Element> = partners.collect();
            !elements.is_empty() && elements.iter().all(|element| matches!(element, Element::C))
        }
        _ => false,
    }
}

fn is_hydrogen_bond_donor(structure: &Structure, atom: usize) -> bool {
    let bonded_hydrogen = structure.covalent_partners(atom).any(|partner| partner.element == Element::H);
    matches!(structure.atoms[atom].element, Element::O | Element::N) && (structure.bonding.implicit_hydrogens[atom] > 0 || bonded_hydrogen)
}

fn is_hydrogen_bond_acceptor(structure: &Structure, atom_properties: &AtomProperties, atom: usize) -> bool {
    match structure.atoms[atom].element {
        Element::O => atom_properties.lone_pairs[atom].is_some_and(|lone_pairs| lone_pairs > 0),
        Element::N => {
            let shares_pair = matches!(atom_properties.hybridization[atom], Some(Hybridization::Sp2)) && structure.bonding.atoms_bonded_sum(atom) == 3;
            atom_properties.lone_pairs[atom].is_some_and(|lone_pairs| lone_pairs > 0) && !shares_pair
        }
        _ => false,
    }
}
