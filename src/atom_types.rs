use crate::{element::Element, mask::Mask, selection::Selection, structure::Structure};

pub struct AtomTypes {
    pub hydrophobic: Selection,
    pub metal: Selection,
    pub metal_binder: Selection,
}

impl AtomTypes {
    pub fn new(structure: &Structure) -> Self {
        AtomTypes {
            hydrophobic: Mask::from_bools((0..structure.atoms.len()).map(|atom| is_hydrophobic(structure, atom)).collect()),
            metal: Selection::by_atom(structure, |atom| atom.element.is_metal()),
            metal_binder: Selection::by_atom(structure, |atom| matches!(atom.element, Element::N | Element::O | Element::S)),
        }
    }
}

fn is_hydrophobic(structure: &Structure, atom: usize) -> bool {
    let present = structure.bonding.bonds.neighbors(atom).iter().map(|&(neighbor, _)| structure.atoms[neighbor].element);
    let missing = structure.bonding.missing_partners[atom].iter().map(|partner| partner.element);
    let mut partners = present.chain(missing);
    match structure.atoms[atom].element {
        Element::C => partners.all(|element| !matches!(element, Element::N | Element::O | Element::F)),
        Element::S => partners.all(|element| matches!(element, Element::C | Element::S)) && structure.bonding.implicit_hydrogens[atom] == 0,
        // A free halide ion has no neighbors and is a charged anion, not a hydrophobe.
        Element::Br | Element::I => {
            let neighbors = structure.bonding.bonds.neighbors(atom);
            !neighbors.is_empty() && neighbors.iter().all(|&(neighbor, _)| matches!(structure.atoms[neighbor].element, Element::C))
        }
        _ => false,
    }
}
