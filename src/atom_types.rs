use crate::{
    atom_properties::{self, AtomProperties, Hybridization},
    bond::BondKind,
    bonding::{covalent_partners, Bonding},
    element::{self, Element},
    mask::Mask,
    selection::Selection,
    structure::{EntityKind, Structure},
};

pub struct AtomTypes {
    pub hydrophobic: Selection,
    pub hydrogen_bond_donor: Selection,
    pub hydrogen_bond_acceptor: Selection,
    pub metal: Selection,
    pub metal_binder: Selection,
    pub cation: Selection,
    pub anion: Selection,
    pub water: Selection,
    pub halogen_bond_donor: Selection,
    pub halogen_bond_acceptor: Selection,
}

impl AtomTypes {
    pub fn new(structure: &Structure, bonding: &Bonding, atom_properties: &AtomProperties) -> Self {
        AtomTypes {
            hydrophobic: Mask::from_bools((0..structure.atoms.len()).map(|atom| is_hydrophobic(structure, bonding, atom)).collect()),
            hydrogen_bond_donor: Mask::from_bools((0..structure.atoms.len()).map(|atom| is_hydrogen_bond_donor(structure, bonding, atom)).collect()),
            hydrogen_bond_acceptor: Mask::from_bools(
                (0..structure.atoms.len()).map(|atom| is_hydrogen_bond_acceptor(structure, bonding, atom_properties, atom)).collect(),
            ),
            metal: Selection::by_atom(structure, |atom| atom.element.is_metal()),
            metal_binder: Selection::by_atom(structure, |atom| matches!(atom.element, Element::N | Element::O | Element::S)),
            cation: Mask::from_bools((0..structure.atoms.len()).map(|atom| shares_charge(structure, bonding, atom_properties, atom, 1)).collect()),
            anion: Mask::from_bools((0..structure.atoms.len()).map(|atom| shares_charge(structure, bonding, atom_properties, atom, -1)).collect()),
            water: Selection::by_instance(structure, |_, entity| entity.kind == EntityKind::Water),
            halogen_bond_donor: Selection::by_atom(structure, |atom| matches!(atom.element, Element::Cl | Element::Br | Element::I)),
            halogen_bond_acceptor: Mask::from_bools(
                (0..structure.atoms.len()).map(|atom| is_halogen_bond_acceptor(structure, bonding, atom_properties, atom)).collect(),
            ),
        }
    }
}

fn is_halogen_bond_acceptor(structure: &Structure, bonding: &Bonding, atom_properties: &AtomProperties, atom: usize) -> bool {
    match structure.atoms[atom].element {
        Element::S => atom_properties.lone_pairs[atom].is_some_and(|lone_pairs| lone_pairs > 0),
        _ => is_hydrogen_bond_acceptor(structure, bonding, atom_properties, atom),
    }
}

fn shares_charge(structure: &Structure, bonding: &Bonding, atom_properties: &AtomProperties, atom: usize, sign: i8) -> bool {
    let charged = |atom: usize| atom_properties.formal_charges[atom].is_some_and(|charge| charge.signum() == sign);
    // Resonance: the two oxygens of a carboxylate are identical, yet the file draws the double bond, and so the
    // charge, on one of them. A twin on the same centre whose bond could swap with this atom's shares its charge.
    let twin_charged = covalent_partners(structure, bonding, atom).any(|to_centre| {
        to_centre.atom.is_some_and(|centre| {
            covalent_partners(structure, bonding, centre).any(|to_twin| {
                to_twin.atom.is_some_and(|twin| twin != atom && charged(twin))
                    && to_twin.element == structure.atoms[atom].element
                    && to_twin.order != to_centre.order
            })
        })
    });
    charged(atom) || twin_charged
}

fn is_hydrophobic(structure: &Structure, bonding: &Bonding, atom: usize) -> bool {
    let mut partners = covalent_partners(structure, bonding, atom).map(|partner| partner.element);
    match structure.atoms[atom].element {
        Element::C => partners.all(|element| !matches!(element, Element::N | Element::O | Element::F)),
        Element::S => {
            // A sulfur holding a metal is not hydrophobic, and the metal is not among its covalent partners.
            let holds_metal =
                bonding.bonds.neighbors(atom).iter().any(|&(_, edge)| matches!(bonding.bonds.edges()[edge].kind, BondKind::Coordination));
            partners.all(|element| matches!(element, Element::C | Element::S)) && bonding.implicit_hydrogens[atom] == 0 && !holds_metal
        }
        // A free halide ion has no neighbors and is a charged anion, not a hydrophobe.
        Element::Br | Element::I => {
            let elements: Vec<Element> = partners.collect();
            !elements.is_empty() && elements.iter().all(|element| matches!(element, Element::C))
        }
        _ => false,
    }
}

fn is_hydrogen_bond_donor(structure: &Structure, bonding: &Bonding, atom: usize) -> bool {
    let bonded_hydrogen = covalent_partners(structure, bonding, atom).any(|partner| partner.element == Element::H);
    matches!(structure.atoms[atom].element, Element::O | Element::N) && (bonding.implicit_hydrogens[atom] > 0 || bonded_hydrogen)
}

fn is_hydrogen_bond_acceptor(structure: &Structure, bonding: &Bonding, atom_properties: &AtomProperties, atom: usize) -> bool {
    match structure.atoms[atom].element {
        Element::O => atom_properties.lone_pairs[atom].is_some_and(|lone_pairs| lone_pairs > 0),
        Element::N => {
            let shares_pair = matches!(atom_properties.hybridization[atom], Some(Hybridization::Sp2)) && bonding.atoms_bonded_sum(atom) == 3;
            atom_properties.lone_pairs[atom].is_some_and(|lone_pairs| lone_pairs > 0) && !shares_pair
        }
        _ => false,
    }
}
