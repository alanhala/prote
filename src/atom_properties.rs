use crate::{
    bond_order::BondOrder,
    bonding::{covalent_partners, Bonding},
    structure::Structure,
};

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Hybridization {
    Sp,
    Sp2,
    Sp3,
}

pub struct AtomProperties {
    pub formal_charges: Vec<Option<i8>>,
    pub lone_pairs: Vec<Option<u8>>,
    pub hybridization: Vec<Option<Hybridization>>,
}

impl AtomProperties {
    pub fn new(structure: &Structure, bonding: &Bonding) -> Self {
        let formal_charges = formal_charges(structure, bonding);
        let lone_pairs: Vec<Option<u8>> = (0..structure.atoms.len()).map(|atom| lone_pair(structure, bonding, atom, &formal_charges)).collect();
        let hybridization = (0..structure.atoms.len()).map(|atom| hybridization(structure, bonding, &lone_pairs, atom)).collect();
        Self { formal_charges, lone_pairs, hybridization }
    }
}

fn formal_charges(structure: &Structure, bonding: &Bonding) -> Vec<Option<i8>> {
    let mut charges = vec![None; structure.atoms.len()];
    for (id, atom) in structure.atoms.iter().enumerate() {
        charges[id] = atom
            .element
            .valences()
            .iter()
            .map(|&valence| bonding.bond_order_sum(id) as i8 - valence as i8)
            .min_by_key(|charge| charge.abs())
            // A metal has no usual number of bonds to count against, so only the file can say its charge.
            .or(atom.stated_charge);
    }
    charges
}

fn hybridization(structure: &Structure, bonding: &Bonding, lone_pairs: &[Option<u8>], atom: usize) -> Option<Hybridization> {
    let lone_pairs = lone_pairs[atom]?;
    let neighbour_has_double_bond = covalent_partners(structure, bonding, atom)
        .filter_map(|partner| partner.atom)
        .any(|neighbour| covalent_partners(structure, bonding, neighbour).any(|partner| partner.order != BondOrder::Single || partner.aromatic));
    match lone_pairs + bonding.atoms_bonded_sum(atom) as u8 {
        4 if lone_pairs > 0 && neighbour_has_double_bond => Some(Hybridization::Sp2),
        4 => Some(Hybridization::Sp3),
        3 => Some(Hybridization::Sp2),
        2 => Some(Hybridization::Sp),
        _ => None,
    }
}

fn lone_pair(structure: &Structure, bonding: &Bonding, atom: usize, charges: &[Option<i8>]) -> Option<u8> {
    let charge = charges[atom]?;
    let valence_electrons = structure.atoms[atom].element.valence_electrons() as i8;
    // More bonds than the element has electrons for means the bonding is off, so there is no honest count.
    let free_electrons = u8::try_from(valence_electrons - bonding.bond_order_sum(atom) as i8 - charge).ok()?;
    Some(free_electrons / 2)
}
