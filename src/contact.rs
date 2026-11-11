use crate::{element::Element, selection::Selection, spatial_index::SpatialIndex, structure::Structure};

#[derive(Debug)]
pub struct Contact {
    pub atom_1: usize,
    pub atom_2: usize,
    pub distance: f64,
    pub kind: Option<ContactType>,
}

#[derive(Debug)]
pub enum ContactType {
    Hydrophobic,
}

impl Contact {
    pub fn new(atom_1: usize, atom_2: usize, distance: f64, kind: Option<ContactType>) -> Self {
        Self { atom_1, atom_2, distance, kind }
    }

    pub fn contacts(structure: &Structure, index: &SpatialIndex, query: &Selection, search: &Selection) -> Vec<Contact> {
        index
            .pairs_within(query, search, |atom_id| structure.atoms[atom_id].element.van_der_waals_radius().unwrap())
            .filter(|&(atom_1, atom_2, _)| !structure.bonding.bonds.within_steps(atom_1, atom_2, 3))
            .map(|(atom_1_id, atom_2_id, distance)| {
                let kind =
                    if is_hydrophobic(structure, atom_1_id) && is_hydrophobic(structure, atom_2_id) { Some(ContactType::Hydrophobic) } else { None };
                Contact::new(atom_1_id, atom_2_id, distance, kind)
            })
            .collect()
    }
}

fn is_hydrophobic(structure: &Structure, atom: usize) -> bool {
    let bonded_elements = || structure.bonding.bonds.neighbors(atom).iter().map(|&(neighbor, _)| structure.atoms[neighbor].element);
    match structure.atoms[atom].element {
        Element::C => bonded_elements().all(|element| !matches!(element, Element::N | Element::O)),
        Element::S => bonded_elements().all(|element| matches!(element, Element::C | Element::S)),
        _ => false,
    }
}
