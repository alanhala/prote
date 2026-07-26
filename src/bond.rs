use crate::{atom::Atom, conformer::Conformer};

#[derive(Debug)]
pub struct Bond {
    pub atom_1: usize,
    pub atom_2: usize,
    // order: BondOrder
}

impl Bond {
    pub fn perceive_bonds(atoms: &[Atom], conformer: &Conformer) -> Vec<Self> {
        let spatial_index = conformer.spatial_index();
        let mut bonds: Vec<Bond> = vec![];
        for (i, atom) in atoms.iter().enumerate() {
            let position = &conformer.positions[i];
            let search_radius = atom.covalent_radius() + 3.0; // TODO: use a better one
            let min = [
                position.x - search_radius,
                position.y - search_radius,
                position.z - search_radius,
            ];
            let max = [
                position.x + search_radius,
                position.y + search_radius,
                position.z + search_radius,
            ];
            for neighbor in spatial_index.search(&conformer.positions, min, max) {
                let in_contact = conformer.positions[i].distance(&conformer.positions[neighbor])
                    <= atom.covalent_radius() + atoms[neighbor].covalent_radius();
                if neighbor != i && in_contact {
                    bonds.push(Self {
                        atom_1: i,
                        atom_2: neighbor,
                    })
                }
            }
        }
        bonds
    }
}
