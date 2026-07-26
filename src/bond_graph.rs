use crate::{atom::Atom, bond::Bond, position::Position, spatial_index::SpatialIndex};

pub struct BondGraph {
    bonds: Vec<Vec<Bond>>,
}

impl BondGraph {
    pub fn new(atoms: &[Atom], positions: &[Position]) -> Self {
        let spatial_index = SpatialIndex::new(positions);
        let mut bonds: Vec<Vec<Bond>> = (0..atoms.len()).map(|_| Vec::new()).collect();
        for (i, atom) in atoms.iter().enumerate() {
            let position = &positions[i];
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
            for neighbor in spatial_index.search(positions, min, max) {
                let in_contact = positions[i].distance(&positions[neighbor])
                    <= atom.covalent_radius() + atoms[neighbor].covalent_radius();
                if neighbor > i && in_contact {
                    bonds[i].push(Bond::new(i, neighbor));
                    bonds[neighbor].push(Bond::new(neighbor, i));
                }
            }
        }
        Self { bonds }
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.iter().map(Vec::len).sum()
    }
}
