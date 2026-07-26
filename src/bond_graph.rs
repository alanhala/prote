use crate::{atom::Atom, bond::Bond, bond_angle::BondAngle, geometry::Point3, spatial_index::SpatialIndex};

pub struct BondGraph {
    bonds: Vec<Vec<Bond>>,
    angles: Vec<Vec<BondAngle>>,
}

impl BondGraph {
    pub fn new(atoms: &[Atom], positions: &[Point3]) -> Self {
        let spatial_index = SpatialIndex::new(positions);
        let mut bonds: Vec<Vec<Bond>> = (0..atoms.len()).map(|_| Vec::new()).collect();
        let mut angles: Vec<Vec<BondAngle>> = (0..atoms.len()).map(|_| Vec::new()).collect();
        for (atom_index, atom) in atoms.iter().enumerate() {
            let position = &positions[atom_index];
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
                let in_contact = positions[atom_index].distance(&positions[neighbor])
                    <= atom.covalent_radius() + atoms[neighbor].covalent_radius();
                if neighbor != atom_index && in_contact {
                    bonds[atom_index].push(Bond::new(atom_index, neighbor));
                }
            }
            match bonds[atom_index].len() {
                0 | 1 => {}
                bonds_len => {
                    for i in 0..(bonds_len - 1) {
                        for j in (i + 1)..bonds_len {
                            angles[atom_index].push(BondAngle::new(
                                bonds[atom_index][i].atom_2,
                                atom_index,
                                bonds[atom_index][j].atom_2,
                            ))
                        }
                    }
                }
            }
        }
        Self { bonds, angles }
    }

    pub fn bond_count(&self) -> usize {
        self.bonds.iter().map(Vec::len).sum()
    }

    pub fn angles(&self) -> impl Iterator<Item = &BondAngle> {
        self.angles.iter().flatten()
    }
}
