use crate::structure::{Atom, Entity, EntityInstance, Residue, Structure};

#[derive(Debug)]
pub struct Selection {
    pub atoms: Vec<bool>,
}

impl Selection {
    pub fn by_instance(structure: &Structure, keep: impl Fn(&EntityInstance, &Entity) -> bool) -> Selection {
        let mut atoms = vec![false; structure.atoms.len()];
        let selected_residues = structure
            .instances
            .iter()
            .filter(|instance| keep(instance, &structure.entities[instance.entity]))
            .flat_map(|instance| &structure.residues[instance.residues.clone()]);
        for residue in selected_residues {
            atoms[residue.atoms.clone()].fill(true);
        }
        Selection { atoms }
    }

    pub fn by_residue(structure: &Structure, keep: impl Fn(&Residue) -> bool) -> Selection {
        let mut atoms = vec![false; structure.atoms.len()];
        let selected_residues = structure.residues.iter().filter(|residue| keep(residue));
        for residue in selected_residues {
            atoms[residue.atoms.clone()].fill(true);
        }
        Selection { atoms }
    }

    pub fn by_atom(structure: &Structure, keep: impl Fn(&Atom) -> bool) -> Selection {
        Selection {
            atoms: structure.atoms.iter().map(keep).collect(),
        }
    }

    pub fn residues_touched(&self, structure: &Structure) -> Vec<usize> {
        assert_eq!(
            self.atoms.len(),
            structure.atoms.len(),
            "selection is from a different structure"
        );

        (0..structure.residues.len())
            .filter(|&residue| structure.residues[residue].atoms.clone().any(|atom| self.atoms[atom]))
            .collect()
    }

    pub fn and(&self, other: &Selection) -> Selection {
        assert_eq!(
            self.atoms.len(),
            other.atoms.len(),
            "selections are from different structures"
        );
        let atoms = self.atoms.iter().zip(&other.atoms).map(|(&a, &b)| a && b).collect();
        Selection { atoms }
    }

    pub fn or(&self, other: &Selection) -> Selection {
        assert_eq!(
            self.atoms.len(),
            other.atoms.len(),
            "selections are from different structures"
        );
        let atoms = self.atoms.iter().zip(&other.atoms).map(|(&a, &b)| a || b).collect();
        Selection { atoms }
    }

    pub fn not(&self) -> Selection {
        let atoms = self.atoms.iter().map(|atom| !atom).collect();
        Selection { atoms }
    }
}
