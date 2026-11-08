use crate::{
    mask::Mask,
    structure::{Atom, Entity, EntityInstance, Residue, Structure},
};

pub type Selection = Mask<Atom>;

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
        Mask::from_bools(atoms)
    }

    pub fn by_residue(structure: &Structure, keep: impl Fn(&Residue) -> bool) -> Selection {
        let mut atoms = vec![false; structure.atoms.len()];
        let selected_residues = structure.residues.iter().filter(|residue| keep(residue));
        for residue in selected_residues {
            atoms[residue.atoms.clone()].fill(true);
        }
        Mask::from_bools(atoms)
    }

    pub fn by_atom(structure: &Structure, keep: impl Fn(&Atom) -> bool) -> Selection {
        Mask::from_bools(structure.atoms.iter().map(keep).collect())
    }

    pub fn residues_touched(&self, structure: &Structure) -> Vec<usize> {
        assert_eq!(self.len(), structure.atoms.len(), "selection is from a different structure");
        (0..structure.residues.len()).filter(|&residue| structure.residues[residue].atoms.clone().any(|atom| self.contains(atom))).collect()
    }
}
