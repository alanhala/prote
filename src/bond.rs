use std::collections::HashMap;

use crate::{
    bond_order::BondOrder,
    structure::{EntityKind, PolymerKind, Residue, Structure},
};

pub struct ResidueTemplateBond {
    pub atom_name_1: String,
    pub atom_name_2: String,
    pub order: BondOrder,
    pub aromatic: bool,
}

#[derive(Debug)]
pub struct Bond {
    pub atom_1: usize,
    pub atom_2: usize,
    pub order: BondOrder,
    pub aromatic: bool,
}

// Polymers link consecutive residues through fixed atoms: peptides C(i)–N(i+1), nucleic acids O3'(i)–P(i+1).
fn backbone_link(entity_kind: EntityKind) -> Option<(&'static str, &'static str)> {
    match entity_kind {
        EntityKind::Polymer(PolymerKind::PolypeptideL | PolymerKind::PolypeptideD) => Some(("C", "N")),
        EntityKind::Polymer(PolymerKind::Dna | PolymerKind::Rna | PolymerKind::DnaRnaHybrid) => Some(("O3'", "P")),
        // No single link we can rely on, so we leave these unlinked rather than guess.
        EntityKind::Polymer(PolymerKind::CyclicPseudoPeptide | PolymerKind::PeptideNucleicAcid | PolymerKind::Other) => None,
        EntityKind::NonPolymer | EntityKind::Branched | EntityKind::Water => None,
    }
}

impl Bond {
    pub fn new(atom_1: usize, atom_2: usize, order: BondOrder, aromatic: bool) -> Self {
        Self { atom_1, atom_2, order, aromatic }
    }

    // pub fn bonds(indexed: &IndexedStructure, atoms: &Selection) -> Vec<Bond> {
    //     indexed
    //         // PMC3557220, eq. 1: 0.8 < d <= r_1 + r_2 + 0.4
    //         .pairs_within(atoms, atoms, |atom| atom.element.covalent_radius() + 0.2)
    //         .filter(|&(_, _, distance)| distance > 0.8)
    //         .map(|(atom_1, atom_2, distance)| Bond::new(atom_1, atom_2, distance))
    //         .collect()
    // }
    //
    pub fn bonds(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>) -> Vec<Bond> {
        let mut bonds = vec![];
        let atom_named = |residue: &Residue, name: &str| residue.atoms.clone().find(|&atom| structure.atoms[atom].name == name);
        for instance in &structure.instances {
            let residues = &structure.residues[instance.residues.clone()];
            for residue in residues {
                for template in templates.get(&residue.name).into_iter().flatten() {
                    if let (Some(atom_1), Some(atom_2)) = (atom_named(residue, &template.atom_name_1), atom_named(residue, &template.atom_name_2)) {
                        bonds.push(Bond::new(atom_1, atom_2, template.order, template.aromatic));
                    }
                }
            }
            if let Some((name_1, name_2)) = backbone_link(structure.entities[instance.entity].kind) {
                for pair in residues.windows(2) {
                    let (previous, next) = (&pair[0], &pair[1]);
                    // label_seq_id counts unmodelled residues, so a gap means a chain break.
                    if !previous.seq_id.is_some_and(|seq_id| Some(seq_id + 1) == next.seq_id) {
                        continue;
                    }
                    if let (Some(atom_1), Some(atom_2)) = (atom_named(previous, name_1), atom_named(next, name_2)) {
                        bonds.push(Bond::new(atom_1, atom_2, BondOrder::Single, false));
                    }
                }
            }
        }
        bonds
    }
}
