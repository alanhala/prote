use std::ops::Range;

use crate::{bond::BondKind, bond_graph::Bonding, bond_order::BondOrder, element::Element, structure};

#[derive(Debug)]
pub struct Structure {
    pub atoms: Vec<Atom>,
    pub residues: Vec<Residue>,
    pub instances: Vec<EntityInstance>,
    pub entities: Vec<Entity>,
    pub bonding: Bonding,
}

pub struct Partner {
    pub element: Element,
    pub order: BondOrder,
    pub aromatic: bool,
    pub atom: Option<usize>, // None when it is missing from the file
}

impl Structure {
    pub fn residue_of(&self, atom: usize) -> usize {
        self.residues.partition_point(|residue| residue.atoms.end <= atom)
    }

    pub fn instance_of(&self, residue: usize) -> usize {
        self.instances.partition_point(|instance| instance.residues.end <= residue)
    }

    pub fn atom_named(&self, residue: &Residue, name: &str) -> Option<usize> {
        residue.atoms.clone().find(|&atom| self.atoms[atom].name == name)
    }

    pub fn covalent_partners(&self, atom: usize) -> impl Iterator<Item = Partner> + '_ {
        let covalent_neighbors =
            self.bonding.bonds.neighbors(atom).iter().filter_map(|&(neighbor, edge)| match self.bonding.bonds.edges()[edge].kind {
                BondKind::Coordination => None,
                BondKind::Covalent { order, aromatic, .. } => {
                    Some(Partner { element: self.atoms[neighbor].element, atom: Some(neighbor), order, aromatic })
                }
            });
        let covalent_missing = self.bonding.missing_partners[atom].iter().map(|missing_partner| Partner {
            element: missing_partner.element,
            atom: None,
            order: missing_partner.order,
            aromatic: missing_partner.aromatic,
        });
        covalent_neighbors.chain(covalent_missing)
    }
}

#[derive(Debug)]
pub struct Atom {
    pub element: Element,
    pub name: String,
    /// The non-zero formal charge the file gives this atom (`pdbx_formal_charge`), if any.
    pub stated_charge: Option<i8>,
}

#[derive(Debug)]
pub struct Residue {
    pub name: String,
    pub number: i64,
    // Position in the entity's sequence (`label_seq_id`), counting residues missing from the model,
    // so a gap between neighbours is a chain break. None for non-polymers.
    pub seq_id: Option<i64>,
    pub atoms: Range<usize>,
}

#[derive(Debug)]
pub struct EntityInstance {
    pub id: String,
    pub chain: String,
    pub entity: usize,
    pub residues: Range<usize>,
}

#[derive(Debug, Clone)]
pub struct Entity {
    pub name: String,
    pub kind: EntityKind,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum EntityKind {
    Polymer(PolymerKind),
    NonPolymer,
    Branched,
    Water,
}

#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PolymerKind {
    PolypeptideL,
    PolypeptideD,
    CyclicPseudoPeptide,
    Dna,
    Rna,
    DnaRnaHybrid,
    PeptideNucleicAcid,
    Other,
}
