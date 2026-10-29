use std::ops::Range;

use crate::{bond::Bond, element::Element};

#[derive(Debug)]
pub struct Structure {
    pub atoms: Vec<Atom>,
    pub residues: Vec<Residue>,
    pub instances: Vec<EntityInstance>,
    pub entities: Vec<Entity>,
}

#[derive(Debug)]
pub struct Atom {
    pub element: Element,
    pub name: String,
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
