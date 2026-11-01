use std::collections::{HashMap, HashSet};

use crate::{
    bond::{Bond, Link, LinkType, ResidueTemplateBond},
    bond_order::BondOrder,
    graph::Graph,
    structure::{self, EntityKind, PolymerKind, Residue, Structure},
};

pub fn new(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>, links: &[Link]) -> Graph<Bond> {
    let mut bonds = vec![];
    let atom_named = |residue: &Residue, name: &str| residue.atoms.clone().find(|&atom| structure.atoms[atom].name == name);
    for instance in &structure.instances {
        let residues = &structure.residues[instance.residues.clone()];
        for residue in residues {
            for template in templates.get(&residue.name).into_iter().flatten() {
                if let (Some(atom_1), Some(atom_2)) = (atom_named(residue, &template.atom_name_1), atom_named(residue, &template.atom_name_2)) {
                    bonds.push(Bond { atom_1, atom_2, order: template.order, aromatic: template.aromatic });
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
                    bonds.push(Bond { atom_1, atom_2, order: BondOrder::Single, aromatic: false });
                }
            }
        }
    }
    // A link can repeat a bond the backbone rule already made, such as the peptide bond next to a non-standard residue.
    let bonded: HashSet<(usize, usize)> = bonds.iter().map(|bond| atom_pair(bond.atom_1, bond.atom_2)).collect();
    for link in links {
        if bonded.contains(&atom_pair(link.atom_1, link.atom_2)) {
            continue;
        }
        if let Some(order) = link_order(link) {
            bonds.push(Bond { atom_1: link.atom_1, atom_2: link.atom_2, order, aromatic: false });
        }
    }
    Graph::new(structure.atoms.len(), bonds)
}

pub fn valence_check(structure: Structure) -> bool {}

// None for links that aren't covalent bonds, and for covalent links whose order is still unknown.
fn link_order(link: &Link) -> Option<BondOrder> {
    match link.link_type {
        // A disulfide is always a single S–S bond.
        LinkType::Disulfide => Some(link.order.unwrap_or(BondOrder::Single)),
        LinkType::Covalent => link.order,
        LinkType::MetalCoordination | LinkType::HydrogenBond => None,
    }
}

fn atom_pair(atom_1: usize, atom_2: usize) -> (usize, usize) {
    (atom_1.min(atom_2), atom_1.max(atom_2))
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
