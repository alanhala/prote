use std::collections::{HashMap, HashSet};

use crate::{
    bond::{Bond, BondSource, Link, LinkType, ResidueTemplateBond},
    bond_order::BondOrder,
    element::Element,
    graph::Graph,
    ring::Ring,
    structure::{self, EntityKind, PolymerKind, Residue, Structure},
};

#[derive(Debug)]
pub struct MissingResiudeTemplateError {
    name: String,
}

#[derive(Debug)]
pub struct MissingPartner {
    pub name: String,
    pub element: Element,
    pub order: BondOrder,
}

#[derive(Debug)]
pub struct LinkIssue {
    pub atom_1: usize,
    pub atom_2: usize,
    pub kind: LinkIssueKind,
}

#[derive(Debug)]
pub enum LinkIssueKind {
    NoRoom {
        atom: usize,
        order: BondOrder,
    },
    /// Several heavy partners are missing, so we can't tell which one the link replaced.
    SeveralMissingHeavy,
    /// Both ends are missing a heavy partner, with different orders, so no single order fills both.
    MismatchedMissingHeavy,
    /// The atom keeps a non-aromatic multiple bond to a present atom: probably an addition across it,
    /// which would make that template bond wrong for this structure.
    AdditionAcrossMultipleBond,
}

enum Payment {
    MissingPartner(usize),
    Hydrogens,
}

#[derive(Debug)]
pub struct Bonding {
    pub bonds: Graph<Bond>,
    pub implicit_hydrogens: Vec<u8>,
    pub missing_partners: Vec<Vec<MissingPartner>>,
    pub coordination: Vec<(usize, usize)>,
}

pub fn with_bonds(
    structure: Structure,
    templates: &HashMap<String, Vec<ResidueTemplateBond>>,
    links: &[Link],
) -> Result<(Structure, Vec<LinkIssue>), MissingResiudeTemplateError> {
    let mut bonding = template_bonds(&structure, templates)?;
    let mut issues = backbone_bonds(&structure, &mut bonding);
    issues.extend(link_bonds(&mut bonding, links));
    Ok((Structure { bonding, ..structure }, issues))
}

fn payment(bonding: &Bonding, atom: usize, order: BondOrder) -> Option<Payment> {
    if let Some(position) = bonding.missing_partners[atom].iter().position(|partner| partner.order == order) {
        Some(Payment::MissingPartner(position))
    } else if bonding.implicit_hydrogens[atom] >= order.as_u8() {
        Some(Payment::Hydrogens)
    } else {
        None
    }
}

fn pay(bonding: &mut Bonding, atom: usize, order: BondOrder, payment: Payment) {
    match payment {
        Payment::MissingPartner(position) => {
            bonding.missing_partners[atom].remove(position);
        }
        Payment::Hydrogens => bonding.implicit_hydrogens[atom] -= order.as_u8(),
    }
}

// Makes room at both atoms and adds the bond, or changes nothing when either atom has no room.
fn add_bond(bonding: &mut Bonding, bond: Bond) -> Result<(), LinkIssueKind> {
    let (atom_1, atom_2, order) = (bond.atom_1, bond.atom_2, bond.order);
    let payment_1 = payment(bonding, atom_1, order).ok_or(LinkIssueKind::NoRoom { atom: atom_1, order })?;
    let payment_2 = payment(bonding, atom_2, order).ok_or(LinkIssueKind::NoRoom { atom: atom_2, order })?;
    pay(bonding, atom_1, order, payment_1);
    pay(bonding, atom_2, order, payment_2);
    bonding.bonds.add_edge(bond);
    Ok(())
}

fn atom_named(structure: &Structure, residue: &Residue, name: &str) -> Option<usize> {
    residue.atoms.clone().find(|&atom| structure.atoms[atom].name == name)
}

fn template_bonds(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>) -> Result<Bonding, MissingResiudeTemplateError> {
    let mut bonds = vec![];
    let mut implicit_hydrogens = vec![0; structure.atoms.len()];
    let mut missing_partners: Vec<Vec<MissingPartner>> = (0..structure.atoms.len()).map(|_| vec![]).collect();
    let mut coordination = vec![];
    for (residue_id, residue) in structure.residues.iter().enumerate() {
        // Single-atom residues, such as metal ions, have no bonds and so no template.
        if residue.atoms.len() == 1 && !templates.contains_key(&residue.name) {
            continue;
        }
        let residue_templates = templates.get(&residue.name).ok_or_else(|| MissingResiudeTemplateError { name: residue.name.clone() })?;
        for template in residue_templates {
            let residue = &structure.residues[residue_id];
            let atom_1 = atom_named(&structure, residue, &template.atom_name_1);
            let atom_2 = atom_named(&structure, residue, &template.atom_name_2);
            // A bond to a metal doesn't use up valence, so it is kept apart from the covalent bonds.
            if template.element_1.is_metal() || template.element_2.is_metal() {
                if let (Some(atom_1), Some(atom_2)) = (atom_1, atom_2) {
                    coordination.push((atom_1, atom_2));
                }
                continue;
            }
            let (present, name, element) = match (atom_1, atom_2) {
                (Some(atom_1), Some(atom_2)) => {
                    bonds.push(Bond { atom_1, atom_2, order: template.order, aromatic: template.aromatic, source: BondSource::Template });
                    continue;
                }
                (Some(present), None) => (present, &template.atom_name_2, template.element_2),
                (None, Some(present)) => (present, &template.atom_name_1, template.element_1),
                (None, None) => continue,
            };
            if element == Element::H {
                implicit_hydrogens[present] += 1;
            } else {
                missing_partners[present].push(MissingPartner { name: name.clone(), element, order: template.order });
            }
        }
    }
    Ok(Bonding { bonds: Graph::new(structure.atoms.len(), bonds), implicit_hydrogens, missing_partners, coordination })
}

fn backbone_bonds(structure: &Structure, bonding: &mut Bonding) -> Vec<LinkIssue> {
    let mut issues = vec![];
    for instance in &structure.instances {
        let Some((name_1, name_2)) = backbone_link(structure.entities[instance.entity].kind) else {
            continue;
        };
        for pair in structure.residues[instance.residues.clone()].windows(2) {
            let (previous, next) = (&pair[0], &pair[1]);
            let (Some(previous_seq), Some(next_seq)) = (previous.seq_id, next.seq_id) else {
                continue;
            };
            // label_seq_id counts unmodelled residues, so a gap means a chain break.
            if previous_seq + 1 != next_seq {
                continue;
            }
            let (Some(atom_1), Some(atom_2)) = (atom_named(&structure, previous, name_1), atom_named(&structure, next, name_2)) else {
                continue;
            };
            let bond = Bond { atom_1, atom_2, order: BondOrder::Single, aromatic: false, source: BondSource::Backbone };
            if let Err(kind) = add_bond(bonding, bond) {
                issues.push(LinkIssue { atom_1, atom_2, kind });
            }
        }
    }
    issues
}

fn link_bonds(bonding: &mut Bonding, links: &[Link]) -> Vec<LinkIssue> {
    let mut issues = vec![];
    for link in links {
        let source = match link.link_type {
            LinkType::MetalCoordination => {
                bonding.coordination.push((link.atom_1, link.atom_2));
                continue;
            }
            LinkType::HydrogenBond => continue,
            LinkType::Disulfide => BondSource::Disulfide,
            LinkType::Covalent => BondSource::CovalentLink,
        };
        let (atom_1, atom_2) = (link.atom_1, link.atom_2);
        // A link can repeat a bond that already exists, such as the peptide bond next to a non-standard residue.
        if bonding.bonds.neighbors(atom_1).iter().any(|&(other, _)| other == atom_2) {
            continue;
        }
        // Each bond is added as soon as it is decided, so a later link sees it.
        let added = link_order(bonding, link).and_then(|order| add_bond(bonding, Bond { atom_1, atom_2, order, aromatic: false, source }));
        if let Err(kind) = added {
            issues.push(LinkIssue { atom_1, atom_2, kind });
        }
    }
    issues
}

fn link_order(bonding: &Bonding, link: &Link) -> Result<BondOrder, LinkIssueKind> {
    if let Some(order) = link.order {
        return Ok(order);
    }
    if link.link_type == LinkType::Disulfide {
        return Ok(BondOrder::Single);
    }
    let order = match (bonding.missing_partners[link.atom_1].as_slice(), bonding.missing_partners[link.atom_2].as_slice()) {
        ([], []) => BondOrder::Single,
        ([partner], []) | ([], [partner]) => partner.order,
        ([partner_1], [partner_2]) if partner_1.order == partner_2.order => partner_1.order,
        ([_], [_]) => return Err(LinkIssueKind::MismatchedMissingHeavy),
        _ => return Err(LinkIssueKind::SeveralMissingHeavy),
    };
    for atom in [link.atom_1, link.atom_2] {
        if matches!(payment(bonding, atom, order), Some(Payment::Hydrogens)) && keeps_multiple_bond(bonding, atom) {
            return Err(LinkIssueKind::AdditionAcrossMultipleBond);
        }
    }
    Ok(order)
}

fn keeps_multiple_bond(bonding: &Bonding, atom: usize) -> bool {
    bonding.bonds.neighbors(atom).iter().any(|&(_, edge)| {
        let bond = &bonding.bonds.edges()[edge];
        bond.order != BondOrder::Single && !bond.aromatic
    })
}

pub fn rings(structure: &Structure) -> Vec<Ring> {
    structure
        .bonding
        .bonds
        .rings_where(|edge| !matches!(structure.bonding.bonds.edges()[edge].source, BondSource::Disulfide | BondSource::CovalentLink))
        .iter()
        .map(|ring| {
            let aromatic = ring.edges.iter().all(|&edge| structure.bonding.bonds.edges()[edge].aromatic);
            Ring { atoms: ring.nodes.clone(), aromatic }
        })
        .collect()
}

fn backbone_link(entity_kind: EntityKind) -> Option<(&'static str, &'static str)> {
    match entity_kind {
        EntityKind::Polymer(PolymerKind::PolypeptideL | PolymerKind::PolypeptideD) => Some(("C", "N")),
        EntityKind::Polymer(PolymerKind::Dna | PolymerKind::Rna | PolymerKind::DnaRnaHybrid) => Some(("O3'", "P")),
        EntityKind::Polymer(PolymerKind::CyclicPseudoPeptide | PolymerKind::PeptideNucleicAcid | PolymerKind::Other) => None,
        EntityKind::NonPolymer | EntityKind::Branched | EntityKind::Water => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::Atom;

    fn bond(atom_1: usize, atom_2: usize, aromatic: bool, source: BondSource) -> Bond {
        Bond { atom_1, atom_2, order: BondOrder::Single, aromatic, source }
    }

    fn structure_with(atom_count: usize, bonds: Vec<Bond>) -> Structure {
        let atoms = (0..atom_count).map(|_| Atom { element: Element::C, name: "C".to_string() }).collect();
        let bonding = Bonding {
            bonds: Graph::new(atom_count, bonds),
            implicit_hydrogens: vec![0; atom_count],
            missing_partners: (0..atom_count).map(|_| vec![]).collect(),
            coordination: vec![],
        };
        Structure { atoms, residues: vec![], instances: vec![], entities: vec![], bonding }
    }

    #[test]
    fn rings_are_aromatic_only_when_every_bond_is() {
        let mut bonds = vec![];
        // Benzene-like ring 0..6: every bond aromatic.
        for i in 0..6 {
            bonds.push(bond(i, (i + 1) % 6, true, BondSource::Template));
        }
        // Five-ring 6..11: aromatic except the last bond, as a check on "every bond".
        for i in 0..5 {
            bonds.push(bond(6 + i, 6 + (i + 1) % 5, i < 4, BondSource::Template));
        }
        let rings = rings(&structure_with(11, bonds));

        let mut found: Vec<(usize, bool)> = rings.iter().map(|ring| (ring.atoms.len(), ring.aromatic)).collect();
        found.sort();
        assert_eq!(found, vec![(5, false), (6, true)]);
    }

    #[test]
    fn loops_closed_by_links_are_not_rings() {
        // Chain 0-1-2-3 closed back to 0 by a disulfide, like Cys22–Cys157 in 3PTB.
        let bonds = vec![
            bond(0, 1, false, BondSource::Template),
            bond(1, 2, false, BondSource::Backbone),
            bond(2, 3, false, BondSource::Template),
            bond(3, 0, false, BondSource::Disulfide),
        ];
        assert!(rings(&structure_with(4, bonds)).is_empty());
    }
}
