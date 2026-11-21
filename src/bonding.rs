use std::{collections::HashMap, vec};

use crate::{
    bond::{Bond, BondKind, BondSource, Link, LinkType, ResidueTemplateBond},
    bond_order::BondOrder,
    element::Element,
    functional_groups::functional_groups,
    graph::Graph,
    protonation::protonate,
    ring::Ring,
    structure::{EntityKind, PolymerKind, Structure},
};

pub struct Partner {
    pub element: Element,
    pub order: BondOrder,
    pub aromatic: bool,
    pub atom: Option<usize>, // None when it is missing from the file
}

#[derive(Debug)]
pub struct MissingResiudeTemplateError {
    name: String,
}

#[derive(Debug)]
pub struct MissingPartner {
    pub name: String,
    pub element: Element,
    pub order: BondOrder,
    pub aromatic: bool,
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
    pub rings: Vec<Ring>,
}

impl Bonding {
    pub fn bond_order_sum(&self, atom: usize) -> u8 {
        let present: u8 = self
            .bonds
            .neighbors(atom)
            .iter()
            .map(|&(_, edge)| match self.bonds.edges()[edge].kind {
                BondKind::Coordination => 0,
                BondKind::Covalent { order, .. } => order.as_u8(),
            })
            .sum();
        let missing: u8 = self.missing_partners[atom].iter().map(|partner| partner.order.as_u8()).sum();
        present + self.implicit_hydrogens[atom] + missing
    }

    pub fn atoms_bonded_sum(&self, atom: usize) -> u8 {
        let present =
            self.bonds.neighbors(atom).iter().filter(|&&(_, edge)| matches!(self.bonds.edges()[edge].kind, BondKind::Covalent { .. })).count() as u8;
        let missing = self.missing_partners[atom].len() as u8;
        present + self.implicit_hydrogens[atom] + missing
    }

    pub fn new(
        structure: &Structure,
        templates: &HashMap<String, Vec<ResidueTemplateBond>>,
        links: &[Link],
        ph: f64,
    ) -> Result<(Self, Vec<LinkIssue>), MissingResiudeTemplateError> {
        let mut bonding = template_bonds(structure, templates)?;
        let mut issues = backbone_bonds(&structure, &mut bonding);
        issues.extend(link_bonds(&mut bonding, links));
        bonding.rings = rings(&bonding);
        let groups = functional_groups(&structure, &bonding);
        protonate(&structure, &mut bonding, &groups, ph);
        Ok((bonding, issues))
    }
}

pub fn covalent_heavy_partners<'a>(structure: &'a Structure, bonding: &'a Bonding, atom: usize) -> impl Iterator<Item = Partner> + 'a {
    covalent_partners(structure, bonding, atom).filter(|partner| partner.element != Element::H)
}

pub fn covalent_partners<'a>(structure: &'a Structure, bonding: &'a Bonding, atom: usize) -> impl Iterator<Item = Partner> + 'a {
    let covalent_neighbors = bonding.bonds.neighbors(atom).iter().filter_map(|&(neighbor, edge)| match bonding.bonds.edges()[edge].kind {
        BondKind::Coordination => None,
        BondKind::Covalent { order, aromatic, .. } => {
            Some(Partner { element: structure.atoms[neighbor].element, atom: Some(neighbor), order, aromatic })
        }
    });
    let covalent_missing = bonding.missing_partners[atom].iter().map(|missing_partner| Partner {
        element: missing_partner.element,
        atom: None,
        order: missing_partner.order,
        aromatic: missing_partner.aromatic,
    });
    covalent_neighbors.chain(covalent_missing)
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
    if let BondKind::Covalent { order, .. } = bond.kind {
        let (atom_1, atom_2) = (bond.atom_1, bond.atom_2);
        let payment_1 = payment(bonding, atom_1, order).ok_or(LinkIssueKind::NoRoom { atom: atom_1, order })?;
        let payment_2 = payment(bonding, atom_2, order).ok_or(LinkIssueKind::NoRoom { atom: atom_2, order })?;
        pay(bonding, atom_1, order, payment_1);
        pay(bonding, atom_2, order, payment_2);
    }
    bonding.bonds.add_edge(bond);
    Ok(())
}

fn template_bonds(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>) -> Result<Bonding, MissingResiudeTemplateError> {
    let mut bonds = vec![];
    let mut implicit_hydrogens = vec![0; structure.atoms.len()];
    let mut missing_partners: Vec<Vec<MissingPartner>> = (0..structure.atoms.len()).map(|_| vec![]).collect();
    for (residue_id, residue) in structure.residues.iter().enumerate() {
        // Single-atom residues, such as metal ions, have no bonds and so no template.
        if residue.atoms.len() == 1 && !templates.contains_key(&residue.name) {
            continue;
        }
        let residue_templates = templates.get(&residue.name).ok_or_else(|| MissingResiudeTemplateError { name: residue.name.clone() })?;
        for template in residue_templates {
            let residue = &structure.residues[residue_id];
            let atom_1 = structure.atom_named(residue, &template.atom_name_1);
            let atom_2 = structure.atom_named(residue, &template.atom_name_2);
            // A bond to a metal doesn't use up valence, so it is kept apart from the covalent bonds.
            if template.element_1.is_metal() || template.element_2.is_metal() {
                if let (Some(atom_1), Some(atom_2)) = (atom_1, atom_2) {
                    bonds.push(Bond { atom_1, atom_2, kind: BondKind::Coordination });
                }
                continue;
            }
            let (present, name, element) = match (atom_1, atom_2) {
                (Some(atom_1), Some(atom_2)) => {
                    bonds.push(Bond {
                        atom_1,
                        atom_2,
                        kind: BondKind::Covalent { order: template.order, aromatic: template.aromatic, source: BondSource::Template },
                    });
                    continue;
                }
                (Some(present), None) => (present, &template.atom_name_2, template.element_2),
                (None, Some(present)) => (present, &template.atom_name_1, template.element_1),
                (None, None) => continue,
            };
            if element == Element::H {
                implicit_hydrogens[present] += 1;
            } else {
                missing_partners[present].push(MissingPartner { name: name.clone(), element, order: template.order, aromatic: template.aromatic });
            }
        }
    }
    Ok(Bonding { bonds: Graph::new(structure.atoms.len(), bonds), implicit_hydrogens, missing_partners, rings: vec![] })
}

fn backbone_bonds(structure: &Structure, bonding: &mut Bonding) -> Vec<LinkIssue> {
    let mut issues = vec![];
    for instance in &structure.instances {
        let Some(((name_1, element_1), name_2)) = backbone_link(structure.entities[instance.entity].kind) else {
            continue;
        };
        let residues = &structure.residues[instance.residues.clone()];
        for (position, next) in residues.iter().enumerate() {
            let (Some(next_seq), Some(atom_2)) = (next.seq_id, structure.atom_named(next, name_2)) else {
                continue;
            };
            // label_seq_id counts unmodelled residues, so a gap means a chain break.
            let previous = position.checked_sub(1).map(|position| &residues[position]).filter(|previous| previous.seq_id == Some(next_seq - 1));
            match previous.and_then(|previous| structure.atom_named(previous, name_1)) {
                Some(atom_1) => {
                    let bond =
                        Bond { atom_1, atom_2, kind: BondKind::Covalent { order: BondOrder::Single, aromatic: false, source: BondSource::Backbone } };
                    if let Err(kind) = add_bond(bonding, bond) {
                        issues.push(LinkIssue { atom_1, atom_2, kind });
                    }
                }
                // The neighbour is in the molecule but not in the file, so it still takes the place of a hydrogen.
                None if next_seq > 1 => {
                    if let Some(Payment::Hydrogens) = payment(bonding, atom_2, BondOrder::Single) {
                        pay(bonding, atom_2, BondOrder::Single, Payment::Hydrogens);
                        bonding.missing_partners[atom_2].push(MissingPartner {
                            name: name_1.to_string(),
                            element: element_1,
                            order: BondOrder::Single,
                            aromatic: false,
                        });
                    }
                }
                None => {}
            }
        }
    }
    issues
}

fn link_bonds(bonding: &mut Bonding, links: &[Link]) -> Vec<LinkIssue> {
    let mut issues = vec![];
    for link in links {
        let source = match link.link_type {
            // We don't consider the metal coordination as bonds since they can change withing different geometries.
            LinkType::MetalCoordination | LinkType::HydrogenBond => continue,
            LinkType::Disulfide => BondSource::Disulfide,
            LinkType::Covalent => BondSource::CovalentLink,
        };
        let (atom_1, atom_2) = (link.atom_1, link.atom_2);
        // A link can repeat a bond that already exists, such as the peptide bond next to a non-standard residue.
        if bonding.bonds.neighbors(atom_1).iter().any(|&(other, _)| other == atom_2) {
            continue;
        }
        // Each bond is added as soon as it is decided, so a later link sees it.
        let added = link_order(bonding, link)
            .and_then(|order| add_bond(bonding, Bond { atom_1, atom_2, kind: BondKind::Covalent { order, aromatic: false, source } }));
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
        match bond.kind {
            BondKind::Covalent { order, aromatic, .. } => order != BondOrder::Single && !aromatic,
            BondKind::Coordination => false,
        }
    })
}

pub fn rings(bonding: &Bonding) -> Vec<Ring> {
    bonding
        .bonds
        .rings_where(|edge| {
            let bond = &bonding.bonds.edges()[edge];
            match bond.kind {
                BondKind::Covalent { source, .. } => !matches!(source, BondSource::Disulfide | BondSource::CovalentLink),
                BondKind::Coordination => false,
            }
        })
        .iter()
        .map(|ring| {
            let aromatic = ring.edges.iter().all(|&edge| {
                let bond = &bonding.bonds.edges()[edge];
                match bond.kind {
                    BondKind::Covalent { order: _, aromatic, source: _ } => aromatic,
                    BondKind::Coordination => false,
                }
            });
            Ring { atoms: ring.nodes.clone(), aromatic }
        })
        .collect()
}

fn backbone_link(entity_kind: EntityKind) -> Option<((&'static str, Element), &'static str)> {
    match entity_kind {
        EntityKind::Polymer(PolymerKind::PolypeptideL | PolymerKind::PolypeptideD) => Some((("C", Element::C), "N")),
        EntityKind::Polymer(PolymerKind::Dna | PolymerKind::Rna | PolymerKind::DnaRnaHybrid) => Some((("O3'", Element::O), "P")),
        EntityKind::Polymer(PolymerKind::CyclicPseudoPeptide | PolymerKind::PeptideNucleicAcid | PolymerKind::Other) => None,
        EntityKind::NonPolymer | EntityKind::Branched | EntityKind::Water => None,
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::structure::Atom;

    fn bond(atom_1: usize, atom_2: usize, aromatic: bool, source: BondSource) -> Bond {
        Bond { atom_1, atom_2, kind: BondKind::Covalent { order: BondOrder::Single, aromatic, source } }
    }

    fn structure_with(atom_count: usize, bonds: Vec<Bond>) -> (Structure, Bonding) {
        let atoms = (0..atom_count).map(|_| Atom { element: Element::C, name: "C".to_string(), stated_charge: None }).collect();
        let bonding = Bonding {
            bonds: Graph::new(atom_count, bonds),
            implicit_hydrogens: vec![0; atom_count],
            missing_partners: (0..atom_count).map(|_| vec![]).collect(),
            rings: vec![],
        };
        (Structure { atoms, residues: vec![], instances: vec![], entities: vec![] }, bonding)
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
        let (_, bonding) = structure_with(11, bonds);
        let rings = rings(&bonding);

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
        let (_, bonding) = structure_with(4, bonds);
        let rings = rings(&bonding);
        assert!(rings.is_empty());
    }
}
