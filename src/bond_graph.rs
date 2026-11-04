use std::collections::{HashMap, HashSet};

use crate::{
    bond::{Bond, Link, LinkType, ResidueTemplateBond},
    bond_order::BondOrder,
    element::Element,
    graph::Graph,
    structure::{EntityKind, PolymerKind, Residue, Structure},
};

pub fn new(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>, links: &[Link]) -> Result<Graph<Bond>, Vec<LinkIssue>> {
    let mut bonds = vec![];
    let mut backbone_atoms = HashSet::new();
    for instance in &structure.instances {
        let residues = &structure.residues[instance.residues.clone()];
        for residue in residues {
            for template in templates.get(&residue.name).into_iter().flatten() {
                if let (Some(atom_1), Some(atom_2)) = (atom_named(structure, residue, &template.atom_name_1), atom_named(structure, residue, &template.atom_name_2)) {
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
                if let (Some(atom_1), Some(atom_2)) = (atom_named(structure, previous, name_1), atom_named(structure, next, name_2)) {
                    bonds.push(Bond { atom_1, atom_2, order: BondOrder::Single, aromatic: false });
                    backbone_atoms.extend([atom_1, atom_2]);
                }
            }
        }
    }
    // A link can repeat a bond the backbone rule already made, such as the peptide bond next to a non-standard residue.
    let bonded: HashSet<(usize, usize)> = bonds.iter().map(|bond| atom_pair(bond.atom_1, bond.atom_2)).collect();
    let new_links: Vec<&Link> = links.iter().filter(|link| !bonded.contains(&atom_pair(link.atom_1, link.atom_2))).collect();
    // Metal and hydrogen bonds don't take the place of template partners, so only covalent links count.
    let mut links_at: HashMap<usize, usize> = HashMap::new();
    for link in new_links.iter().filter(|link| matches!(link.link_type, LinkType::Covalent | LinkType::Disulfide)) {
        *links_at.entry(link.atom_1).or_default() += 1;
        *links_at.entry(link.atom_2).or_default() += 1;
    }
    let mut issues = vec![];
    for link in new_links {
        match link_order(structure, templates, &backbone_atoms, &links_at, link) {
            Ok(Some(order)) => bonds.push(Bond { atom_1: link.atom_1, atom_2: link.atom_2, order, aromatic: false }),
            Ok(None) => {}
            Err(kind) => issues.push(LinkIssue { atom_1: link.atom_1, atom_2: link.atom_2, kind }),
        }
    }
    if !issues.is_empty() {
        return Err(issues);
    }
    Ok(Graph::new(structure.atoms.len(), bonds))
}

// A covalent link whose order we can't infer yet. Each kind is a case to study, not one to guess at.
#[derive(Debug)]
pub struct LinkIssue {
    pub atom_1: usize,
    pub atom_2: usize,
    pub kind: LinkIssueKind,
}

#[derive(Debug)]
pub enum LinkIssueKind {
    // The atom also has a backbone bond, which may have taken one of its missing partners.
    SharesBackbone { atom: usize },
    // Several links start at the atom, so its missing partners are shared between them.
    SeveralLinks { atom: usize, links: usize },
    // Several heavy partners are missing, so we can't tell which one the link replaced.
    SeveralMissingHeavy { atom: usize, orders: Vec<BondOrder> },
    // Both ends lost a heavy partner, with different orders, so no single order fills both.
    MismatchedMissingHeavy { orders: (BondOrder, BondOrder) },
    // The end that gives up hydrogens has fewer missing ones than the link needs.
    NotEnoughHydrogens { atom: usize, needed: u8, available: usize },
    // The atom keeps a non-aromatic multiple bond to a present atom: probably an addition across it,
    // which would make that template bond wrong for this structure.
    AdditionAcrossMultipleBond { atom: usize, partner: String },
}

// Assumes links form by condensation: each end either lost one heavy partner, whose bond the link takes
// over, or lost hydrogens. The template's other bonds and hydrogens are taken as unchanged.
fn link_order(
    structure: &Structure,
    templates: &HashMap<String, Vec<ResidueTemplateBond>>,
    backbone_atoms: &HashSet<usize>,
    links_at: &HashMap<usize, usize>,
    link: &Link,
) -> Result<Option<BondOrder>, LinkIssueKind> {
    match (link.link_type, link.order) {
        (LinkType::MetalCoordination | LinkType::HydrogenBond, _) => Ok(None),
        (_, Some(order)) => Ok(Some(order)),
        (LinkType::Disulfide, None) => Ok(Some(BondOrder::Single)),
        (LinkType::Covalent, None) => covalent_order(structure, templates, backbone_atoms, links_at, link).map(Some),
    }
}

fn covalent_order(
    structure: &Structure,
    templates: &HashMap<String, Vec<ResidueTemplateBond>>,
    backbone_atoms: &HashSet<usize>,
    links_at: &HashMap<usize, usize>,
    link: &Link,
) -> Result<BondOrder, LinkIssueKind> {
    // Past these checks, every missing partner at either end belongs to this link.
    for atom in [link.atom_1, link.atom_2] {
        if backbone_atoms.contains(&atom) {
            return Err(LinkIssueKind::SharesBackbone { atom });
        }
        let links = links_at[&atom];
        if links > 1 {
            return Err(LinkIssueKind::SeveralLinks { atom, links });
        }
    }
    let ends = [link.atom_1, link.atom_2].map(|atom| (atom, missing_partners(structure, templates, atom)));
    for (atom, (heavy, _)) in &ends {
        if heavy.len() > 1 {
            return Err(LinkIssueKind::SeveralMissingHeavy { atom: *atom, orders: heavy.clone() });
        }
    }
    let order = match (ends[0].1 .0.first(), ends[1].1 .0.first()) {
        (Some(&order_1), Some(&order_2)) if order_1 == order_2 => order_1,
        (Some(&order_1), Some(&order_2)) => return Err(LinkIssueKind::MismatchedMissingHeavy { orders: (order_1, order_2) }),
        (Some(&order), None) | (None, Some(&order)) => order,
        (None, None) => BondOrder::Single,
    };
    // An end that lost no heavy partner gives up hydrogens instead.
    for (atom, (heavy, hydrogens)) in &ends {
        if !heavy.is_empty() {
            continue;
        }
        if *hydrogens < order.as_u8() as usize {
            return Err(LinkIssueKind::NotEnoughHydrogens { atom: *atom, needed: order.as_u8(), available: *hydrogens });
        }
        if let Some(partner) = kept_multiple_bond(structure, templates, *atom) {
            return Err(LinkIssueKind::AdditionAcrossMultipleBond { atom: *atom, partner });
        }
    }
    Ok(order)
}

// Orders of the atom's template bonds to absent heavy partners, and how many of its template hydrogens are absent.
fn missing_partners(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>, atom: usize) -> (Vec<BondOrder>, usize) {
    let residue = &structure.residues[structure.residue_of(atom)];
    let (hydrogens, heavy): (Vec<_>, Vec<_>) = template_partners(structure, templates, atom)
        .filter(|(_, partner, _)| atom_named(structure, residue, partner).is_none())
        .partition(|(_, _, element)| *element == Element::H);
    (heavy.iter().map(|(bond, _, _)| bond.order).collect(), hydrogens.len())
}

// A present partner the atom has a non-aromatic double or triple bond to, if any.
fn kept_multiple_bond(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>, atom: usize) -> Option<String> {
    let residue = &structure.residues[structure.residue_of(atom)];
    template_partners(structure, templates, atom)
        .find(|(bond, partner, _)| bond.order != BondOrder::Single && !bond.aromatic && atom_named(structure, residue, partner).is_some())
        .map(|(_, partner, _)| partner.to_string())
}

// The atom's template bonds, each with the partner's name and element.
fn template_partners<'a>(
    structure: &'a Structure,
    templates: &'a HashMap<String, Vec<ResidueTemplateBond>>,
    atom: usize,
) -> impl Iterator<Item = (&'a ResidueTemplateBond, &'a str, Element)> {
    let residue = &structure.residues[structure.residue_of(atom)];
    let name = &structure.atoms[atom].name;
    templates.get(&residue.name).into_iter().flatten().filter_map(move |bond| {
        if &bond.atom_name_1 == name {
            Some((bond, bond.atom_name_2.as_str(), bond.element_2))
        } else if &bond.atom_name_2 == name {
            Some((bond, bond.atom_name_1.as_str(), bond.element_1))
        } else {
            None
        }
    })
}

fn atom_named(structure: &Structure, residue: &Residue, name: &str) -> Option<usize> {
    residue.atoms.clone().find(|&atom| structure.atoms[atom].name == name)
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
