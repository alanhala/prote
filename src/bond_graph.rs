use std::collections::{HashMap, HashSet};

use crate::{
    bond::{Bond, Link, LinkType, ResidueTemplateBond},
    bond_order::BondOrder,
    element::Element,
    graph::Graph,
    structure::{EntityKind, PolymerKind, Residue, Structure},
};

pub fn new(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>, links: &[Link]) -> Result<Graph<Bond>, Vec<LinkIssue>> {
    let mut bonds = template_bonds(structure, templates);
    let backbone = backbone_bonds(structure);
    let backbone_atoms: HashSet<usize> = backbone.iter().flat_map(|bond| [bond.atom_1, bond.atom_2]).collect();
    bonds.extend(backbone);
    let link_bonds = link_bonds(structure, templates, &bonds, &backbone_atoms, links)?;
    bonds.extend(link_bonds);
    Ok(Graph::new(structure.atoms.len(), bonds))
}

/// A covalent link whose order we can't infer yet. Each kind is a case to study, not one to guess at.
#[derive(Debug)]
pub struct LinkIssue {
    pub atom_1: usize,
    pub atom_2: usize,
    pub kind: LinkIssueKind,
}

#[derive(Debug)]
pub enum LinkIssueKind {
    /// The atom also has a backbone bond, which may have taken one of its missing partners.
    SharesBackbone { atom: usize },
    /// Several links start at the atom, so its missing partners are shared between them.
    SeveralLinks { atom: usize, links: usize },
    /// The templates have no bonds for the atom, so there is nothing to compare it with.
    NoTemplate { atom: usize },
    /// Several heavy partners are missing, so we can't tell which one the link replaced.
    SeveralMissingHeavy { atom: usize, orders: Vec<BondOrder> },
    /// Both ends are missing a heavy partner, with different orders, so no single order fills both.
    MismatchedMissingHeavy { orders: (BondOrder, BondOrder) },
    /// The end with no missing heavy partner has fewer missing hydrogens than the link needs.
    NotEnoughHydrogens { atom: usize, needed: u8, available: usize },
    /// The atom keeps a non-aromatic multiple bond to a present atom: probably an addition across it,
    /// which would make that template bond wrong for this structure.
    AdditionAcrossMultipleBond { atom: usize, partner: String },
}

fn template_bonds(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>) -> Vec<Bond> {
    structure
        .residues
        .iter()
        .flat_map(|residue| {
            templates.get(&residue.name).into_iter().flatten().filter_map(move |template| {
                let atom_1 = atom_named(structure, residue, &template.atom_name_1)?;
                let atom_2 = atom_named(structure, residue, &template.atom_name_2)?;
                Some(Bond { atom_1, atom_2, order: template.order, aromatic: template.aromatic })
            })
        })
        .collect()
}

fn backbone_bonds(structure: &Structure) -> Vec<Bond> {
    let mut bonds = vec![];
    for instance in &structure.instances {
        let Some((name_1, name_2)) = backbone_link(structure.entities[instance.entity].kind) else {
            continue;
        };
        for pair in structure.residues[instance.residues.clone()].windows(2) {
            let (previous, next) = (&pair[0], &pair[1]);
            // label_seq_id counts unmodelled residues, so a gap means a chain break.
            if !previous.seq_id.is_some_and(|seq_id| Some(seq_id + 1) == next.seq_id) {
                continue;
            }
            if let (Some(atom_1), Some(atom_2)) = (atom_named(structure, previous, name_1), atom_named(structure, next, name_2)) {
                bonds.push(Bond { atom_1, atom_2, order: BondOrder::Single, aromatic: false });
            }
        }
    }
    bonds
}

// Fails with every link we can't order, not just the first, so one run shows all the cases to study.
fn link_bonds(
    structure: &Structure,
    templates: &HashMap<String, Vec<ResidueTemplateBond>>,
    bonds: &[Bond],
    backbone_atoms: &HashSet<usize>,
    links: &[Link],
) -> Result<Vec<Bond>, Vec<LinkIssue>> {
    // A link can repeat a bond the backbone rule already made, such as the peptide bond next to a non-standard residue.
    let bonded: HashSet<(usize, usize)> = bonds.iter().map(|bond| atom_pair(bond.atom_1, bond.atom_2)).collect();
    let unbonded_links: Vec<&Link> = links.iter().filter(|link| !bonded.contains(&atom_pair(link.atom_1, link.atom_2))).collect();
    // Metal and hydrogen bonds don't take the place of template partners, so only covalent links count.
    let mut covalent_links_at: HashMap<usize, usize> = HashMap::new();
    for link in unbonded_links.iter().filter(|link| matches!(link.link_type, LinkType::Covalent | LinkType::Disulfide)) {
        *covalent_links_at.entry(link.atom_1).or_default() += 1;
        *covalent_links_at.entry(link.atom_2).or_default() += 1;
    }
    let mut link_bonds = vec![];
    let mut issues = vec![];
    for link in unbonded_links {
        match link_order(structure, templates, backbone_atoms, &covalent_links_at, link) {
            Ok(Some(order)) => link_bonds.push(Bond { atom_1: link.atom_1, atom_2: link.atom_2, order, aromatic: false }),
            Ok(None) => {}
            Err(kind) => issues.push(LinkIssue { atom_1: link.atom_1, atom_2: link.atom_2, kind }),
        }
    }
    if issues.is_empty() {
        Ok(link_bonds)
    } else {
        Err(issues)
    }
}

// Assumes links form by condensation: each end is missing either one heavy partner, whose bond the link
// takes over, or hydrogens. The template's other bonds and hydrogens are taken as unchanged.
fn link_order(
    structure: &Structure,
    templates: &HashMap<String, Vec<ResidueTemplateBond>>,
    backbone_atoms: &HashSet<usize>,
    covalent_links_at: &HashMap<usize, usize>,
    link: &Link,
) -> Result<Option<BondOrder>, LinkIssueKind> {
    match (link.link_type, link.order) {
        (LinkType::MetalCoordination | LinkType::HydrogenBond, _) => Ok(None),
        (_, Some(order)) => Ok(Some(order)),
        (LinkType::Disulfide, None) => Ok(Some(BondOrder::Single)),
        (LinkType::Covalent, None) => covalent_order(structure, templates, backbone_atoms, covalent_links_at, link).map(Some),
    }
}

fn covalent_order(
    structure: &Structure,
    templates: &HashMap<String, Vec<ResidueTemplateBond>>,
    backbone_atoms: &HashSet<usize>,
    covalent_links_at: &HashMap<usize, usize>,
    link: &Link,
) -> Result<BondOrder, LinkIssueKind> {
    for atom in [link.atom_1, link.atom_2] {
        if backbone_atoms.contains(&atom) {
            return Err(LinkIssueKind::SharesBackbone { atom });
        }
        let links = covalent_links_at[&atom];
        if links > 1 {
            return Err(LinkIssueKind::SeveralLinks { atom, links });
        }
    }
    // Past these checks, every missing partner at either end belongs to this link.
    let ends = [link_end(structure, templates, link.atom_1)?, link_end(structure, templates, link.atom_2)?];
    let order = replaced_order(&ends)?;
    // An end with no missing heavy partner must be missing hydrogens instead.
    for end in ends.iter().filter(|end| end.missing_heavy.is_empty()) {
        if end.missing_hydrogens < order.as_u8() as usize {
            return Err(LinkIssueKind::NotEnoughHydrogens { atom: end.atom, needed: order.as_u8(), available: end.missing_hydrogens });
        }
        if let Some(partner) = &end.kept_multiple_bond {
            return Err(LinkIssueKind::AdditionAcrossMultipleBond { atom: end.atom, partner: partner.clone() });
        }
    }
    Ok(order)
}

/// One end of a link, compared with its residue's template.
struct LinkEnd {
    atom: usize,
    /// Orders of the template bonds to heavy partners missing from the structure.
    missing_heavy: Vec<BondOrder>,
    /// Template hydrogens missing from the structure.
    missing_hydrogens: usize,
    /// A present partner the atom keeps a non-aromatic double or triple bond to.
    kept_multiple_bond: Option<String>,
}

fn link_end(structure: &Structure, templates: &HashMap<String, Vec<ResidueTemplateBond>>, atom: usize) -> Result<LinkEnd, LinkIssueKind> {
    let residue = &structure.residues[structure.residue_of(atom)];
    let partners: Vec<_> = template_partners(templates, residue, &structure.atoms[atom].name).collect();
    if partners.is_empty() {
        return Err(LinkIssueKind::NoTemplate { atom });
    }
    let mut end = LinkEnd { atom, missing_heavy: vec![], missing_hydrogens: 0, kept_multiple_bond: None };
    for (bond, partner, element) in partners {
        match (atom_named(structure, residue, partner).is_some(), element) {
            (false, Element::H) => end.missing_hydrogens += 1,
            (false, _) => end.missing_heavy.push(bond.order),
            (true, _) if bond.order != BondOrder::Single && !bond.aromatic => {
                end.kept_multiple_bond.get_or_insert_with(|| partner.to_string());
            }
            (true, _) => {}
        }
    }
    Ok(end)
}

// The order of the template bond the link took over, or Single when neither end is missing a heavy partner.
fn replaced_order([end_1, end_2]: &[LinkEnd; 2]) -> Result<BondOrder, LinkIssueKind> {
    match (end_1.missing_heavy.as_slice(), end_2.missing_heavy.as_slice()) {
        ([], []) => Ok(BondOrder::Single),
        ([order], []) | ([], [order]) => Ok(*order),
        ([order_1], [order_2]) if order_1 == order_2 => Ok(*order_1),
        ([order_1], [order_2]) => Err(LinkIssueKind::MismatchedMissingHeavy { orders: (*order_1, *order_2) }),
        _ => {
            let end = if end_1.missing_heavy.len() > 1 { end_1 } else { end_2 };
            Err(LinkIssueKind::SeveralMissingHeavy { atom: end.atom, orders: end.missing_heavy.clone() })
        }
    }
}

// The template bonds of the atom with this name, each with the partner's name and element.
fn template_partners<'a>(
    templates: &'a HashMap<String, Vec<ResidueTemplateBond>>,
    residue: &Residue,
    name: &'a str,
) -> impl Iterator<Item = (&'a ResidueTemplateBond, &'a str, Element)> {
    templates.get(&residue.name).into_iter().flatten().filter_map(move |bond| {
        if bond.atom_name_1 == name {
            Some((bond, bond.atom_name_2.as_str(), bond.element_2))
        } else if bond.atom_name_2 == name {
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
