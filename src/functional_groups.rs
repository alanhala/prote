use crate::{
    bond,
    bond_order::BondOrder,
    bonding::{covalent_heavy_partners, Bonding, Partner},
    element::Element,
    structure::Structure,
};

#[derive(Debug)]
pub enum FunctionalGroup {
    Carboxyl { carbon: usize, double_oxygen: usize, single_oxygen: usize },
    Sulfonate { sulfur: usize, double_oxygens: [usize; 2], single_oxygen: usize },
    InorganicPhosphate { phosphate: usize, double_oxygen: usize, single_oxygens: [usize; 3] },
    PhosphateMonoester { phosphate: usize, double_oxygen: usize, single_oxygens: [usize; 2] },
    PhosphateDiester { phosphate: usize, double_oxygen: usize, single_oxygen: usize },
    AromaticThiol { sulfur: usize, carbon: usize },
    Thiol { sulfur: usize, carbon: usize },
    Guanidine { carbon: usize, double_nitrogen: usize, single_nitrogens: [usize; 2] },
    Amidine { carbon: usize, double_nitrogen: usize, single_nitrogen: usize },
    Imidazole { carbon: usize, double_nitrogen: usize, single_nitrogen: usize },
    Pyridine { nitrogen: usize },
    Amine { nitrogen: usize },
    Phenol { oxygen: usize, carbon: usize },
}

pub fn functional_groups(structure: &Structure, bonding: &Bonding) -> Vec<FunctionalGroup> {
    let mut functional_groups = vec![];
    for atom in 0..structure.atoms.iter().len() {
        match structure.atoms[atom].element {
            Element::C => {
                if let Some(carboxyl) = carboxyl(structure, bonding, atom) {
                    functional_groups.push(carboxyl);
                } else if let Some(amidine) = amidine(structure, bonding, atom) {
                    functional_groups.push(amidine);
                }
            }
            Element::S => {
                if let Some(sulfonate) = sulfonate(structure, bonding, atom) {
                    functional_groups.push(sulfonate);
                } else if let Some(thiol) = thiol(structure, bonding, atom) {
                    functional_groups.push(thiol);
                }
            }
            Element::N => {
                if let Some(imidazole) = imidazole(structure, bonding, atom) {
                    functional_groups.push(imidazole);
                } else if let Some(pyridine) = pyridine(structure, bonding, atom) {
                    functional_groups.push(pyridine);
                } else if let Some(amine) = amine(structure, bonding, atom) {
                    functional_groups.push(amine);
                }
            }
            Element::P => {
                if let Some(phosphate) = phosphate(structure, bonding, atom) {
                    functional_groups.push(phosphate);
                }
            }
            Element::O => {
                if let Some(phenol) = phenol(structure, bonding, atom) {
                    functional_groups.push(phenol);
                }
            }
            _ => {}
        }
    }
    functional_groups
}

fn phenol(structure: &Structure, bonding: &Bonding, oxygen: usize) -> Option<FunctionalGroup> {
    let [partner]: [Partner; 1] = covalent_heavy_partners(structure, bonding, oxygen).collect::<Vec<Partner>>().try_into().ok()?;
    let carbon = partner.atom?;
    if partner.element != Element::C || partner.order != BondOrder::Single {
        return None;
    }
    covalent_heavy_partners(structure, bonding, carbon).any(|partner| partner.aromatic).then_some(FunctionalGroup::Phenol { oxygen, carbon })
}

fn imidazole(structure: &Structure, bonding: &Bonding, double_nitrogen: usize) -> Option<FunctionalGroup> {
    let ring = bonding.rings.iter().find(|ring| ring.aromatic && ring.atoms.len() == 5 && ring.atoms.contains(&double_nitrogen))?;
    let [single_nitrogen]: [usize; 1] = ring
        .atoms
        .iter()
        .copied()
        .filter(|&atom| atom != double_nitrogen && structure.atoms[atom].element == Element::N)
        .collect::<Vec<usize>>()
        .try_into()
        .ok()?;
    let double_nitrogen_partners: [Partner; 2] =
        covalent_heavy_partners(structure, bonding, double_nitrogen).collect::<Vec<Partner>>().try_into().ok()?;
    let carbon = double_nitrogen_partners.iter().find(|partner| partner.element == Element::C && partner.order == BondOrder::Double)?.atom?;
    covalent_heavy_partners(structure, bonding, carbon)
        .any(|partner| partner.atom == Some(single_nitrogen) && partner.order == BondOrder::Single)
        .then_some(FunctionalGroup::Imidazole { carbon, double_nitrogen, single_nitrogen })
}

fn pyridine(structure: &Structure, bonding: &Bonding, nitrogen: usize) -> Option<FunctionalGroup> {
    bonding.rings.iter().find(|ring| {
        ring.aromatic
            && ring.atoms.len() == 6
            && ring.atoms.contains(&nitrogen)
            && ring.atoms.iter().all(|&atom| atom == nitrogen || structure.atoms[atom].element != Element::N)
    })?;
    let partners: [Partner; 2] = covalent_heavy_partners(structure, bonding, nitrogen).collect::<Vec<Partner>>().try_into().ok()?;
    partners.iter().any(|partner| partner.order == BondOrder::Double).then_some(FunctionalGroup::Pyridine { nitrogen })
}

fn amine(structure: &Structure, bonding: &Bonding, nitrogen: usize) -> Option<FunctionalGroup> {
    let partners: Vec<Partner> = covalent_heavy_partners(structure, bonding, nitrogen).collect();
    if partners.is_empty() || partners.len() > 3 {
        return None;
    }
    partners
        .iter()
        .all(|partner| {
            partner.element == Element::C
                && partner.order == BondOrder::Single
                && partner.atom.is_some_and(|carbon| {
                    covalent_heavy_partners(structure, bonding, carbon).all(|partner| partner.order == BondOrder::Single && !partner.aromatic)
                })
        })
        .then_some(FunctionalGroup::Amine { nitrogen })
}

fn amidine(structure: &Structure, bonding: &Bonding, carbon: usize) -> Option<FunctionalGroup> {
    let nitrogens: Vec<_> = covalent_heavy_partners(structure, bonding, carbon).filter(|partner| partner.element == Element::N).collect();
    let [double_nitrogen]: [usize; 1] =
        nitrogens.iter().filter(|n| n.order == BondOrder::Double).map(|n| n.atom).collect::<Option<Vec<usize>>>()?.try_into().ok()?;
    let single_nitrogens = nitrogens.iter().filter(|n| n.order == BondOrder::Single).map(|n| n.atom).collect::<Option<Vec<usize>>>()?;
    // Nitrogens with another partner are ring members, as in histidine, where the pK is far lower.
    let free = |nitrogen: usize| covalent_heavy_partners(structure, bonding, nitrogen).count() == 1;
    match single_nitrogens[..] {
        [single_nitrogen] if free(double_nitrogen) && free(single_nitrogen) => {
            Some(FunctionalGroup::Amidine { carbon, double_nitrogen, single_nitrogen })
        }
        [n1, n2] if [double_nitrogen, n1, n2].into_iter().filter(|&nitrogen| free(nitrogen)).count() >= 2 => {
            Some(FunctionalGroup::Guanidine { carbon, double_nitrogen, single_nitrogens: [n1, n2] })
        }
        _ => None,
    }
}

fn thiol(structure: &Structure, bonding: &Bonding, sulfur: usize) -> Option<FunctionalGroup> {
    let [partner]: [Partner; 1] = covalent_heavy_partners(structure, bonding, sulfur).collect::<Vec<Partner>>().try_into().ok()?;
    let carbon = partner.atom?;
    if partner.element != Element::C || partner.order != BondOrder::Single {
        return None;
    }
    if covalent_heavy_partners(structure, bonding, carbon).any(|partner| partner.aromatic) {
        return Some(FunctionalGroup::AromaticThiol { sulfur, carbon });
    }
    Some(FunctionalGroup::Thiol { sulfur, carbon })
}

fn phosphate(structure: &Structure, bonding: &Bonding, phosphate: usize) -> Option<FunctionalGroup> {
    let oxygens: Vec<_> = covalent_heavy_partners(structure, bonding, phosphate).filter(|partner| partner.element == Element::O).collect();
    let [double_oxygen]: [usize; 1] =
        oxygens.iter().filter(|o| o.order == BondOrder::Double).filter_map(|o| o.atom).collect::<Vec<usize>>().try_into().ok()?;
    let single_oxygens: Vec<usize> = oxygens
        .iter()
        .filter(|o| o.order == BondOrder::Single)
        .filter_map(|o| o.atom)
        .filter(|&atom| covalent_heavy_partners(structure, bonding, atom).count() == 1)
        .collect::<Vec<usize>>();
    match single_oxygens[..] {
        [single_oxygen] => Some(FunctionalGroup::PhosphateDiester { phosphate, double_oxygen, single_oxygen }),
        [o1, o2] => Some(FunctionalGroup::PhosphateMonoester { phosphate, double_oxygen, single_oxygens: [o1, o2] }),
        [o1, o2, o3] => Some(FunctionalGroup::InorganicPhosphate { phosphate, double_oxygen, single_oxygens: [o1, o2, o3] }),
        _ => None,
    }
}

fn sulfonate(structure: &Structure, bonding: &Bonding, sulfur: usize) -> Option<FunctionalGroup> {
    let oxygens: Vec<_> = covalent_heavy_partners(structure, bonding, sulfur).filter(|partner| partner.element == Element::O).collect();
    let double_oxygens: [usize; 2] =
        oxygens.iter().filter(|o| o.order == BondOrder::Double).map(|o| o.atom).collect::<Option<Vec<usize>>>()?.try_into().ok()?;
    // Two free single-bonded oxygens is free sulfate, which is not a sulfonate.
    let [single_oxygen]: [usize; 1] = oxygens
        .iter()
        .filter(|o| o.order == BondOrder::Single)
        .filter_map(|o| o.atom)
        .filter(|&atom| covalent_heavy_partners(structure, bonding, atom).count() == 1)
        .collect::<Vec<usize>>()
        .try_into()
        .ok()?;
    Some(FunctionalGroup::Sulfonate { sulfur, single_oxygen, double_oxygens })
}

fn carboxyl(structure: &Structure, bonding: &Bonding, carbon: usize) -> Option<FunctionalGroup> {
    let oxygens: Vec<_> = covalent_heavy_partners(structure, bonding, carbon).filter(|partner| partner.element == Element::O).collect();
    let double_oxygen = oxygens.iter().filter(|o| o.order == BondOrder::Double).find_map(|o| o.atom)?;
    let single_oxygen = oxygens
        .iter()
        .filter(|o| o.order == BondOrder::Single)
        .filter_map(|o| o.atom)
        .find(|&atom| covalent_heavy_partners(structure, bonding, atom).count() == 1)?;
    Some(FunctionalGroup::Carboxyl { carbon, single_oxygen, double_oxygen })
}
