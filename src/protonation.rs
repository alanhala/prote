use crate::{
    bonding::{covalent_partners, Bonding},
    functional_groups::FunctionalGroup,
    structure::{EntityKind, PolymerKind, Structure},
};

/// Measured in water for the compound named, from R. Williams's pKa compilation; sulfonate is from Dimorphite-DL.
pub fn group_pks(group: &FunctionalGroup) -> &'static [f64] {
    match group {
        FunctionalGroup::Carboxyl { .. } => &[4.76],
        FunctionalGroup::Sulfonate { .. } => &[-2.4],
        FunctionalGroup::InorganicPhosphate { .. } => &[1.97, 6.82, 12.5],
        FunctionalGroup::PhosphateMonoester { .. } => &[1.54, 6.31],
        FunctionalGroup::PhosphateDiester { .. } => &[1.72],
        FunctionalGroup::AromaticThiol { .. } => &[6.52],
        FunctionalGroup::Thiol { .. } => &[10.50],
        FunctionalGroup::Guanidine { .. } => &[13.71],
        FunctionalGroup::Amidine { .. } => &[11.6],
        FunctionalGroup::Imidazole { .. } => &[6.95],
        FunctionalGroup::Pyridine { .. } => &[5.17],
        FunctionalGroup::Amine { .. } => &[10.62],
        FunctionalGroup::Phenol { .. } => &[9.95],
    }
}

fn in_peptide(structure: &Structure, atom: usize) -> bool {
    let entity = structure.instances[structure.instance_of(structure.residue_of(atom))].entity;
    matches!(structure.entities[entity].kind, EntityKind::Polymer(PolymerKind::PolypeptideD | PolymerKind::PolypeptideL))
}

/// The pK measured on the residue itself, where the atom is a standard residue's own.
fn residue_pk(structure: &Structure, atom: usize) -> Option<f64> {
    if !in_peptide(structure, atom) {
        return None;
    }
    let residue = &structure.residues[structure.residue_of(atom)];
    match (residue.name.as_str(), structure.atoms[atom].name.as_str()) {
        (_, "N") if residue.seq_id == Some(1) => Some(8.0),
        (_, "OXT") => Some(3.67),
        ("ASP", "OD1" | "OD2") => Some(3.67),
        ("GLU", "OE1" | "OE2") => Some(4.25),
        ("HIS", "ND1" | "NE2") => Some(6.54),
        ("CYS", "SG") => Some(8.55),
        ("LYS", "NZ") => Some(10.40),
        ("ARG", "NE" | "NH1" | "NH2") => Some(12.0),
        ("TYR", "OH") => Some(9.84),
        _ => None,
    }
}

fn set_hydrogens(structure: &Structure, bonding: &mut Bonding, atom: usize, change: i8) {
    let used: i8 = covalent_partners(structure, bonding, atom).map(|partner| partner.order.as_u8() as i8).sum();
    let hydrogens = structure.atoms[atom].element.valences().iter().map(|&valence| valence as i8 + change - used).min().unwrap();
    bonding.implicit_hydrogens[atom] = hydrogens.max(0) as u8;
}

pub fn protonate(structure: &Structure, bonding: &mut Bonding, groups: &[FunctionalGroup], ph: f64) {
    groups.iter().for_each(|group| {
        // The atoms that change, one per pK of the group, and the charge each takes:
        // +1 gains an H below its pK, -1 loses its H from its pK up.
        let (atoms, charge): (Vec<usize>, i8) = match *group {
            FunctionalGroup::Carboxyl { single_oxygen, .. } => (vec![single_oxygen], -1),
            FunctionalGroup::Sulfonate { single_oxygen, .. } => (vec![single_oxygen], -1),
            FunctionalGroup::InorganicPhosphate { single_oxygens, .. } => (single_oxygens.to_vec(), -1),
            FunctionalGroup::PhosphateMonoester { single_oxygens, .. } => (single_oxygens.to_vec(), -1),
            FunctionalGroup::PhosphateDiester { single_oxygen, .. } => (vec![single_oxygen], -1),
            FunctionalGroup::AromaticThiol { sulfur, .. } => (vec![sulfur], -1),
            FunctionalGroup::Thiol { sulfur, .. } => (vec![sulfur], -1),
            FunctionalGroup::Guanidine { double_nitrogen, .. } => (vec![double_nitrogen], 1),
            FunctionalGroup::Amidine { double_nitrogen, .. } => (vec![double_nitrogen], 1),
            FunctionalGroup::Imidazole { double_nitrogen, .. } => (vec![double_nitrogen], 1),
            FunctionalGroup::Pyridine { nitrogen } => (vec![nitrogen], 1),
            FunctionalGroup::Amine { nitrogen } => (vec![nitrogen], 1),
            FunctionalGroup::Phenol { oxygen, .. } => (vec![oxygen], -1),
        };
        for (&atom, &group_pk) in atoms.iter().zip(group_pks(group)) {
            let pk = residue_pk(structure, atom).unwrap_or(group_pk);
            let charged = if charge > 0 { ph < pk } else { ph >= pk };
            set_hydrogens(structure, bonding, atom, if charged { charge } else { 0 });
        }
    });
}
