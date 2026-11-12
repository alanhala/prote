use crate::{bond::BondKind, structure::Structure};

pub fn formal_charges(structure: &Structure) -> Vec<Option<i8>> {
    let mut charges = vec![None; structure.atoms.len()];
    for (id, atom) in structure.atoms.iter().enumerate() {
        let bonds: u8 = structure
            .bonding
            .bonds
            .neighbors(id)
            .iter()
            .map(|&(_, edge)| match structure.bonding.bonds.edges()[edge].kind {
                BondKind::Coordination => 0,
                BondKind::Covalent { order, .. } => order.as_u8(),
            })
            .sum();
        let missing: u8 = structure.bonding.missing_partners[id].iter().map(|partner| partner.order.as_u8()).sum();
        let value = bonds + structure.bonding.implicit_hydrogens[id] + missing;
        let value = value as i8;
        charges[id] = atom.element.valences().iter().map(|&valence| value - valence as i8).min_by_key(|charge| charge.abs());
    }
    charges
}
