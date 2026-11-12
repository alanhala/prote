use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    fs,
    path::Path,
};

use prote::{
    bond::BondKind,
    charge,
    contact::{Contact, ContactType},
    dot::{self, edge},
    mmcif::MmCIF,
    selection::Selection,
    spatial_index::SpatialIndex,
    structure::{EntityKind, PolymerKind},
    Cif,
};

fn main() {
    let cif = Cif::new(&fs::read_to_string("./cifs/4HHB.cif").unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    let hem = Selection::by_residue(&structure, |residue| residue.name == "HEM");
    let mut nodes = vec![];
    let mut edges = vec![];
    for bond in structure.bonding.bonds.edges() {
        if hem.contains(bond.atom_1) || hem.contains(bond.atom_2) {
            nodes.push(dot::Node { id: bond.atom_1.to_string(), label: Some(structure.atoms[bond.atom_1].name.to_string()), ..Default::default() });
            nodes.push(dot::Node { id: bond.atom_2.to_string(), label: Some(structure.atoms[bond.atom_2].name.to_string()), ..Default::default() });
            match bond.kind {
                BondKind::Coordination => edges.push(dot::Edge {
                    from: bond.atom_1.to_string(),
                    to: bond.atom_2.to_string(),
                    color: Some("red".to_string()),
                    dashed: true,
                    ..Default::default()
                }),
                BondKind::Covalent { order, .. } => {
                    edges.push(dot::Edge { from: bond.atom_1.to_string(), to: bond.atom_2.to_string(), lines: order.as_u8(), ..Default::default() })
                }
            }
        }
    }
    let graph = dot::Graph { title: "Test".to_string(), nodes, edges };
    dot::write_svg(&graph, dot::Layout::Sfdp, Path::new("test.svg")).unwrap();
    charge::formal_charges(&structure).iter().enumerate().for_each(|(atom_id, charge)| {
        if charge.is_some_and(|charge| charge != 0) {
            let residue = &structure.residues[structure.residue_of(atom_id)];
            let bonds = structure.bonding.bonds.neighbors(atom_id).len();
            println!(
                "atom {:?} {:?} {:?} {:?} {:?}",
                residue, structure.atoms[atom_id], charge, structure.bonding.implicit_hydrogens[atom_id], bonds
            );
        }
    });
}
