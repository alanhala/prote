use std::{collections::HashMap, fs, hash::Hash};

use prote::{
    mmcif::MmCIF,
    selection::Selection,
    structure::{self, EntityKind, Residue},
    Cif,
};

fn main() {
    let cif = Cif::new(&fs::read_to_string("./cifs/3PTB.cif").unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    let protein = Selection::by_instance(&structure, |instance, entity| entity.kind == EntityKind::Polymer);
    let tyr = Selection::by_residue(&structure, |residue| residue.name == "TYR");
    let lys = Selection::by_residue(&structure, |residue| residue.name == "LYS");
    let not_c = Selection::by_atom(&structure, |atom| atom.name == "C").not();
    let selection = protein.and(tyr.or(lys)).and(not_c);
    for &ri in &selection.residues_touched(&structure) {
        let residue = &structure.residues[ri];
        println!("Resnum: {:?}, Resname: {:?}", residue.number, residue.name);
        for i in residue.atoms.clone() {
            if selection.atoms[i] {
                println!("Atom: {:?}", structure.atoms[i].name);
            }
        }
    }
}
