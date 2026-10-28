use std::{collections::HashMap, fs, hash::Hash};

use prote::{
    bond::Bond,
    conformer::Conformer,
    contact::Contact,
    mmcif::MmCIF,
    selection::Selection,
    spatial_index::IndexedStructure,
    structure::{self, EntityKind, Residue},
    Cif,
};

fn main() {
    let cif = Cif::new(&fs::read_to_string("./cifs/3PTB.cif").unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    let indexed_structure = IndexedStructure::new(Conformer::new(&structure, &frames[0]).unwrap());
    println!("{:?}", structure.bonds)
    // let protein = Selection::by_instance(&structure, |instance, entity| entity.kind == EntityKind::Polymer);
    // let ligands = Selection::by_instance(&structure, |instance, entity| entity.kind == EntityKind::NonPolymer);
    // let contacts = Contact::contacts(&indexed_structure, &protein, &ligands);
    // for contact in contacts {
    //     let atom_1 = &structure.atoms[contact.atom_1];
    //     let resiude = structure.residues.iter().find(|residue| residue.atoms.contains(&contact.atom_1)).unwrap();
    //     let atom_2 = &structure.atoms[contact.atom_2];
    //     println!(
    //         "Distance from protein residue {:?} atom {:?} to ligand atom {:?} is {:?}",
    //         resiude.name, atom_1.name, atom_2.name, contact.distance
    //     );
    // }
}
