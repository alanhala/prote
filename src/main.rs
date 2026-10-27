use std::{collections::HashMap, fs, hash::Hash};

use prote::{mmcif::MmCIF, Cif};

fn main() {
    let cif = Cif::new(&fs::read_to_string("./cifs/3PTB.cif").unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    println!("{:?}", structure);
    // for ring_pointers in bond_graph.atom_rings.get(1561) {
    //     for ring_pointer in ring_pointers {
    //         println!("{:?}", bond_graph.rings.get(*ring_pointer).unwrap().atom_pointers.len());
    //         for atom_pointer in &bond_graph.rings.get(*ring_pointer).unwrap().atom_pointers {
    //             println!("{:?}", atom_pointer.index());
    //         }
    //     }
    // }
}
