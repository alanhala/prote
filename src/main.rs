use std::fs;

use prote::{mmcif::MmCIF, Cif};

fn main() {
    let cif_file = fs::read_to_string("./cifs/4d1m.cif").unwrap();
    let cif = Cif::new(&cif_file);
    let mmcif = MmCIF::new(cif).unwrap();
    println!("{:#?}", mmcif.build_topologies());
}
