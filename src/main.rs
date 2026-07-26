use std::fs;

use prote::{mmcif::MmCIF, molecule::Molecule, Cif};

fn main() {
    let cif_file = fs::read_to_string("./cifs/3PTB.cif").unwrap();
    let cif = Cif::new(&cif_file);
    let mmcif = MmCIF::new(cif).unwrap();
    let molecules = mmcif.build_molecules();
    let proteins: Vec<&Molecule> = molecules
        .iter()
        .filter(|molecule| molecule.topology.name == "BETA-TRYPSIN")
        .collect();
    for protein in &proteins {
        println!("{:?}", protein.topology);
    }
}
