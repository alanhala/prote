use std::fs;

use prote::{ensemble::Ensemble, mmcif::MmCIF, Cif};

fn main() {
    let cif_file = fs::read_to_string("./cifs/3PTB.cif").unwrap();
    let cif = Cif::new(&cif_file);
    let mmcif = MmCIF::new(cif).unwrap();
    let ensembles = mmcif.build_ensembles();
    let proteins: Vec<&Ensemble> = ensembles
        .iter()
        .filter(|ensemble| ensemble.topology.name == "BETA-TRYPSIN")
        .collect();
    for protein in &proteins {
        println!("{:?}", protein.topology);
    }
}
