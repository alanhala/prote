use std::fs;

use prote::{bond_graph::BondGraph, mmcif::MmCIF, spatial_index::SpatialIndex, Cif};

fn main() {
    let cif = Cif::new(&fs::read_to_string("./cifs/3PTB.cif").unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    let proteins = mmcif.build_proteins();
    let molecules = mmcif.build_molecules();
    for protein in &proteins {
        println!("Protein: {:?}", protein.name)
    }
    for molecule in &molecules {
        println!("Molecule {:?}", molecule.name)
    }
    let protein = proteins.first().unwrap();
    let protein_conformer = protein.conformer(0).unwrap();
    let spatial_index = SpatialIndex::new(&protein_conformer.frame.positions);
    let bond_graph = BondGraph::new(&protein_conformer, &spatial_index, &protein_conformer);
    println!("{:?}", bond_graph)
}
