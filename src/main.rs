use std::fs;

use prote::{bond_graph::BondGraph, ensemble::Ensemble, mmcif::MmCIF, spatial_index::SpatialIndex, Cif};

fn main() {
    let cif_file = fs::read_to_string("./cifs/3PTB.cif").unwrap();
    let cif = Cif::new(&cif_file);
    let mmcif = MmCIF::new(cif).unwrap();
    let ensembles = mmcif.build_ensembles();

    let proteins: Vec<&Ensemble> = ensembles
        .iter()
        .filter(|ensemble| ensemble.topology.name == "BETA-TRYPSIN")
        .collect();
    let protein = proteins.get(0).unwrap().molecule(0);
    let spatial_index = SpatialIndex::new(&protein.conformer.positions);
    let bond_graph = BondGraph::new(&protein, &spatial_index, &protein);
    println!(
        "{} atoms, {} bonds",
        protein.topology.atoms.len(),
        bond_graph.bonds.iter().map(Vec::len).sum::<usize>()
    );
}
