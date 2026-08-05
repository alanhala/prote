use std::fs;
use std::path::Path;

use prote::{
    bond,
    bond_graph::BondGraph,
    cif,
    ensemble::{resolve, Ensemble},
    mmcif::MmCIF,
    molecule,
    spatial_index::SpatialIndex,
    Cif,
};

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
    let mut bond_graph = BondGraph::new(&protein, &spatial_index, &protein);
    bond_graph.assign_orders(&protein.topology, &[&mmcif]);
    bond_graph.assign_remaining_orders(&ensembles);
    for atom_bonds in bond_graph.bonds.iter() {
        for bond in atom_bonds.iter() {
            if let Some(order) = bond.order {
                let (atom_1, _) = resolve(&ensembles, &bond.atom_1);
                let (atom_2, _) = resolve(&ensembles, &bond.atom_2);
                println!(
                    "{:?}{} -- {:?}{}: {:?}",
                    atom_1.element(),
                    atom_1.name,
                    atom_2.element(),
                    atom_2.name,
                    order
                );
            }
        }
    }
    // println!("{:#?}", bond_graph);
    // let residue_index = pocket::residue_index_by_auth_seq_id(&mmcif, "A");
    // let pockets = pocket::load_pockets(Path::new("./cifs/3PTB_out/pockets"), "A", &residue_index);

    // for (i, pocket) in pockets.iter().enumerate() {
    //     let names: Vec<&str> = pocket
    //         .residues
    //         .iter()
    //         .map(|&residue_index| protein.topology.residues[residue_index].name.as_str())
    //         .collect();
    //     println!("Pocket {}: {} residues -> {:?}", i + 1, pocket.residues.len(), names);
    // }

    // println!("{:?}", );
}
