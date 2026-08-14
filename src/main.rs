use std::fs;

use prote::{
    atom::Atom,
    bond_graph::{self, BondGraph, Hybridization},
    mmcif::MmCIF,
    spatial_index::SpatialIndex,
    Cif,
};

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
    for (i, atom) in protein.topology.atoms.iter().enumerate() {
        match bond_graph.hybridization_at(i) {
            Some(Hybridization::SP2) => println!("Found SP2 {:?}", atom.name),
            _ => {}
        }
    }
    bond_graph.dfs();
}

// fn visit(node: usize, neighbors: &Vec<Vec<usize>>, visited: &mut Vec<usize>, cycles: &mut Vec<Vec<usize>>) {
//     for &neighbor in &neighbors[node] {
//         if visited.last() == Some(&neighbor) {
//             continue;
//         } else if let Some(idx) = visited.iter().position(|&x| x == neighbor) {
//             let mut cycle = visited[idx..].to_vec();
//             cycle.push(node);
//             cycles.push(cycle);
//         } else {
//             visited.push(node);
//             visit(neighbor, neighbors, visited, cycles);
//             visited.pop();
//         }
//     }
// }

// fn main() {
//     let neighbors = vec![vec![1, 2], vec![0, 2], vec![0, 1]];
//     let mut visited = vec![];
//     let mut cycles: Vec<Vec<usize>> = vec![];
//     visit(0, &neighbors, &mut visited, &mut cycles);
//     println!("finished: {} cycles found", cycles.len());
// }
