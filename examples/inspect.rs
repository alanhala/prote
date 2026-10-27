//! Prints geometry and hybridization around atoms: cargo run --example inspect -- file.cif COMP ATOM [min_bond_length_to_flag]
use std::fs;
use prote::{bond_graph::BondGraph, mmcif::MmCIF, spatial_index::SpatialIndex, Cif};
fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    let (path, comp, name) = (&args[0], &args[1], &args[2]);
    let flag: f64 = args.get(3).and_then(|s| s.parse().ok()).unwrap_or(0.0);
    let cif = Cif::new(&fs::read_to_string(path).unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    let proteins = mmcif.build_proteins();
    let molecules = mmcif.build_molecules();
    let conformers = proteins.iter().map(|p| p.conformer(0).unwrap()).chain(molecules.iter().map(|m| m.conformer(0).unwrap()));
    for conformer in conformers {
        let spatial_index = SpatialIndex::new(&conformer.frame.positions);
        let mut graph = BondGraph::new(&conformer, &spatial_index, &conformer);
        graph.perceive_orders();
        let topology = conformer.topology;
        for (i, atom) in topology.atoms.iter().enumerate() {
            let component = &topology.components[topology.component_index_for(i)];
            if component.name != *comp || atom.name != *name { continue; }
            let longest = graph.bonds[i].iter().map(|b| b.distance).fold(0.0, f64::max);
            if longest < flag { continue; }
            println!("{} {} (graph index {i}) hybridization {:?} rings {:?}", component.name, atom.name, graph.hybridization_at(i), graph.atom_rings[i]);
            for bond in &graph.bonds[i] {
                let other = &topology.atoms[bond.atom_2.index()];
                let ratio = bond.distance / (atom.covalent_radius() + other.covalent_radius());
                println!("   - {:4} {:.2} A ratio {:.3} {:?}", other.name, bond.distance, ratio, bond.order);
            }
            for angle in &graph.angles[i] {
                println!("   angle {}-{}-{} {:.1}", topology.atoms[angle.atom_1.index()].name, atom.name, topology.atoms[angle.atom_2.index()].name, angle.value);
            }
        }
    }
}
