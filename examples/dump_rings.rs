//! Prints the rings found in every non-water molecule of a CIF, with atom names, coplanarity and mean bond ratio.
use std::fs;
use prote::{bond_graph::BondGraph, mmcif::MmCIF, spatial_index::SpatialIndex, Cif};
fn main() {
    let mut args = std::env::args().skip(1);
    let path = args.next().unwrap();
    let only: Option<String> = args.next();
    let cif = Cif::new(&fs::read_to_string(&path).unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    for molecule in mmcif.build_molecules() {
        if molecule.name == "water" { continue; }
        let comp = &molecule.topology.components[0].name;
        if only.as_deref().map_or(false, |o| o != comp) { continue; }
        let conformer = molecule.conformer(0).unwrap();
        let spatial_index = SpatialIndex::new(&conformer.frame.positions);
        let graph = BondGraph::new(&conformer, &spatial_index, &conformer);
        println!("{} ({} atoms, {} rings)", comp, molecule.topology.atoms.len(), graph.rings.len());
        for ring in &graph.rings {
            let names: Vec<&str> = ring.atom_pointers.iter().map(|p| molecule.topology.atoms[p.index()].name.as_str()).collect();
            let n = ring.atom_pointers.len();
            let mean: f64 = (0..n).map(|k| {
                let (i, j) = (ring.atom_pointers[k].index(), ring.atom_pointers[(k + 1) % n].index());
                let bond = graph.bond_between(i, j).unwrap();
                bond.distance / (molecule.topology.atoms[i].covalent_radius() + molecule.topology.atoms[j].covalent_radius())
            }).sum::<f64>() / n as f64;
            println!("  size {n:2} coplanar={} mean_ratio={mean:.3} {}", ring.coplanar, names.join("-"));
        }
        break;
    }
}
