//! Compares geometric bond-order perception against each CIF's `_chem_comp_bond` dictionary.
//! Usage: cargo run --example validate_orders -- cifs/3PTB.cif [more.cif ...]

use std::{collections::BTreeMap, fs};

use prote::{bond_graph::BondGraph, bond_order::BondOrder, mmcif::MmCIF, spatial_index::SpatialIndex, Cif};

fn main() {
    let mut total = 0usize;
    let mut agreed = 0usize;
    let mut swaps = 0usize;
    let mut mismatches: BTreeMap<String, usize> = BTreeMap::new();
    let mut inter: BTreeMap<String, usize> = BTreeMap::new();
    for path in std::env::args().skip(1) {
        let cif = Cif::new(&fs::read_to_string(&path).unwrap());
        let mmcif = MmCIF::new(cif).unwrap();
        let proteins = mmcif.build_proteins();
        let molecules = mmcif.build_molecules();
        let conformers = proteins
            .iter()
            .map(|protein| protein.conformer(0).unwrap())
            .chain(molecules.iter().filter(|m| m.name != "water").map(|molecule| molecule.conformer(0).unwrap()));
        for conformer in conformers {
            let spatial_index = SpatialIndex::new(&conformer.frame.positions);
            let mut graph = BondGraph::new(&conformer, &spatial_index, &conformer);
            graph.perceive_orders();
            let topology = conformer.topology;
            // Per (atom, neighbor element): the expected and perceived orders of those bonds. A
            // mismatch whose orders reappear, permuted, among same-element neighbors of one of its
            // atoms is a swap between resonance-equivalent partners (carboxylate O's, guanidinium
            // N's, phosphate O's), not a wrong perception.
            let mut around: BTreeMap<(usize, u8), (Vec<u8>, Vec<u8>)> = BTreeMap::new();
            let mut pending: Vec<((usize, u8), (usize, u8), String)> = vec![];
            let value = |order: BondOrder| match order {
                BondOrder::Aromatic => 9,
                other => other.as_u8(),
            };
            for bonds in &graph.bonds {
                for bond in bonds {
                    let (i, j) = (bond.atom_1.index(), bond.atom_2.index());
                    if i >= j {
                        continue;
                    }
                    let (a, b) = (&topology.atoms[i], &topology.atoms[j]);
                    let (ci, cj) = (topology.component_index_for(i), topology.component_index_for(j));
                    let got = bond.order.unwrap();
                    if ci != cj {
                        *inter.entry(format!("{}-{} {:?}", a.name, b.name, got)).or_default() += 1;
                        continue;
                    }
                    let comp = &topology.components[ci].name;
                    let key = (comp.clone(), a.name.clone(), b.name.clone());
                    let rev = (comp.clone(), b.name.clone(), a.name.clone());
                    let Some(entry) = mmcif.chem_comp_bonds.get(&key).or_else(|| mmcif.chem_comp_bonds.get(&rev)) else {
                        *mismatches.entry(format!("{comp} {}-{}: not in dictionary, got {got:?}", a.name, b.name)).or_default() += 1;
                        continue;
                    };
                    let expected = if entry.aromatic_flag { BondOrder::Aromatic } else { entry.value_order };
                    total += 1;
                    let key_i = (i, b.atomic_number());
                    let key_j = (j, a.atomic_number());
                    for key in [key_i, key_j] {
                        let entry = around.entry(key).or_default();
                        entry.0.push(value(expected));
                        entry.1.push(value(got));
                    }
                    if expected == got {
                        agreed += 1;
                    } else {
                        pending.push((key_i, key_j, format!("{comp} {}-{}: expected {expected:?}, got {got:?} ({:.2} A)", a.name, b.name, bond.distance)));
                    }
                }
            }
            for (key_i, key_j, message) in pending {
                let same_multiset = |key: &(usize, u8)| {
                    let (mut expected, mut got) = around[key].clone();
                    expected.sort();
                    got.sort();
                    expected == got
                };
                if same_multiset(&key_i) || same_multiset(&key_j) {
                    swaps += 1;
                } else {
                    *mismatches.entry(message).or_default() += 1;
                }
            }
        }
    }
    println!("intra-residue bonds with a dictionary entry: {total}, agreeing: {agreed}, resonance swaps: {swaps}");
    println!("\nmismatches:");
    for (key, count) in &mismatches {
        println!("  {count:4}  {key}");
    }
    println!("\ninter-residue bonds:");
    for (key, count) in &inter {
        println!("  {count:4}  {key}");
    }
}
