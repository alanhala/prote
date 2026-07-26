use std::fs;

use prote::{kdtree::SpatialPoint, mmcif::MmCIF, Cif};

fn main() {
    let cif_file = fs::read_to_string("./cifs/4d1m.cif").unwrap();
    let cif = Cif::new(&cif_file);
    let mmcif = MmCIF::new(cif).unwrap();
    let molecules = mmcif.build_molecules();

    let zincs: Vec<_> = molecules.iter().filter(|m| m.topology.name == "ZINC ION").collect();
    let proteins: Vec<_> = molecules
        .iter()
        .filter(|m| m.topology.name == "CELLULAR TUMOR ANTIGEN P53")
        .collect();

    let cutoff = 3.0;

    for (zn_idx, zn) in zincs.iter().enumerate() {
        let zn_point = zn.conformer.positions[0].point();

        for (protein_idx, protein) in proteins.iter().enumerate() {
            let tree = protein.conformer.spatial_index();
            let min = [zn_point[0] - cutoff, zn_point[1] - cutoff, zn_point[2] - cutoff];
            let max = [zn_point[0] + cutoff, zn_point[1] + cutoff, zn_point[2] + cutoff];

            for atom_idx in tree.search(&protein.conformer.positions, min, max) {
                let atom_point = protein.conformer.positions[atom_idx].point();
                let dist = ((zn_point[0] - atom_point[0]).powi(2)
                    + (zn_point[1] - atom_point[1]).powi(2)
                    + (zn_point[2] - atom_point[2]).powi(2))
                .sqrt();

                if dist <= cutoff {
                    let atom = &protein.topology.atoms()[atom_idx];
                    let residue = protein
                        .topology
                        .residues()
                        .iter()
                        .find(|r| r.atom_range.contains(&atom_idx));
                    let residue_name = residue.map(|r| r.name.as_str()).unwrap_or("?");

                    println!(
                        "Zn #{zn_idx} <-> protein #{protein_idx}, {residue_name} {}: {dist:.2} \u{c5}",
                        atom.name,
                    );
                }
            }
        }
    }
}
