//! Draws the contacts between one ligand and everything around it, one node per residue.
use std::{collections::BTreeMap, fs, path::Path};

use prote::{
    atom_properties::AtomProperties,
    atom_types::AtomTypes,
    contact::{self, Contact},
    dot,
    mask::Mask,
    mmcif::MmCIF,
    spatial_index::SpatialIndex,
    structure::Structure,
    Cif,
};

// Most specific first: an edge takes the colour of the most specific kind between its two residues.
const KINDS: [(&str, &str); 5] =
    [("coordination", "#eb6834"), ("face-to-face", "#1baf7a"), ("edge-to-face", "#eda100"), ("hydrophobic", "#2a78d6"), ("van der Waals", "#8a8984")];

/// The kind's position in `KINDS`, one atom from each side, and the distance.
fn describe(structure: &Structure, contact: &Contact) -> (usize, usize, usize, f64) {
    let ring_atom = |ring: usize| structure.bonding.rings[ring].atoms[0];
    match *contact {
        Contact::MetalCoordination { atom_1, atom_2, distance } => (0, atom_1, atom_2, distance),
        Contact::FaceToFace { ring_1, ring_2, distance, .. } => (1, ring_atom(ring_1), ring_atom(ring_2), distance),
        Contact::EdgeToFace { ring_1, ring_2, distance, .. } => (2, ring_atom(ring_1), ring_atom(ring_2), distance),
        Contact::Hydrophobic { atom_1, atom_2, distance } => (3, atom_1, atom_2, distance),
        Contact::VanDerWaals { atom_1, atom_2, distance } => (4, atom_1, atom_2, distance),
    }
}

fn main() {
    // let (path, ligand_name, output) = ("./cifs/1ASN.cif", "PLP", "test.svg");
    let (path, ligand_name, output) = ("./cifs/4HHB.cif", "HEM", "test.svg");
    let mmcif = MmCIF::new(Cif::new(&fs::read_to_string(path).unwrap())).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    let index = SpatialIndex::new(&frames[0].positions);

    // Only the first copy of the ligand, so the picture shows one binding site.
    let ligand = structure.residues.iter().position(|residue| residue.name == ligand_name).expect("no residue with that name");
    let in_ligand: Vec<bool> = (0..structure.atoms.len()).map(|atom| structure.residues[ligand].atoms.contains(&atom)).collect();
    let query = Mask::from_bools(in_ligand.clone());
    let search = Mask::from_bools(in_ligand.iter().map(|&inside| !inside).collect());

    // (other residue, kind) -> (number of contacts, shortest distance)
    let mut found: BTreeMap<(usize, usize), (usize, f64)> = BTreeMap::new();
    let types = AtomTypes::new(&structure, &AtomProperties::new(&structure));
    for contact in contact::find(&structure, &types, &index, &query, &search) {
        let (kind, atom_1, atom_2, distance) = describe(&structure, &contact);
        let (residue_1, residue_2) = (structure.residue_of(atom_1), structure.residue_of(atom_2));
        let other = if residue_1 == ligand { residue_2 } else { residue_1 };
        let entry = found.entry((other, kind)).or_insert((0, f64::INFINITY));
        *entry = (entry.0 + 1, entry.1.min(distance));
    }

    let label = |residue: usize| format!("{} {}", structure.residues[residue].name, structure.residues[residue].number);
    let mut nodes = vec![dot::Node { id: ligand.to_string(), label: Some(label(ligand)), fill: Some("#e8e6df".to_string()), shape: dot::Shape::Box }];
    let mut edges = vec![];
    let mut by_residue: BTreeMap<usize, Vec<(usize, usize, f64)>> = BTreeMap::new();
    for (&(other, kind), &(count, distance)) in &found {
        by_residue.entry(other).or_default().push((kind, count, distance));
    }
    for (other, kinds) in &by_residue {
        // The counts go in the residue's node; as edge labels they pile up around the ligand.
        let mut lines = vec![label(*other)];
        lines.extend(kinds.iter().map(|&(kind, count, distance)| format!("{} x{} ({:.1} A)", KINDS[kind].0, count, distance)));
        nodes.push(dot::node(other, lines.join("\n")));
        edges.push(dot::Edge { from: ligand.to_string(), to: other.to_string(), color: Some(KINDS[kinds[0].0].1.to_string()), ..Default::default() });
    }
    // The legend: one filled box per kind that appears.
    for (kind, (name, color)) in KINDS.iter().enumerate() {
        if found.keys().any(|&(_, found_kind)| found_kind == kind) {
            let count: usize = found.iter().filter(|((_, found_kind), _)| *found_kind == kind).map(|(_, (count, _))| count).sum();
            nodes.push(dot::Node {
                id: format!("legend-{kind}"),
                label: Some(format!("{name}: {count}")),
                fill: Some(color.to_string()),
                shape: dot::Shape::Box,
            });
            println!("{name}: {count}");
        }
    }
    let title = format!("{} contacts in {} (count, shortest distance)", label(ligand), Path::new(path).file_stem().unwrap().to_str().unwrap());
    let graph = dot::Graph { title, nodes, edges };
    dot::write_svg(&graph, dot::Layout::Twopi, Path::new(output)).unwrap();
}
