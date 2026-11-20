use std::{collections::HashSet, fs, path::Path, str::FromStr};

use prote::{
    atom_properties::{self, AtomProperties},
    atom_types::AtomTypes,
    bond::BondKind,
    bond_graph,
    cif::Cif,
    contact::{self, Contact},
    dot,
    mmcif::MmCIF,
    selection::Selection,
    spatial_index::SpatialIndex,
};

fn main() {
    let mmcif = MmCIF::new(Cif::new(&fs::read_to_string("./cifs/3PTB.cif").unwrap())).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    let mut nodes = vec![];
    let mut edges = vec![];
    let mut seen = HashSet::new();
    let mut water_halves = HashSet::new();
    let res = Selection::by_residue(&structure, |residue| residue.name == "BEN");
    let protein = Selection::by_instance(&structure, |instance, entity| entity.name == "BETA-TRYPSIN");
    let index = SpatialIndex::new(&frames[0].positions);
    let properties = AtomProperties::new(&structure);
    let types = AtomTypes::new(&structure, &properties);
    let mut in_contact = vec![];
    for contact in contact::find(&structure, &types, &properties, &index, &protein, &res) {
        let ring_atom = |ring: usize| structure.bonding.rings[ring].atoms[0];
        let (kind, color, atom_1, atom_2, distance) = match contact {
            Contact::MetalCoordination { atom_1, atom_2, distance } => ("coordination", "#eb6834", atom_1, atom_2, distance),
            Contact::HydrogenBond { donor, acceptor, distance } => ("hydrogen bond", "#c2255c", donor, acceptor, distance),
            Contact::FaceToFace { ring_1, ring_2, distance, .. } => ("face-to-face", "#1baf7a", ring_atom(ring_1), ring_atom(ring_2), distance),
            Contact::EdgeToFace { ring_1, ring_2, distance, .. } => ("edge-to-face", "#eda100", ring_atom(ring_1), ring_atom(ring_2), distance),
            Contact::Hydrophobic { atom_1, atom_2, distance } => ("hydrophobic", "#2a78d6", atom_1, atom_2, distance),
            Contact::VanDerWaals { atom_1, atom_2, distance } => ("van der Waals", "#8a8984", atom_1, atom_2, distance),
            Contact::Ionic { anion, cation, distance } => ("ionic", "#02f02a", anion, cation, distance),
            Contact::SaltBridge { anion, cation, distance } => ("salt bridge", "#b70cf0", anion, cation, distance),
            Contact::CationPi { cation, ring, distance, .. } => ("cation-pi", "#7048e8", cation, ring_atom(ring), distance),
            Contact::Halogen { halogen, acceptor, distance } => ("halogen bond", "#0b7285", halogen, acceptor, distance),
            // Drawn through the water: one node for it and an edge to each end.
            Contact::WaterBridge { atom_1, atom_2, water, distance_1_to_water, distance_2_to_water } => {
                if seen.insert(water) {
                    let residue = &structure.residues[structure.residue_of(water)];
                    nodes.push(dot::Node {
                        fill: Some("#c5f6fa".to_string()),
                        ..dot::node(water, format!("{} {}\n{}", residue.name, residue.number, structure.atoms[water].name))
                    });
                }
                for (atom, distance) in [(atom_1, distance_1_to_water), (atom_2, distance_2_to_water)] {
                    in_contact.extend(structure.residues[structure.residue_of(atom)].atoms.clone());
                    // Two bridges through one water share the half at their common end.
                    if water_halves.insert((atom, water)) {
                        edges.push(dot::Edge {
                            label: Some(format!("water bridge\n{distance:.1} A")),
                            color: Some("#15aabf".to_string()),
                            ..dot::edge(atom, water)
                        });
                    }
                }
                edges.push(dot::Edge { color: Some("#15aabf".to_string()), dashed: true, ..dot::edge(atom_1, atom_2) });
                continue;
            }
        };
        for atom in [atom_1, atom_2] {
            in_contact.extend(structure.residues[structure.residue_of(atom)].atoms.clone());
        }
        edges.push(dot::Edge { label: Some(format!("{kind}\n{distance:.1} A")), color: Some(color.to_string()), ..dot::edge(atom_1, atom_2) });
    }
    let drawn = res.or(&Selection::from_indices(structure.atoms.len(), in_contact));
    for bond in structure.bonding.bonds.edges().iter() {
        if drawn.contains(bond.atom_1) && drawn.contains(bond.atom_2) {
            for atom in [bond.atom_1, bond.atom_2] {
                if seen.insert(atom) {
                    let residue = &structure.residues[structure.residue_of(atom)];
                    let fill = if res.contains(atom) { "#ffd8a8" } else { "#d0ebff" };
                    nodes.push(dot::Node {
                        fill: Some(fill.to_string()),
                        ..dot::node(atom, format!("{} {}\n{}", residue.name, residue.number, structure.atoms[atom].name))
                    });
                }
            }
            match bond.kind {
                BondKind::Coordination => edges.push(dot::Edge { dashed: true, ..dot::edge(bond.atom_1, bond.atom_2) }),
                BondKind::Covalent { order, .. } => edges.push(dot::Edge { lines: order.as_u8(), ..dot::edge(bond.atom_1, bond.atom_2) }),
            }
        }
    }
    let graph = dot::Graph { title: String::from("Test"), nodes, edges };
    dot::write_svg(&graph, dot::Layout::Neato, Path::new("test.svg")).unwrap();
}
