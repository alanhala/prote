// Writes Graphviz files for one structure: one graph per entity instance (atoms and their bonds), and
// graphs of the contacts from polymers to non-polymers and to waters.
//
//   cargo run --example export_dot -- cifs/3PTB.cif target/dot
//   neato -Tsvg target/dot/3PTB_A.dot -o 3PTB_A.svg

use std::{collections::BTreeSet, fs, path::Path};

use prote::{
    bond::Bond,
    bond_order::BondOrder,
    conformer::Conformer,
    contact::Contact,
    element::Element,
    graph,
    mmcif::MmCIF,
    selection::Selection,
    spatial_index::SpatialIndex,
    structure::{EntityKind, Structure},
    Cif,
};

fn main() {
    let mut args = std::env::args().skip(1);
    let cif_path = args.next().unwrap_or("cifs/3PTB.cif".to_string());
    let out_dir = args.next().unwrap_or("target/dot".to_string());
    let name = Path::new(&cif_path).file_stem().unwrap().to_string_lossy().to_string();
    fs::create_dir_all(&out_dir).unwrap();

    let mmcif = MmCIF::new(Cif::new(&fs::read_to_string(&cif_path).unwrap())).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    let frame = &frames[0];

    for instance in &structure.instances {
        let residues = &structure.residues[instance.residues.clone()];
        let atoms = residues[0].atoms.start..residues[residues.len() - 1].atoms.end;
        let entity = &structure.entities[instance.entity];

        let mut dot = format!("graph \"{name} {} {}\" {{\n  node [style=filled, fontsize=10];\n", instance.id, escape(&entity.name));
        for atom in atoms.clone() {
            dot += &format!("  {atom} [label=\"{}\", fillcolor=\"{}\"];\n", atom_label(&structure, atom), color(structure.atoms[atom].element));
        }
        for bond in structure.bonds.iter().filter(|bond| atoms.contains(&bond.atom_1) && atoms.contains(&bond.atom_2)) {
            let between_residues = residue_of(&structure, bond.atom_1) != residue_of(&structure, bond.atom_2);
            dot += &format!(
                "  {} -- {} [label=\"{}\", color=\"{}\"{}];\n",
                bond.atom_1,
                bond.atom_2,
                bond_label(bond, bond.length(frame)),
                bond_color(bond.order),
                if between_residues { ", penwidth=3" } else { "" },
            );
        }
        dot += "}\n";
        let path = format!("{out_dir}/{name}_{}.dot", instance.id);
        fs::write(&path, dot).unwrap();
        println!("{path}: instance {} ({}, {:?})", instance.id, entity.name, entity.kind);
    }

    let index = SpatialIndex::new(&frame.positions);
    let polymers = Selection::by_instance(&structure, |_, entity| matches!(entity.kind, EntityKind::Polymer(_)));
    let non_polymers = Selection::by_instance(&structure, |_, entity| entity.kind == EntityKind::NonPolymer);
    let waters = Selection::by_instance(&structure, |_, entity| entity.kind == EntityKind::Water);
    for (label, search) in [("non-polymers", &non_polymers), ("waters", &waters)] {
        let contacts = Contact::contacts(&structure, &index, &polymers, search);
        let path = format!("{out_dir}/{name}_contacts_{label}.dot");
        fs::write(&path, contacts_dot(&structure, &format!("{name} polymers - {label}"), &contacts)).unwrap();
        println!("{path}: {} contacts between polymers and {label}", contacts.len());
    }
}

fn contacts_dot(structure: &Structure, title: &str, contacts: &[Contact]) -> String {
    let mut dot = format!("graph \"{title}\" {{\n  node [style=filled, fontsize=10];\n");
    let atoms: BTreeSet<usize> = contacts.iter().flat_map(|contact| [contact.atom_1, contact.atom_2]).collect();
    for &atom in &atoms {
        let instance = &structure.instances[instance_of(structure, residue_of(structure, atom))];
        dot += &format!(
            "  {atom} [label=\"{}\\n{}\", fillcolor=\"{}\"];\n",
            instance.id,
            atom_label(structure, atom),
            color(structure.atoms[atom].element)
        );
    }
    for contact in contacts {
        dot += &format!("  {} -- {} [label=\"{:.2} Å\", style=dashed];\n", contact.atom_1, contact.atom_2, contact.distance);
    }
    dot + "}\n"
}

fn residue_of(structure: &Structure, atom: usize) -> usize {
    structure.residues.partition_point(|residue| residue.atoms.end <= atom)
}

fn instance_of(structure: &Structure, residue: usize) -> usize {
    structure.instances.partition_point(|instance| instance.residues.end <= residue)
}

fn atom_label(structure: &Structure, atom: usize) -> String {
    let residue = &structure.residues[residue_of(structure, atom)];
    format!("{} {}\\n{}", residue.name, residue.number, escape(&structure.atoms[atom].name))
}

fn bond_label(bond: &Bond, length: f64) -> String {
    let order = match bond.order {
        BondOrder::Single => "sing",
        BondOrder::Double => "doub",
        BondOrder::Triple => "trip",
        BondOrder::Quadruple => "quad",
        BondOrder::Pi => "pi",
        BondOrder::Aromatic => "arom",
    };
    format!("{order} {length:.2} Å{}", if bond.aromatic { " ar" } else { "" })
}

// Graphviz draws a colon-separated colour list as parallel lines, so a double bond shows as two.
fn bond_color(order: BondOrder) -> &'static str {
    match order {
        BondOrder::Double => "black:invis:black",
        BondOrder::Triple => "black:invis:black:invis:black",
        _ => "black",
    }
}

// CPK colours.
fn color(element: Element) -> &'static str {
    match element.symbol() {
        "C" => "#c8c8c8",
        "N" => "#8f8fff",
        "O" => "#f00000",
        "S" => "#ffc832",
        "P" => "#ffa500",
        "H" => "#ffffff",
        _ => "#ff1493",
    }
}

// Atom names like O5' and entity names can hold quotes, which would end a DOT string.
fn escape(text: &str) -> String {
    text.replace('\\', "\\\\").replace('"', "\\\"")
}
