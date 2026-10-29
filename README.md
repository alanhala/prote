# prote

[![Rust](https://github.com/alanhala/prote/actions/workflows/rust.yml/badge.svg)](https://github.com/alanhala/prote/actions/workflows/rust.yml)

A Rust library for protein structure analysis.

## What it does

- Parses mmCIF with its own lexer and parser.
- Builds a flat `Structure` (atoms, residues, entity instances, entities, bonds) and one `Frame` of positions per model and altloc.
- Bonds come from the file's `_chem_comp_bond` templates, plus backbone links between consecutive polymer residues (C–N for peptides, O3'–P for nucleic acids).
- Atom selections by instance, residue or atom.
- Spatial queries over a frame with a KD-tree.
- Contacts between two selections, excluding atoms within three bonds of each other.

## Example

```rust
use std::fs;
use prote::{contact::Contact, mmcif::MmCIF, selection::Selection, spatial_index::SpatialIndex, structure::EntityKind, Cif};

let mmcif = MmCIF::new(Cif::new(&fs::read_to_string("3PTB.cif").unwrap())).unwrap();
let (structure, frames) = mmcif.build_structure().unwrap();
let index = SpatialIndex::new(&frames[0].positions);

let protein = Selection::by_instance(&structure, |_, entity| matches!(entity.kind, EntityKind::Polymer(_)));
let ligand = Selection::by_residue(&structure, |residue| residue.name == "BEN");

for contact in Contact::contacts(&structure, &index, &ligand, &protein) {
    let residue = &structure.residues[structure.residue_of(contact.atom_2)];
    println!("{} {} {:.2} Å", residue.name, residue.number, contact.distance);
}
```

`examples/export_dot.rs` writes Graphviz files of a structure's bonds and contacts.

## Not yet

- Disulfides and other links from `_struct_conn`.
- Interaction types (H-bonds, salt bridges, metal coordination, hydrophobic, π-stacking).
- Rings.
- Structures whose file has no `_chem_comp_bond` (needs the Chemical Component Dictionary).

## License

Copyright (c) 2026 Alan Halatian. Licensed under the GNU Affero General Public License v3.0 only; see `LICENSE`.
