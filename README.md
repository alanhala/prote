# prote

[![Rust](https://github.com/alanhala/prote/actions/workflows/rust.yml/badge.svg)](https://github.com/alanhala/prote/actions/workflows/rust.yml)

prote aims to be the Rust tool for working with proteins. One project you can reach for whatever protein analysis you need, instead of a different library for each piece.

## What it does today

- Its own mmCIF lexer, parser, and typed category layer, not a wrapper around an existing library.
- Geometric bond perception from raw coordinates. It does not depend on a file declaring its own bonds.
- A generic KD-tree for spatial queries.
- Bond angles.
- Lennard-Jones interaction energy between two molecules.

The whole project is built on the Rust standard library. No runtime dependencies.

## Example

```rust
use std::fs;
use prote::{mmcif::MmCIF, Cif};

let cif_file = fs::read_to_string("some_structure.cif").unwrap();
let cif = Cif::new(&cif_file);
let mmcif = MmCIF::new(cif).unwrap();
let ensembles = mmcif.build_ensembles();

for ensemble in &ensembles {
    println!("{}: {} atoms", ensemble.topology.name, ensemble.topology.atoms.len());
}
```

## Design

The domain model is deliberately decoupled from mmCIF. A `Topology` (atoms, bonds, residues) describes chemical identity, and a `Conformer` (positions, occupancy, B-factor) describes one spatial snapshot of it. A `Molecule` pairs the two. This split exists so the same topology can be shared across many conformers: NMR ensembles, MD trajectory frames, and crystallographic alternate locations are all the same underlying idea, just with more or fewer snapshots.

## Status

Working: mmCIF parsing, bond and angle perception, spatial indexing, Lennard-Jones interaction energy.

Not yet built: bond order, partial charges, electrostatics, torsion angles, non-covalent contact detection.

## Building

```
cargo build
```

## License

Licensed under either of MIT or Apache-2.0 at your option.
