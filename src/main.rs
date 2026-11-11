use std::{
    collections::{BTreeMap, BTreeSet, HashSet},
    fs,
    path::Path,
};

use prote::{
    charge,
    contact::{Contact, ContactType},
    dot,
    mmcif::MmCIF,
    selection::Selection,
    spatial_index::SpatialIndex,
    structure::{EntityKind, PolymerKind},
    Cif,
};

fn main() {
    let cif = Cif::new(&fs::read_to_string("./cifs/7RNP.cif").unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    println!("{:?}", charge::formal_charges(&structure));
}
