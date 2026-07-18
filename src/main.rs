use std::fs;

use prote::Cif;

fn main() {
    let cif_file = fs::read_to_string("4d1m.cif").unwrap();
    let cif = Cif::new(&cif_file);
    let block = cif.block("4D1M").unwrap();
    println!("{:#?}", block.get("_entry.id"));
    // for x in block.column("_atom_site.Cartn_x").unwrap() {
    //     // all x-coords, no intermediate Vec
    //     println!("{:#?}", x);
    // }
}
