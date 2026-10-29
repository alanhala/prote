use std::fs;

use prote::{mmcif::MmCIF, selection::Selection, spatial_index::SpatialIndex, Cif};

fn main() {
    let cif = Cif::new(&fs::read_to_string("./cifs/3PTB.cif").unwrap());
    let mmcif = MmCIF::new(cif).unwrap();
    let (structure, frames) = mmcif.build_structure().unwrap();
    let ben_atoms = Selection::by_residue(&structure, |residue| residue.name == "BEN");
    let nitrogen = Selection::by_atom(&structure, |atom| atom.name == "N1" || atom.name == "N2");
    let protein_atoms = Selection::by_instance(&structure, |_, entity| entity.name == "BETA-TRYPSIN");
    let asp_189 = Selection::by_residue(&structure, |residue| residue.name == "ASP" && residue.number == 189);
    let carboxylate_oxygens = Selection::by_atom(&structure, |atom| atom.name == "OD1" || atom.name == "OD2");

    let query = &ben_atoms.and(&nitrogen);
    let search = &protein_atoms.and(&asp_189).and(&carboxylate_oxygens);
    let index = SpatialIndex::new(&frames[0].positions);
    for (atom_1, atom_2, distance) in index.pairs_closer_than(query, search, 4.0) {
        let ben_atom = &structure.atoms[atom_1];
        let protein_atom = &structure.atoms[atom_2];
        let residue = &structure.residues[structure.residue_of(atom_2)];
        println!("atom {:?} from {:?}{:?} is at {:?} from atom {:?}", protein_atom, residue.name, residue.number, distance, ben_atom);
    }
}
