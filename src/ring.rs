use crate::{
    atom_pointer::AtomPointer,
    conformer::{self, Conformer},
    geometry::{Plane3, Point3, Vector3},
};

#[derive(Debug)]
pub struct Ring {
    pub atom_pointers: Vec<AtomPointer>,
    pub plane: Plane3,
    pub coplanar: bool,
}

impl Ring {
    pub fn new(atom_pointers: Vec<AtomPointer>, conformers: &[&Conformer]) -> Option<Self> {
        if !(5..=8).contains(&atom_pointers.len()) {
            return None;
        }
        let positions: Vec<Point3> = atom_pointers
            .iter()
            .map(|atom_pointer| conformers[atom_pointer.conformer_index].frame.positions[atom_pointer.atom_index])
            .collect();
        let (plane, deviation) = Point3::plane(&positions);
        Some(Self {
            atom_pointers,
            plane,
            coplanar: deviation <= 0.2,
        })
    }
}
