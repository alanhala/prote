use crate::geometry::Plane3;

#[derive(Debug)]
pub struct Ring {
    pub atoms: Vec<usize>,
    pub plane: Plane3,
    pub coplanar: bool,
}

impl Ring {
    // pub fn new(atom_pointers: Vec<AtomPointer>, conformers: &[&Conformer]) -> Option<Self> {
    //     if !(5..20).contains(&atom_pointers.len()) {
    //         return None;
    //     }
    //     let positions: Vec<Point3> = atom_pointers
    //         .iter()
    //         .map(|atom_pointer| conformers[atom_pointer.conformer_index].frame.positions[atom_pointer.atom_index])
    //         .collect();
    //     let (plane, deviation) = Point3::plane(&positions);
    //     Some(Self {
    //         atom_pointers,
    //         plane,
    //         coplanar: deviation <= 0.2,
    //     })
    // }
}
