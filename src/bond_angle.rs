use crate::{atom_pointer::AtomPointer, geometry::Point3};

#[derive(Debug)]
pub struct BondAngle<'a> {
    atom_1: AtomPointer<'a>,
    vertex: AtomPointer<'a>,
    atom_2: AtomPointer<'a>,
    value: f64,
}

impl<'a> BondAngle<'a> {
    pub fn new(
        atom_1: (AtomPointer<'a>, &Point3),
        vertex: (AtomPointer<'a>, &Point3),
        atom_2: (AtomPointer<'a>, &Point3),
    ) -> Self {
        let (atom_1, atom_1_position) = atom_1;
        let (vertex, vertex_position) = vertex;
        let (atom_2, atom_2_position) = atom_2;
        let arm_1 = *vertex_position - *atom_1_position;
        let arm_2 = *vertex_position - *atom_2_position;
        let value = arm_1.angle(&arm_2);
        Self {
            atom_1,
            vertex,
            atom_2,
            value,
        }
    }

    pub fn atoms(&self) -> (&AtomPointer<'_>, &AtomPointer<'_>, &AtomPointer<'_>) {
        (&self.atom_1, &self.vertex, &self.atom_2)
    }
}
