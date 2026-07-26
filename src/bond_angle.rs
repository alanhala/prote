use crate::geometry::Point3;

pub struct BondAngle {
    atom_1: usize,
    vertex: usize,
    atom_2: usize,
}

impl BondAngle {
    pub fn new(atom_1: usize, vertex: usize, atom_2: usize) -> Self {
        Self { atom_1, vertex, atom_2 }
    }

    pub fn atoms(&self) -> (usize, usize, usize) {
        (self.atom_1, self.vertex, self.atom_2)
    }

    pub fn value(&self, positions: &[Point3]) -> f64 {
        let vertex = positions[self.vertex];
        let arm_1 = positions[self.atom_1] - vertex;
        let arm_2 = positions[self.atom_2] - vertex;
        arm_1.angle(&arm_2)
    }
}
