use crate::geometry::{Plane3, Point3};

#[derive(Debug)]
pub struct Ring {
    pub atoms: Vec<usize>,
    pub aromatic: bool,
}

impl Ring {
    pub fn plane(&self, positions: &[Point3]) -> Plane3 {
        let positions: Vec<Point3> = self.atoms.iter().map(|&atom| positions[atom]).collect();
        let (plane, _deviation) = Point3::plane(&positions);
        plane
    }
}
