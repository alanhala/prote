use crate::{
    frame::Frame,
    geometry::{Plane3, Point3},
};

#[derive(Debug)]
pub struct Ring {
    pub atoms: Vec<usize>,
    pub aromatic: bool,
}

impl Ring {
    pub fn plane(&self, frame: &Frame) -> Plane3 {
        let positions: Vec<Point3> = self.atoms.iter().map(|&atom| frame.positions[atom]).collect();
        let (plane, _deviation) = Point3::plane(&positions);
        plane
    }
}
