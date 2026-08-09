use crate::geometry::Point3;

#[derive(Debug)]
pub struct Frame {
    pub id: usize,
    pub positions: Vec<Point3>,
    pub occupancies: Vec<f32>,
    pub b_factors: Vec<f32>,
}

impl PartialEq for Frame {
    fn eq(&self, other: &Self) -> bool {
        self.id == other.id
    }
}
impl Frame {
    pub fn new(id: usize, positions: Vec<Point3>, occupancies: Vec<f32>, b_factors: Vec<f32>) -> Self {
        Self {
            id,
            positions,
            occupancies,
            b_factors,
        }
    }
}
