use crate::geometry::Point3;

#[derive(Debug)]
pub struct Conformer {
    pub positions: Vec<Point3>,
    pub occupancies: Vec<f32>,
    pub b_factors: Vec<f32>,
}

impl Conformer {
    pub fn new(positions: Vec<Point3>, occupancies: Vec<f32>, b_factors: Vec<f32>) -> Self {
        Self {
            positions,
            occupancies,
            b_factors,
        }
    }
}
