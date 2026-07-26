use crate::position::Position;

#[derive(Debug)]
pub struct Conformer {
    pub positions: Vec<Position>,
    pub occupancies: Vec<f32>,
    pub b_factors: Vec<f32>,
}

impl Conformer {
    pub fn new(positions: Vec<Position>, occupancies: Vec<f32>, b_factors: Vec<f32>) -> Self {
        Self {
            positions,
            occupancies,
            b_factors,
        }
    }
}
