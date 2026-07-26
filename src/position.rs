use crate::kdtree::SpatialPoint;

#[derive(Debug)]
pub struct Position {
    x: f64,
    y: f64,
    z: f64,
}

impl Position {
    pub fn new(x: f64, y: f64, z: f64) -> Self {
        Self { x, y, z }
    }
}

impl SpatialPoint<f64, 3> for Position {
    fn point(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}
