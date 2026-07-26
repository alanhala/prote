use crate::{
    geometry::Point3,
    kdtree::{KdTree, SpatialPoint},
};

pub type SpatialIndex = KdTree<f64, 3>;

impl SpatialPoint<f64, 3> for Point3 {
    fn point(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}
