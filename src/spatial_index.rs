use crate::{
    geometry::Point3,
    kdtree::{KdTree, SpatialPoint},
};

impl SpatialPoint<f64, 3> for Point3 {
    fn point(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

pub struct SpatialIndex {
    tree: KdTree<f64, 3>,
}

impl SpatialIndex {
    pub fn new(positions: &[Point3]) -> Self {
        Self {
            tree: KdTree::new(positions),
        }
    }

    pub fn candidates_within<F>(
        &self,
        positions: &[Point3],
        center: &Point3,
        half_size: f64,
        mut candidate_check: F,
    ) -> Vec<(usize, f64)>
    where
        F: FnMut(usize, f64) -> bool,
    {
        let min = [center.x - half_size, center.y - half_size, center.z - half_size];
        let max = [center.x + half_size, center.y + half_size, center.z + half_size];
        self.tree
            .search(positions, min, max)
            .into_iter()
            .map(|id| (id, positions[id].distance(center)))
            .filter(|&(id, distance)| candidate_check(id, distance))
            .collect()
    }
}
