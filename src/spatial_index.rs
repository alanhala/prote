use crate::{
    conformer::Conformer,
    geometry::Point3,
    kdtree::{KdTree, SpatialPoint},
};

pub struct IndexedStructure<'a> {
    conformer: Conformer<'a>,
    tree: KdTree<f64, 3>,
}

impl SpatialPoint<f64, 3> for Point3 {
    fn point(&self) -> [f64; 3] {
        [self.x, self.y, self.z]
    }
}

impl<'a> IndexedStructure<'a> {
    pub fn new(conformer: Conformer<'a>) -> Self {
        let tree = KdTree::new(&conformer.frame().positions);
        Self { conformer, tree }
    }

    pub fn conformer(&self) -> Conformer<'a> {
        self.conformer
    }

    pub fn candidates_within(&self, center: &Point3, half_size: f64) -> Vec<(usize, f64)> {
        let min = [center.x - half_size, center.y - half_size, center.z - half_size];
        let max = [center.x + half_size, center.y + half_size, center.z + half_size];
        self.tree
            .search(&self.conformer.frame().positions, min, max)
            .into_iter()
            .map(|id| (id, self.conformer.frame().positions[id].distance(center)))
            .collect()
    }
}
