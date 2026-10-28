use crate::{
    conformer::Conformer,
    geometry::Point3,
    kdtree::{KdTree, SpatialPoint},
    selection::Selection,
    structure::Atom,
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

    /// Pairs whose spheres overlap: `distance <= radius(atom_1) + radius(atom_2)`.
    /// For a fixed cutoff `d`, pass `|_| d / 2.0`.
    pub fn pairs_within<'s>(
        &'s self,
        query: &'s Selection,
        search: &'s Selection,
        radius: impl Fn(&Atom) -> f64 + 's,
    ) -> impl Iterator<Item = (usize, usize, f64)> + 's {
        let structure = self.conformer.structure();
        let positions = &self.conformer.frame().positions;
        assert_eq!(query.atoms.len(), structure.atoms.len(), "selection is from a different structure");
        assert_eq!(search.atoms.len(), structure.atoms.len(), "selection is from a different structure");
        // A query atom's partners lie within its radius plus the largest radius in `search`.
        let largest = search.indices().map(|atom| radius(&structure.atoms[atom])).fold(0.0, f64::max);
        query.indices().flat_map(move |atom_1| {
            let r_1 = radius(&structure.atoms[atom_1]);
            let mut pairs = self.candidates_within(&positions[atom_1], r_1 + largest);
            pairs.retain(|&(atom_2, distance)| {
                // Skip the pair if atom_2's own search reports it first.
                let found_from_atom_2 = atom_2 < atom_1 && query.contains(atom_2) && search.contains(atom_1);
                atom_2 != atom_1 && search.contains(atom_2) && distance <= r_1 + radius(&structure.atoms[atom_2]) && !found_from_atom_2
            });
            pairs.into_iter().map(move |(atom_2, distance)| (atom_1, atom_2, distance))
        })
    }

    fn candidates_within(&self, center: &Point3, half_size: f64) -> Vec<(usize, f64)> {
        let min = [center.x - half_size, center.y - half_size, center.z - half_size];
        let max = [center.x + half_size, center.y + half_size, center.z + half_size];
        self.tree
            .search(&self.conformer.frame().positions, min, max)
            .into_iter()
            .map(|id| (id, self.conformer.frame().positions[id].distance(center)))
            .collect()
    }

    pub fn pairs_closer_than<'s>(
        &'s self,
        query: &'s Selection,
        search: &'s Selection,
        distance: f64,
    ) -> impl Iterator<Item = (usize, usize, f64)> + 's {
        self.pairs_within(query, search, move |_| distance / 2.0)
    }
}
