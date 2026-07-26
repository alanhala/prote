use std::collections::HashMap;

type Id = usize;
type Point<T, const D: usize> = [T; D];

#[derive(Debug)]
struct Arena<T> {
    nodes: Vec<Node<T>>,
}

#[derive(Debug)]
struct Node<T> {
    id: Id,
    value: T,
    left: Option<usize>,
    right: Option<usize>,
}

impl<T> Node<T> {
    fn is_leaf(&self) -> bool {
        self.left.is_none() && self.right.is_none()
    }
}

impl<T> Default for Arena<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T> Arena<T> {
    pub fn new() -> Self {
        Arena { nodes: vec![] }
    }

    pub fn add(&mut self, id: Id, value: T) -> usize {
        self.nodes.push(Node {
            id,
            value,
            left: None,
            right: None,
        });
        self.nodes.len() - 1
    }
}

#[derive(Debug)]
pub struct KdTree<T, const D: usize>
where
    T: PartialOrd,
    T: Clone,
    T: Copy,
{
    arena: Arena<Point<T, D>>,
    root: usize,
}

impl<const D: usize, T> KdTree<T, D>
where
    T: PartialOrd,
    T: Clone,
    T: Copy,
{
    pub fn new(points: HashMap<Id, Point<T, D>>) -> Self {
        assert!(!points.is_empty(), "points must not be empty");
        let mut arena = Arena::new();
        let sorted_by_axes: Vec<Vec<Id>> = (0..D)
            .map(|i| {
                let mut sorted_points: Vec<(&Id, &Point<T, D>)> = points.iter().collect();
                sorted_points.sort_by(|a, b| {
                    a.1[i]
                        .partial_cmp(&b.1[i])
                        .unwrap()
                        .then_with(|| a.1.iter().partial_cmp(b.1).unwrap())
                });
                sorted_points.iter().map(|p| *p.0).collect()
            })
            .collect();
        let root = Self::build(&mut arena, &points, sorted_by_axes, 0).unwrap();
        Self { arena, root }
    }

    pub fn search(&self, min: Point<T, D>, max: Point<T, D>) -> Vec<Id> {
        let mut reported_nodes: Vec<Id> = vec![];
        Self::report_tree(&self.arena, &mut reported_nodes, self.root, 0, min, max);
        reported_nodes
    }

    fn report_tree(
        arena: &Arena<Point<T, D>>,
        reported_nodes: &mut Vec<Id>,
        node: usize,
        depth: usize,
        min: Point<T, D>,
        max: Point<T, D>,
    ) {
        let axis = depth % D;
        if arena.nodes[node].is_leaf() {
            if Self::point_contained(arena.nodes[node].value, min, max) {
                reported_nodes.push(arena.nodes[node].id);
            }
        } else {
            if arena.nodes[node].value[axis] >= min[axis] {
                Self::report_tree(
                    arena,
                    reported_nodes,
                    arena.nodes[node].left.unwrap(),
                    depth + 1,
                    min,
                    max,
                );
            }
            if arena.nodes[node].value[axis] <= max[axis] {
                Self::report_tree(
                    arena,
                    reported_nodes,
                    arena.nodes[node].right.unwrap(),
                    depth + 1,
                    min,
                    max,
                );
            }
        }
    }

    fn point_contained(point: Point<T, D>, min: Point<T, D>, max: Point<T, D>) -> bool {
        point
            .iter()
            .zip(min.iter())
            .zip(max.iter())
            .all(|((vi, lo), hi)| vi >= lo && vi <= hi)
    }

    fn build(
        arena: &mut Arena<Point<T, D>>,
        points: &HashMap<Id, Point<T, D>>,
        sorted_by_axes: Vec<Vec<Id>>,
        depth: usize,
    ) -> Option<usize> {
        let axis = depth % D;
        match sorted_by_axes[axis].len() {
            0 => None,
            1 => Some(arena.add(sorted_by_axes[axis][0], points[&sorted_by_axes[axis][0]])),
            n => {
                let split = (n - 1) / 2;
                let median = &sorted_by_axes[axis][split];
                let node = arena.add(sorted_by_axes[axis][split], points[median]);
                let (left, right): (Vec<Vec<Id>>, Vec<Vec<Id>>) = (0..D)
                    .map(|i| {
                        if i == axis {
                            let left = &sorted_by_axes[i][0..=split];
                            let right = &sorted_by_axes[i][split + 1..];
                            (left.to_vec(), right.to_vec())
                        } else {
                            let (left, right): (Vec<Id>, Vec<Id>) = sorted_by_axes[i].iter().partition(|p| {
                                points[p][axis]
                                    .partial_cmp(&points[median][axis])
                                    .unwrap()
                                    .then_with(|| points[p].partial_cmp(&points[median]).unwrap())
                                    .is_le()
                            });
                            (left, right)
                        }
                    })
                    .unzip();
                arena.nodes[node].left = Self::build(arena, &points, left, depth + 1);
                arena.nodes[node].right = Self::build(arena, &points, right, depth + 1);
                Some(node)
            }
        }
    }
}
