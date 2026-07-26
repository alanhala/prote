use std::marker::PhantomData;

type Id = usize;
type Point<T, const D: usize> = [T; D];

#[derive(Debug)]
struct Arena {
    nodes: Vec<Node>,
}

#[derive(Debug)]
enum Node {
    Leaf { id: Id },
    Split { split: Id, left: usize, right: usize },
}

impl Default for Arena {
    fn default() -> Self {
        Self::new()
    }
}

impl Arena {
    pub fn new() -> Self {
        Arena { nodes: vec![] }
    }

    pub fn add(&mut self, node: Node) -> usize {
        self.nodes.push(node);
        self.nodes.len() - 1
    }
}

#[derive(Debug)]
pub struct KdTree<T, const D: usize> {
    arena: Arena,
    root: usize,
    _marker: PhantomData<T>,
}

impl<const D: usize, T> KdTree<T, D>
where
    T: PartialOrd,
    T: Clone,
    T: Copy,
{
    pub fn new(points: &[Point<T, D>]) -> Self {
        assert!(!points.is_empty(), "points must not be empty");
        let mut arena = Arena::new();
        let sorted_by_axes: Vec<Vec<Id>> = (0..D)
            .map(|axis| {
                let mut ids: Vec<Id> = (0..points.len()).collect();
                ids.sort_by(|&a, &b| {
                    points[a][axis]
                        .partial_cmp(&points[b][axis])
                        .unwrap()
                        .then_with(|| points[a].iter().partial_cmp(points[b].iter()).unwrap())
                });
                ids
            })
            .collect();
        let root = Self::build(&mut arena, points, sorted_by_axes, 0).unwrap();
        Self {
            arena,
            root,
            _marker: PhantomData,
        }
    }

    pub fn search(&self, points: &[Point<T, D>], min: Point<T, D>, max: Point<T, D>) -> Vec<Id> {
        let mut reported_nodes: Vec<Id> = vec![];
        Self::report_tree(points, &self.arena, &mut reported_nodes, self.root, 0, min, max);
        reported_nodes
    }

    fn report_tree(
        points: &[Point<T, D>],
        arena: &Arena,
        reported_nodes: &mut Vec<Id>,
        node: usize,
        depth: usize,
        min: Point<T, D>,
        max: Point<T, D>,
    ) {
        let axis = depth % D;
        match arena.nodes[node] {
            Node::Leaf { id } => {
                if Self::point_contained(points[id], min, max) {
                    reported_nodes.push(id);
                }
            }
            Node::Split { split, left, right } => {
                if points[split][axis] >= min[axis] {
                    Self::report_tree(points, arena, reported_nodes, left, depth + 1, min, max);
                }
                if points[split][axis] <= max[axis] {
                    Self::report_tree(points, arena, reported_nodes, right, depth + 1, min, max);
                }
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

    fn build(arena: &mut Arena, points: &[Point<T, D>], sorted_by_axes: Vec<Vec<Id>>, depth: usize) -> Option<usize> {
        let axis = depth % D;
        match sorted_by_axes[axis].len() {
            0 => None,
            1 => Some(arena.add(Node::Leaf {
                id: sorted_by_axes[axis][0],
            })),
            n => {
                let split = (n - 1) / 2;
                let median = sorted_by_axes[axis][split];
                let (left, right): (Vec<Vec<Id>>, Vec<Vec<Id>>) = (0..D)
                    .map(|i| {
                        if i == axis {
                            let left = &sorted_by_axes[i][0..=split];
                            let right = &sorted_by_axes[i][split + 1..];
                            (left.to_vec(), right.to_vec())
                        } else {
                            let (left, right): (Vec<Id>, Vec<Id>) =
                                sorted_by_axes[i].iter().copied().partition(|&id| {
                                    points[id][axis]
                                        .partial_cmp(&points[median][axis])
                                        .unwrap()
                                        .then_with(|| points[id].partial_cmp(&points[median]).unwrap())
                                        .is_le()
                                });
                            (left, right)
                        }
                    })
                    .unzip();
                let left_idx = Self::build(arena, points, left, depth + 1).unwrap();
                let right_idx = Self::build(arena, points, right, depth + 1).unwrap();
                let node = arena.add(Node::Split {
                    split: median,
                    left: left_idx,
                    right: right_idx,
                });
                Some(node)
            }
        }
    }
}
