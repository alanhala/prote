#[derive(Debug)]
pub struct Graph<E> {
    edges: Vec<E>,
    neighbors: Vec<Vec<(usize, usize)>>,
}

pub trait Edge {
    fn nodes(&self) -> (usize, usize);
}

impl<E: Edge> Graph<E> {
    pub fn empty() -> Self {
        Graph { edges: vec![], neighbors: vec![] }
    }

    pub fn new(node_count: usize, edges: Vec<E>) -> Self {
        let mut neighbors = vec![vec![]; node_count];
        for (edge_index, edge) in edges.iter().enumerate() {
            let (a, b) = edge.nodes();
            neighbors[a].push((b, edge_index));
            neighbors[b].push((a, edge_index));
        }
        Graph { edges, neighbors }
    }

    pub fn edges(&self) -> &[E] {
        &self.edges
    }

    pub fn neighbors(&self, node: usize) -> &[(usize, usize)] {
        &self.neighbors[node]
    }

    pub fn neighbors_within_steps(&self, from: usize, steps: usize) -> Vec<Vec<(usize, usize)>> {
        let mut visited = vec![false; self.neighbors.len()];
        visited[from] = true;
        let mut outermost = vec![from];
        let mut by_distance = Vec::with_capacity(steps);
        for _ in 0..steps {
            let mut at_distance = vec![];
            for &node in &outermost {
                for &(next, edge) in self.neighbors(node) {
                    if !visited[next] {
                        visited[next] = true;
                        at_distance.push((next, edge));
                    }
                }
            }
            outermost = at_distance.iter().map(|&(node, _)| node).collect();
            by_distance.push(at_distance);
        }
        by_distance
    }

    pub fn within_steps(&self, from: usize, to: usize, steps: usize) -> bool {
        from == to || steps > 0 && self.neighbors[from].iter().any(|&(next, _)| self.within_steps(next, to, steps - 1))
    }
}
