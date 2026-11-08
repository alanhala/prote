use std::collections::VecDeque;

use crate::mask::{Mask, XorBasis};

#[derive(Debug)]
pub struct Graph<E> {
    edges: Vec<E>,
    neighbors: Vec<Vec<(usize, usize)>>,
}

pub trait Edge {
    fn nodes(&self) -> (usize, usize);
}

impl Edge for (usize, usize) {
    fn nodes(&self) -> (usize, usize) {
        *self
    }
}

struct CandidateRing<E> {
    pub nodes: Vec<usize>,
    pub edges: Mask<E>,
}

pub struct Ring {
    pub nodes: Vec<usize>,
    pub edges: Vec<usize>,
}

struct ShortestPathTree {
    came_from: Vec<Option<(usize, usize)>>, // (parent, edge); None for the root
}

#[derive(Debug)]
pub struct Component {
    pub nodes: Vec<usize>,
    pub edges: Vec<usize>,
}

impl ShortestPathTree {
    fn new<E: Edge>(graph: &Graph<E>, root: usize) -> Self {
        let mut came_from = vec![None; graph.neighbors.len()];
        let mut visited = vec![false; graph.neighbors.len()];
        visited[root] = true;
        let mut queue = VecDeque::from([root]);
        while let Some(node) = queue.pop_front() {
            for &(next, edge) in graph.neighbors(node) {
                if !visited[next] {
                    visited[next] = true;
                    came_from[next] = Some((node, edge));
                    queue.push_back(next);
                }
            }
        }
        ShortestPathTree { came_from }
    }

    fn is_tree_edge(&self, edge: usize, a: usize, b: usize) -> bool {
        self.came_from[a] == Some((b, edge)) || self.came_from[b] == Some((a, edge))
    }

    fn path_to_root(&self, node: usize) -> (Vec<usize>, Vec<usize>) {
        let mut nodes = vec![node];
        let mut edges = vec![];
        let mut current = node;
        while let Some((parent, edge)) = self.came_from[current] {
            nodes.push(parent);
            edges.push(edge);
            current = parent;
        }
        nodes.reverse();
        (nodes, edges)
    }
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

    pub fn rings(&self) -> Vec<Ring> {
        self.rings_where(|_| true)
    }

    pub fn rings_where(&self, keep_edge: impl Fn(usize) -> bool) -> Vec<Ring> {
        let bridges = self.bridges_where(&keep_edge);
        let components: Vec<Component> = self
            .components_where(|edge| keep_edge(edge) && !bridges.contains(edge))
            .into_iter()
            .filter(|component| component.nodes.len() > 1)
            .collect();

        let mut local_index = vec![0; self.neighbors.len()];
        for component in &components {
            for (i, &node) in component.nodes.iter().enumerate() {
                local_index[node] = i;
            }
        }
        let mut rings = vec![];
        for component in &components {
            let edges: Vec<(usize, usize)> = component
                .edges
                .iter()
                .map(|&edge| {
                    let (a, b) = self.edges[edge].nodes();
                    (local_index[a], local_index[b])
                })
                .collect();
            let subgraph = Graph::new(component.nodes.len(), edges);
            for ring in subgraph.smallest_rings() {
                rings.push(Ring {
                    nodes: ring.nodes.iter().map(|&i| component.nodes[i]).collect(),
                    edges: ring.edges.iter().map(|&j| component.edges[j]).collect(),
                });
            }
        }
        rings
    }

    pub fn bridges(&self) -> Mask<E> {
        self.bridges_where(|_| true)
    }

    /// Bridges of the graph without the edges that fail `keep_edge`. Those edges are never bridges.
    pub fn bridges_where(&self, keep_edge: impl Fn(usize) -> bool) -> Mask<E> {
        let mut discovery_time = vec![None; self.neighbors.len()];
        let mut lowest_reachable_discovery_time = vec![0; self.neighbors.len()];
        let mut is_bridge = vec![false; self.edges.len()];
        let mut timer = 0;
        for node_id in 0..self.neighbors.len() {
            if discovery_time[node_id].is_none() {
                self.bridges_step(node_id, None, &keep_edge, &mut discovery_time, &mut lowest_reachable_discovery_time, &mut timer, &mut is_bridge)
            }
        }
        Mask::from_bools(is_bridge)
    }

    fn bridges_step(
        &self,
        node_id: usize,
        parent_edge: Option<usize>,
        keep_edge: &impl Fn(usize) -> bool,
        discovery_time: &mut [Option<usize>],
        lowest_reachable_discovery_time: &mut [usize],
        timer: &mut usize,
        is_bridge: &mut [bool],
    ) {
        *timer += 1;
        let node_time = *timer;
        discovery_time[node_id] = Some(node_time);
        lowest_reachable_discovery_time[node_id] = node_time;
        for &(neighbor_id, edge) in self.neighbors(node_id) {
            if Some(edge) == parent_edge || !keep_edge(edge) {
                continue;
            }
            match discovery_time[neighbor_id] {
                None => {
                    self.bridges_step(neighbor_id, Some(edge), keep_edge, discovery_time, lowest_reachable_discovery_time, timer, is_bridge);
                    lowest_reachable_discovery_time[node_id] =
                        lowest_reachable_discovery_time[node_id].min(lowest_reachable_discovery_time[neighbor_id]);
                    if lowest_reachable_discovery_time[neighbor_id] > node_time {
                        is_bridge[edge] = true;
                    }
                }
                Some(neighbor_time) => {
                    lowest_reachable_discovery_time[node_id] = lowest_reachable_discovery_time[node_id].min(neighbor_time);
                }
            }
        }
    }

    pub fn components(&self) -> Vec<Component> {
        self.components_where(|_| true)
    }

    pub fn components_where(&self, keep_edge: impl Fn(usize) -> bool) -> Vec<Component> {
        let mut used = vec![false; self.neighbors.len()];
        let mut found_components = vec![];
        for node_id in 0..self.neighbors.len() {
            if !used[node_id] {
                let mut component = Component { nodes: vec![], edges: vec![] };
                self.component_step(node_id, &keep_edge, &mut used, &mut component);
                found_components.push(component);
            }
        }
        found_components
    }

    fn component_step(&self, node: usize, keep_edge: &impl Fn(usize) -> bool, used: &mut [bool], component: &mut Component) {
        used[node] = true;
        component.nodes.push(node);
        for &(neighbor, edge) in self.neighbors(node) {
            if !keep_edge(edge) {
                continue;
            }
            if node < neighbor {
                component.edges.push(edge);
            }
            if !used[neighbor] {
                self.component_step(neighbor, keep_edge, used, component);
            }
        }
    }

    pub fn within_steps(&self, from: usize, to: usize, steps: usize) -> bool {
        from == to || steps > 0 && self.neighbors[from].iter().any(|&(next, _)| self.within_steps(next, to, steps - 1))
    }

    fn smallest_rings(&self) -> Vec<Ring> {
        let k = self.edges.len() - self.neighbors.len() + 1;
        let mut candidates: Vec<CandidateRing<E>> = self.ring_roots().into_iter().flat_map(|root| self.candidate_rings(root)).collect();
        // Same-size rings are ordered by their edges too, so duplicates end up next to each other for dedup_by.
        candidates.sort_by_cached_key(|ring| (ring.nodes.len(), ring.edges.indices().collect::<Vec<_>>()));
        candidates.dedup_by(|a, b| a.edges == b.edges);

        let mut rings = vec![];
        let mut basis = XorBasis::new();
        for candidate in candidates {
            if basis.len() == k {
                break;
            }
            if basis.insert(&candidate.edges) {
                rings.push(Ring { nodes: candidate.nodes, edges: candidate.edges.indices().collect() });
            }
        }
        rings
    }

    fn ring_roots(&self) -> Vec<usize> {
        let branch_nodes: Vec<usize> = (0..self.neighbors.len()).filter(|&node| self.neighbors[node].len() >= 3).collect();
        // A lone ring has no branch nodes, and then every node is on it.
        if branch_nodes.is_empty() && !self.neighbors.is_empty() {
            vec![0]
        } else {
            branch_nodes
        }
    }

    fn candidate_rings(&self, root: usize) -> Vec<CandidateRing<E>> {
        let shortest_path_tree = ShortestPathTree::new(self, root);
        let mut candidates = vec![];
        for (edge_id, edge) in self.edges.iter().enumerate() {
            let (x, y) = edge.nodes();
            if shortest_path_tree.is_tree_edge(edge_id, x, y) {
                continue;
            }
            let (nodes_x, edges_x) = shortest_path_tree.path_to_root(x);
            let (nodes_y, edges_y) = shortest_path_tree.path_to_root(y);
            let shares_only_root = !nodes_x[1..].iter().any(|node| nodes_y[1..].contains(node));
            if !shares_only_root {
                continue;
            }
            let mut nodes = nodes_x;
            nodes.extend(nodes_y[1..].iter().rev());
            let mut edges = edges_x;
            edges.extend(edges_y);
            edges.push(edge_id);
            candidates.push(CandidateRing { nodes, edges: Mask::from_indices(self.edges.len(), edges) });
        }
        candidates
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    // Each ring's nodes sorted, then the rings sorted: compares which rings were found, not their walk order.
    fn sorted(rings: Vec<Ring>) -> Vec<Vec<usize>> {
        let mut rings: Vec<Vec<usize>> = rings
            .into_iter()
            .map(|ring| {
                let mut nodes = ring.nodes;
                nodes.sort();
                nodes
            })
            .collect();
        rings.sort();
        rings
    }

    fn sorted_rings(graph: &Graph<(usize, usize)>) -> Vec<Vec<usize>> {
        sorted(graph.smallest_rings())
    }

    #[test]
    fn fused_triangles_give_two_triangles_not_the_square() {
        // Graph 5 in docs/bridges_exercise.md: A=0 B=1 C=2 D=3.
        let graph = Graph::new(4, vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3)]);
        assert_eq!(sorted_rings(&graph), vec![vec![0, 1, 2], vec![1, 2, 3]]);
    }

    #[test]
    fn square_is_skipped_as_a_combination_before_reaching_k() {
        // H2 in docs/bridges_exercise.md: A=0 B=1 C=2 D=3 E=4 F=5 G=6.
        let graph = Graph::new(7, vec![(0, 1), (0, 2), (1, 2), (1, 3), (2, 3), (2, 6), (3, 4), (4, 5), (5, 6)]);
        assert_eq!(sorted_rings(&graph), vec![vec![0, 1, 2], vec![1, 2, 3], vec![2, 3, 4, 5, 6]]);
    }

    #[test]
    fn rings_across_a_bridge_keep_their_original_nodes() {
        // Triangles 0-1-2 and 3-4-5, joined by the bridge 2–3.
        let graph = Graph::new(6, vec![(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5), (5, 3)]);
        assert_eq!(sorted(graph.rings()), vec![vec![0, 1, 2], vec![3, 4, 5]]);
    }

    #[test]
    fn rings_where_drops_loops_closed_by_a_skipped_edge() {
        // Triangle 0-1-2 and square 2-3-4-5, which only edge 6 (5–2) closes.
        let graph = Graph::new(6, vec![(0, 1), (1, 2), (2, 0), (2, 3), (3, 4), (4, 5), (5, 2)]);
        assert_eq!(sorted(graph.rings()), vec![vec![0, 1, 2], vec![2, 3, 4, 5]]);
        assert_eq!(sorted(graph.rings_where(|edge| edge != 6)), vec![vec![0, 1, 2]]);
    }

    #[test]
    fn lone_ring_without_branch_nodes_is_found() {
        let graph = Graph::new(6, vec![(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0)]);
        assert_eq!(sorted_rings(&graph), vec![vec![0, 1, 2, 3, 4, 5]]);
    }
}
