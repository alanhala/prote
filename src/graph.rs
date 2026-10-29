pub fn adjacency(node_count: usize, edges: impl IntoIterator<Item = (usize, usize)>) -> Vec<Vec<(usize, usize)>> {
    let mut neighbors = vec![vec![]; node_count];
    for (edge, (a, b)) in edges.into_iter().enumerate() {
        neighbors[a].push((b, edge));
        neighbors[b].push((a, edge));
    }
    neighbors
}

pub fn within_steps(neighbors: &[Vec<(usize, usize)>], from: usize, to: usize, steps: usize) -> bool {
    from == to || steps > 0 && neighbors[from].iter().any(|&(next, _)| within_steps(neighbors, next, to, steps - 1))
}
