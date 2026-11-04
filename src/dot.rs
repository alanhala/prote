// Graphviz graphs as plain data: build a `Graph`, then `to_dot` or `write_svg`. Callers never write DOT.

use std::{
    io::{self, Write},
    path::Path,
    process::{Command, Stdio},
};

#[derive(Debug, Clone, Default)]
pub struct Graph {
    pub title: String,
    pub nodes: Vec<Node>,
    pub edges: Vec<Edge>,
}

#[derive(Debug, Clone, Default)]
pub struct Node {
    // Must be unique in the graph; edges refer to nodes by it.
    pub id: String,
    // None shows the id.
    pub label: Option<String>,
    // A colour name or "#rrggbb".
    pub fill: Option<String>,
    pub shape: Shape,
}

#[derive(Debug, Clone, Copy, Default)]
pub enum Shape {
    #[default]
    Ellipse,
    Box,
    Circle,
    Point,
}

#[derive(Debug, Clone, Default)]
pub struct Edge {
    pub from: String,
    pub to: String,
    pub label: Option<String>,
    pub color: Option<String>,
    // Parallel lines, e.g. 2 for a double bond. 0 is drawn as 1.
    pub lines: u8,
    pub dashed: bool,
}

// Graphviz's layout programs. Dot draws hierarchies, Neato/Fdp/Sfdp spring layouts (molecules, contacts),
// Circo circles, Twopi radial trees.
#[derive(Debug, Clone, Copy)]
pub enum Layout {
    Dot,
    Neato,
    Fdp,
    Sfdp,
    Circo,
    Twopi,
}

pub fn node(id: impl ToString, label: impl ToString) -> Node {
    Node { id: id.to_string(), label: Some(label.to_string()), ..Default::default() }
}

pub fn edge(from: impl ToString, to: impl ToString) -> Edge {
    Edge { from: from.to_string(), to: to.to_string(), ..Default::default() }
}

pub fn to_dot(graph: &Graph) -> String {
    let mut dot = format!("graph \"{}\" {{\n", escape(&graph.title));
    dot += &format!("  label=\"{}\"; labelloc=t; overlap=false;\n", escape(&graph.title));
    dot += "  node [style=filled, fillcolor=white, fontsize=10];\n  edge [fontsize=8];\n";
    for node in &graph.nodes {
        let mut attributes = vec![format!("shape={}", shape_name(node.shape))];
        attributes.extend(node.label.as_ref().map(|label| format!("label=\"{}\"", escape(label))));
        attributes.extend(node.fill.as_ref().map(|fill| format!("fillcolor=\"{}\"", escape(fill))));
        dot += &format!("  \"{}\" [{}];\n", escape(&node.id), attributes.join(", "));
    }
    for edge in &graph.edges {
        // Graphviz draws a colon-separated colour list as parallel lines.
        let color = escape(edge.color.as_deref().unwrap_or("black"));
        let mut attributes = vec![format!("color=\"{}\"", vec![color; edge.lines.max(1) as usize].join(":invis:"))];
        attributes.extend(edge.label.as_ref().map(|label| format!("label=\"{}\"", escape(label))));
        if edge.dashed {
            attributes.push("style=dashed".to_string());
        }
        dot += &format!("  \"{}\" -- \"{}\" [{}];\n", escape(&edge.from), escape(&edge.to), attributes.join(", "));
    }
    dot + "}\n"
}

// Runs the Graphviz program for `layout`, which must be installed and on PATH.
pub fn write_svg(graph: &Graph, layout: Layout, path: &Path) -> io::Result<()> {
    let mut child = Command::new(layout_program(layout))
        .arg("-Tsvg")
        .arg("-o")
        .arg(path)
        .stdin(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()?;
    child.stdin.take().expect("stdin is piped").write_all(to_dot(graph).as_bytes())?;
    let output = child.wait_with_output()?;
    if !output.status.success() {
        return Err(io::Error::other(String::from_utf8_lossy(&output.stderr).into_owned()));
    }
    Ok(())
}

fn shape_name(shape: Shape) -> &'static str {
    match shape {
        Shape::Ellipse => "ellipse",
        Shape::Box => "box",
        Shape::Circle => "circle",
        Shape::Point => "point",
    }
}

fn layout_program(layout: Layout) -> &'static str {
    match layout {
        Layout::Dot => "dot",
        Layout::Neato => "neato",
        Layout::Fdp => "fdp",
        Layout::Sfdp => "sfdp",
        Layout::Circo => "circo",
        Layout::Twopi => "twopi",
    }
}

// Atom names like O5' and entity names can hold quotes, which would end a DOT string. Newlines become
// DOT's own line break.
fn escape(text: &str) -> String {
    text.replace('\\', "\\\\").replace('"', "\\\"").replace('\n', "\\n")
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn writes_nodes_and_edges() {
        let graph = Graph { title: "t".to_string(), nodes: vec![node(1, "C4A"), node(2, "NZ")], edges: vec![edge(1, 2)] };
        let dot = to_dot(&graph);
        assert!(dot.contains("\"1\" [shape=ellipse, label=\"C4A\"];"));
        assert!(dot.contains("\"1\" -- \"2\" [color=\"black\"];"));
    }

    #[test]
    fn draws_a_double_bond_as_two_lines() {
        let graph = Graph { edges: vec![Edge { lines: 2, ..edge(1, 2) }], ..Default::default() };
        assert!(to_dot(&graph).contains("color=\"black:invis:black\""));
    }

    #[test]
    fn escapes_quotes_and_keeps_line_breaks() {
        let graph = Graph { nodes: vec![node("a", "PLP 411\nO5\"")], ..Default::default() };
        assert!(to_dot(&graph).contains("label=\"PLP 411\\nO5\\\"\""));
    }

    #[test]
    fn dashed_edges_carry_their_label() {
        let graph = Graph { edges: vec![Edge { dashed: true, label: Some("3.1 Å".to_string()), ..edge("a", "b") }], ..Default::default() };
        assert!(to_dot(&graph).contains("[color=\"black\", label=\"3.1 Å\", style=dashed]"));
    }
}
