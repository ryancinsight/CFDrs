use crate::domain::model::NetworkBlueprint;
use crate::geometry::metadata::NodeLayoutMetadata;
use crate::geometry::Point2D;
use serde::Serialize;
use std::collections::{HashMap, VecDeque};

/// Returns `(all_positions, auto_layout_applied)` where `auto_layout_applied` contains
/// only the entries whose coordinates were determined by the topological auto-layout
/// fallback (i.e. the node had no explicit position). The second map can be serialised
/// to JSON to record and later replay the computed layout.
pub(super) fn blueprint_node_positions(
    blueprint: &NetworkBlueprint,
    box_dims: (f64, f64),
) -> (HashMap<String, Point2D>, HashMap<String, Point2D>) {
    let auto_layout = auto_layout_positions(blueprint, box_dims);
    let mut node_points = HashMap::with_capacity(blueprint.nodes.len());
    let mut auto_layout_used: HashMap<String, Point2D> = HashMap::new();
    for node_spec in &blueprint.nodes {
        let explicit = node_spec
            .layout
            .as_ref()
            .map(|layout| (layout.x_mm, layout.y_mm))
            .or_else(|| {
                node_spec
                    .metadata
                    .as_ref()
                    .and_then(|meta| meta.get::<NodeLayoutMetadata>())
                    .map(|layout| (layout.x_mm, layout.y_mm))
            })
            .or_else(|| {
                if node_spec.point == (0.0, 0.0) {
                    None
                } else {
                    Some(node_spec.point)
                }
            });
        let point = explicit.unwrap_or_else(|| {
            let fallback = auto_layout
                .get(node_spec.id.as_str())
                .copied()
                .unwrap_or((box_dims.0 * 0.5, box_dims.1 * 0.5));
            auto_layout_used.insert(node_spec.id.to_string(), fallback);
            fallback
        });
        node_points.insert(node_spec.id.as_str().to_string(), point);
    }
    (node_points, auto_layout_used)
}

#[derive(Serialize)]
struct AutoLayoutRecord<'a> {
    source: &'static str,
    box_dims: [f64; 2],
    layout: HashMap<&'a str, NodePositionRecord>,
}

#[derive(Serialize)]
struct NodePositionRecord {
    x_mm: f64,
    y_mm: f64,
}

pub(super) fn save_auto_layout_json(
    positions: &HashMap<String, Point2D>,
    box_dims: (f64, f64),
    output_path: &str,
) {
    let stem = std::path::Path::new(output_path)
        .file_stem()
        .and_then(|s| s.to_str())
        .unwrap_or("layout");
    let dir = std::path::Path::new(output_path)
        .parent()
        .unwrap_or_else(|| std::path::Path::new("."));
    let json_path = dir.join(format!("{stem}_layout.json"));

    let record = AutoLayoutRecord {
        source: "auto_layout",
        box_dims: [box_dims.0, box_dims.1],
        layout: positions
            .iter()
            .map(|(id, &(x_mm, y_mm))| (id.as_str(), NodePositionRecord { x_mm, y_mm }))
            .collect(),
    };
    if let Ok(json) = serde_json::to_string_pretty(&record) {
        let _ = std::fs::write(&json_path, json);
    }
}

fn auto_layout_positions(
    blueprint: &NetworkBlueprint,
    box_dims: (f64, f64),
) -> HashMap<String, Point2D> {
    let mut indegree: HashMap<&str, usize> = blueprint
        .nodes
        .iter()
        .map(|node| (node.id.as_str(), 0usize))
        .collect();
    let mut outgoing: HashMap<&str, Vec<&str>> = blueprint
        .nodes
        .iter()
        .map(|node| (node.id.as_str(), Vec::new()))
        .collect();

    for channel in &blueprint.channels {
        *indegree.entry(channel.to.as_str()).or_insert(0) += 1;
        outgoing
            .entry(channel.from.as_str())
            .or_default()
            .push(channel.to.as_str());
    }

    let mut queue = VecDeque::new();
    let mut depth: HashMap<&str, usize> = HashMap::new();
    for node in &blueprint.nodes {
        if indegree.get(node.id.as_str()).copied().unwrap_or(0) == 0 {
            queue.push_back(node.id.as_str());
            depth.insert(node.id.as_str(), 0);
        }
    }

    while let Some(node_id) = queue.pop_front() {
        let current_depth = depth.get(node_id).copied().unwrap_or(0);
        if let Some(next_nodes) = outgoing.get(node_id) {
            for next in next_nodes {
                let next_depth = current_depth + 1;
                depth
                    .entry(*next)
                    .and_modify(|existing| *existing = (*existing).max(next_depth))
                    .or_insert(next_depth);
                if let Some(entry) = indegree.get_mut(next) {
                    *entry = entry.saturating_sub(1);
                    if *entry == 0 {
                        queue.push_back(next);
                    }
                }
            }
        }
    }

    let max_depth = depth.values().copied().max().unwrap_or(1).max(1);
    let x_min = box_dims.0 * 0.08;
    let x_span = box_dims.0 * 0.84;
    let y_min = box_dims.1 * 0.12;
    let y_span = box_dims.1 * 0.76;

    let mut columns: HashMap<usize, Vec<&str>> = HashMap::new();
    for node in &blueprint.nodes {
        let node_depth = depth.get(node.id.as_str()).copied().unwrap_or(0);
        columns.entry(node_depth).or_default().push(node.id.as_str());
    }

    let mut positions = HashMap::with_capacity(blueprint.nodes.len());
    for (column, node_ids) in columns {
        let x = if max_depth == 0 {
            box_dims.0 * 0.5
        } else {
            x_min + x_span * (column as f64 / max_depth as f64)
        };
        let count = node_ids.len().max(1);
        for (row, node_id) in node_ids.iter().enumerate() {
            let y = if count == 1 {
                box_dims.1 * 0.5
            } else {
                y_min + y_span * (row as f64 / (count - 1) as f64)
            };
            positions.insert((*node_id).to_string(), (x, y));
        }
    }

    positions
}