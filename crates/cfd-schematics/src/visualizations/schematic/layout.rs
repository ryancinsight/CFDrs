use crate::domain::model::NetworkBlueprint;
use crate::geometry::metadata::NodeLayoutMetadata;
use crate::geometry::Point2D;
use serde::Serialize;
use std::collections::{HashMap, VecDeque};

/// Borrowed node-layout cache for schematic rendering and blueprint materialization.
///
/// # Theorem - Indexed Blueprint Layout Cache
///
/// For a fixed blueprint and box dimensions, each node receives exactly one
/// stored position in authored node order. Auto-layout positions are recorded
/// as indices into the same flat array, so consumers reuse the cache without
/// cloning node IDs or position records.
///
/// **Proof sketch**: The builder performs a single pass over the authored node
/// list, resolves each node's position once, and records the indices that fell
/// back to the topological layout. All subsequent consumers read from the same
/// borrowed cache by index or borrowed ID lookup.
pub(super) struct BlueprintNodeLayout<'a> {
    node_ids: Vec<&'a str>,
    positions: Vec<Point2D>,
    auto_layout_indices: Vec<usize>,
    index_by_id: HashMap<&'a str, usize>,
}

impl<'a> BlueprintNodeLayout<'a> {
    #[must_use]
    pub(super) fn position(&self, node_id: &str) -> Option<Point2D> {
        self.index_by_id
            .get(node_id)
            .copied()
            .and_then(|idx| self.positions.get(idx).copied())
    }

    #[must_use]
    pub(super) fn positions(&self) -> &[Point2D] {
        &self.positions
    }

    #[must_use]
    pub(super) fn index_of(&self, node_id: &str) -> Option<usize> {
        self.index_by_id.get(node_id).copied()
    }

    #[must_use]
    pub(super) fn auto_layout_indices(&self) -> &[usize] {
        &self.auto_layout_indices
    }

    #[must_use]
    pub(super) fn auto_layout_entries(&self) -> impl Iterator<Item = (&'a str, Point2D)> + '_ {
        self.auto_layout_indices
            .iter()
            .copied()
            .map(move |idx| (self.node_ids[idx], self.positions[idx]))
    }
}

/// Build the indexed node-layout cache for a blueprint.
pub(super) fn blueprint_node_positions(
    blueprint: &NetworkBlueprint,
    box_dims: (f64, f64),
) -> BlueprintNodeLayout<'_> {
    let node_count = blueprint.nodes.len();
    let mut node_ids = Vec::with_capacity(node_count);
    let mut index_by_id = HashMap::with_capacity(node_count);
    for (idx, node) in blueprint.nodes.iter().enumerate() {
        let id = node.id.as_str();
        node_ids.push(id);
        index_by_id.insert(id, idx);
    }

    let depths = auto_layout_depths(blueprint, &index_by_id);
    let auto_positions = auto_layout_positions(&depths, box_dims);

    let mut positions = Vec::with_capacity(node_count);
    let mut auto_layout_indices = Vec::new();
    for (idx, node_spec) in blueprint.nodes.iter().enumerate() {
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
            auto_layout_indices.push(idx);
            auto_positions[idx]
        });
        positions.push(point);
    }

    BlueprintNodeLayout {
        node_ids,
        positions,
        auto_layout_indices,
        index_by_id,
    }
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
    layout: &BlueprintNodeLayout<'_>,
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
        layout: layout
            .auto_layout_entries()
            .map(|(id, (x_mm, y_mm))| (id, NodePositionRecord { x_mm, y_mm }))
            .collect(),
    };
    if let Ok(json) = serde_json::to_string_pretty(&record) {
        let _ = std::fs::write(&json_path, json);
    }
}

fn auto_layout_depths(
    blueprint: &NetworkBlueprint,
    index_by_id: &HashMap<&str, usize>,
) -> Vec<usize> {
    let node_count = blueprint.nodes.len();
    let mut indegree = vec![0usize; node_count];
    let mut outgoing = vec![Vec::new(); node_count];

    for channel in &blueprint.channels {
        let Some(&from_idx) = index_by_id.get(channel.from.as_str()) else {
            continue;
        };
        let Some(&to_idx) = index_by_id.get(channel.to.as_str()) else {
            continue;
        };
        indegree[to_idx] = indegree[to_idx].saturating_add(1);
        outgoing[from_idx].push(to_idx);
    }

    let mut queue = VecDeque::new();
    let mut depth = vec![0usize; node_count];
    for idx in 0..node_count {
        if indegree[idx] == 0 {
            queue.push_back(idx);
        }
    }

    while let Some(node_idx) = queue.pop_front() {
        let current_depth = depth[node_idx];
        for next_idx in &outgoing[node_idx] {
            let next_depth = current_depth + 1;
            depth[*next_idx] = depth[*next_idx].max(next_depth);
            indegree[*next_idx] = indegree[*next_idx].saturating_sub(1);
            if indegree[*next_idx] == 0 {
                queue.push_back(*next_idx);
            }
        }
    }

    depth
}

fn auto_layout_positions(depths: &[usize], box_dims: (f64, f64)) -> Vec<Point2D> {
    let max_depth = depths.iter().copied().max().unwrap_or(0);
    let x_min = box_dims.0 * 0.08;
    let x_span = box_dims.0 * 0.84;
    let y_min = box_dims.1 * 0.12;
    let y_span = box_dims.1 * 0.76;

    let mut columns: Vec<Vec<usize>> = vec![Vec::new(); max_depth.saturating_add(1)];
    for (idx, depth) in depths.iter().copied().enumerate() {
        columns[depth].push(idx);
    }

    let mut positions = vec![(box_dims.0 * 0.5, box_dims.1 * 0.5); depths.len()];
    for (depth, node_indices) in columns.into_iter().enumerate() {
        if node_indices.is_empty() {
            continue;
        }
        let x = if max_depth == 0 {
            box_dims.0 * 0.5
        } else {
            x_min + x_span * (depth as f64 / max_depth as f64)
        };
        let count = node_indices.len();
        for (row, idx) in node_indices.into_iter().enumerate() {
            let y = if count == 1 {
                box_dims.1 * 0.5
            } else {
                y_min + y_span * (row as f64 / (count - 1) as f64)
            };
            positions[idx] = (x, y);
        }
    }

    positions
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
    use crate::geometry::metadata::NodeLayoutMetadata;

    fn layout_blueprint() -> NetworkBlueprint {
        NetworkBlueprint {
            name: "layout-cache".to_string(),
            box_dims: (100.0, 50.0),
            box_outline: Vec::new(),
            nodes: vec![
                NodeSpec::new_at("inlet", NodeKind::Inlet, (10.0, 20.0)).with_layout(
                    NodeLayoutMetadata {
                        x_mm: 10.0,
                        y_mm: 20.0,
                    },
                ),
                NodeSpec::new_at("mid", NodeKind::Junction, (0.0, 0.0)),
                NodeSpec::new_at("outlet", NodeKind::Outlet, (90.0, 30.0)),
            ],
            channels: vec![
                ChannelSpec::new_pipe("c0", "inlet", "mid", 1.0, 1.0, 1.0, 0.0),
                ChannelSpec::new_pipe("c1", "mid", "outlet", 1.0, 1.0, 1.0, 0.0),
            ],
            render_hints: None,
            topology: None,
            lineage: None,
            metadata: None,
            geometry_authored: false,
        }
    }

    #[test]
    fn blueprint_node_positions_use_flat_indexed_layout_cache() {
        let blueprint = layout_blueprint();
        let layout = blueprint_node_positions(&blueprint, blueprint.box_dims);

        assert_eq!(layout.position("inlet"), Some((10.0, 20.0)));
        assert_eq!(layout.position("mid"), Some((50.0, 25.0)));
        assert_eq!(layout.position("outlet"), Some((90.0, 30.0)));
        assert_eq!(layout.positions().get(1).copied(), Some((50.0, 25.0)));
        assert_eq!(layout.index_of("mid"), Some(1));
        assert_eq!(layout.auto_layout_indices(), &[1]);

        let auto_entries: Vec<_> = layout.auto_layout_entries().collect();
        assert_eq!(auto_entries, vec![("mid", (50.0, 25.0))]);
    }
}
