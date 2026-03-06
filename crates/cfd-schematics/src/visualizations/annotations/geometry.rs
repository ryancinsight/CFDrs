use crate::domain::model::{ChannelSpec, NetworkBlueprint};
use crate::domain::therapy_metadata::{TherapyZone, TherapyZoneMetadata};
use crate::geometry::metadata::{ChannelVenturiSpec, VenturiGeometryMetadata};
use crate::geometry::{ChannelSystem, Point2D};
use petgraph::algo::astar;
use petgraph::{Directed, Graph};
use std::collections::HashMap;

use super::MarkerRole;

/// Infer inlet/outlet terminals from graph-degree-1 nodes at extreme X.
#[must_use]
pub fn infer_terminal_nodes_by_x(system: &ChannelSystem) -> Option<(usize, usize)> {
    let degrees = node_degrees(system);
    let mut terminals: Vec<usize> = degrees
        .iter()
        .enumerate()
        .filter_map(|(idx, (in_deg, out_deg))| {
            if in_deg + out_deg == 1 {
                Some(idx)
            } else {
                None
            }
        })
        .collect();

    if terminals.len() < 2 {
        return None;
    }

    terminals.sort_by(|a, b| {
        let ax = system.nodes[*a].point.0;
        let bx = system.nodes[*b].point.0;
        ax.total_cmp(&bx)
    });

    Some((terminals[0], terminals[terminals.len() - 1]))
}

/// Compute `(in_degree, out_degree)` for each node index.
#[must_use]
pub fn node_degrees(system: &ChannelSystem) -> Vec<(usize, usize)> {
    let mut degrees = vec![(0_usize, 0_usize); system.nodes.len()];
    for channel in &system.channels {
        if channel.from_node < degrees.len() {
            degrees[channel.from_node].1 += 1;
        }
        if channel.to_node < degrees.len() {
            degrees[channel.to_node].0 += 1;
        }
    }
    degrees
}

/// Classify each node into a visualization marker role.
#[must_use]
pub fn classify_node_roles(system: &ChannelSystem) -> HashMap<usize, MarkerRole> {
    let degrees = node_degrees(system);
    let terminals = infer_terminal_nodes_by_x(system);

    let mut roles = HashMap::with_capacity(system.nodes.len());
    for (idx, &(in_deg, out_deg)) in degrees.iter().enumerate() {
        let role = match terminals {
            Some((inlet_idx, _)) if idx == inlet_idx => MarkerRole::Inlet,
            Some((_, outlet_idx)) if idx == outlet_idx => MarkerRole::Outlet,
            _ if out_deg > 1 => MarkerRole::Split,
            _ if in_deg > 1 => MarkerRole::Merge,
            _ => MarkerRole::Internal,
        };
        roles.insert(idx, role);
    }

    roles
}

/// Extract venturi throat count from blueprint metadata (`ChannelVenturiSpec`).
#[must_use]
pub fn throat_count_from_blueprint_metadata(blueprint: &NetworkBlueprint) -> usize {
    let total = blueprint
        .channels
        .iter()
        .filter_map(|channel| channel_venturi_spec(channel))
        .filter(|spec| spec.is_ctc_stream)
        .map(|spec| usize::from(spec.n_throats))
        .sum::<usize>();
    if total > 0 {
        total
    } else {
        blueprint
            .channels
            .iter()
            .filter(|channel| {
                channel.venturi_geometry.is_some()
                    || channel
                        .metadata
                        .as_ref()
                        .is_some_and(|meta| meta.contains::<VenturiGeometryMetadata>())
            })
            .count()
    }
}

/// Detect whether the blueprint carries therapy-zone channel metadata.
#[must_use]
pub fn therapy_zone_presence(blueprint: &NetworkBlueprint) -> (bool, bool) {
    let mut has_target = false;
    let mut has_bypass = false;

    for channel in &blueprint.channels {
        let Some(metadata) = channel.metadata.as_ref() else {
            continue;
        };
        let Some(zone_meta) = metadata.get::<TherapyZoneMetadata>() else {
            continue;
        };

        match zone_meta.zone {
            TherapyZone::CancerTarget => has_target = true,
            TherapyZone::HealthyBypass => has_bypass = true,
            TherapyZone::MixedFlow => {}
        }
    }

    (has_target, has_bypass)
}

/// Compute a center-biased inlet→outlet path on the node graph.
#[must_use]
pub fn center_biased_main_path(system: &ChannelSystem) -> Vec<Point2D> {
    let Some((inlet_idx, outlet_idx)) = infer_terminal_nodes_by_x(system) else {
        return Vec::new();
    };

    let mut graph = Graph::<usize, f64, Directed>::new();
    let mut graph_nodes = Vec::with_capacity(system.nodes.len());
    for idx in 0..system.nodes.len() {
        graph_nodes.push(graph.add_node(idx));
    }

    let box_h = system.box_dims.1.max(1e-9);
    let y_mid = system.box_dims.1 * 0.5;

    for channel in &system.channels {
        if channel.from_node >= graph_nodes.len() || channel.to_node >= graph_nodes.len() {
            continue;
        }

        let p_from = system.nodes[channel.from_node].point;
        let p_to = system.nodes[channel.to_node].point;
        let edge_len = segment_length(p_from, p_to).max(1e-9);
        let y_penalty = (((p_from.1 + p_to.1) * 0.5 - y_mid).abs() / box_h).clamp(0.0, 2.0);

        let forward_penalty = if p_to.0 + 1e-9 < p_from.0 {
            edge_len * 2.5
        } else {
            0.0
        };
        let reverse_penalty = if p_from.0 + 1e-9 < p_to.0 {
            edge_len * 2.5
        } else {
            0.0
        };

        let forward_weight = edge_len * (1.0 + y_penalty) + forward_penalty;
        let reverse_weight = edge_len * (1.0 + y_penalty) + reverse_penalty;

        graph.add_edge(
            graph_nodes[channel.from_node],
            graph_nodes[channel.to_node],
            forward_weight,
        );
        graph.add_edge(
            graph_nodes[channel.to_node],
            graph_nodes[channel.from_node],
            reverse_weight,
        );
    }

    let start = graph_nodes[inlet_idx];
    let goal = graph_nodes[outlet_idx];
    if let Some((_cost, path)) = astar(
        &graph,
        start,
        |node| node == goal,
        |edge| *edge.weight(),
        |_| 0.0,
    ) {
        return path
            .iter()
            .filter_map(|graph_idx| {
                let node_idx = graph[*graph_idx];
                system.nodes.get(node_idx).map(|node| node.point)
            })
            .collect();
    }

    vec![
        system.nodes[inlet_idx].point,
        system.nodes[outlet_idx].point,
    ]
}

/// Project evenly spaced markers along a path, constrained to a zone in X.
#[must_use]
pub fn project_markers_along_path(
    path: &[Point2D],
    marker_count: usize,
    zone_x_window: (f64, f64),
) -> Vec<Point2D> {
    if marker_count == 0 || path.len() < 2 {
        return Vec::new();
    }

    let cumulative = cumulative_lengths(path);
    let total_len = *cumulative.last().unwrap_or(&0.0);
    if total_len <= 0.0 {
        return vec![path[0]; marker_count];
    }

    let (zone_start_x, zone_end_x) = if zone_x_window.0 <= zone_x_window.1 {
        zone_x_window
    } else {
        (zone_x_window.1, zone_x_window.0)
    };

    let mut zone_s_min = f64::INFINITY;
    let mut zone_s_max = f64::NEG_INFINITY;

    for (seg_idx, segment) in path.windows(2).enumerate() {
        let p0 = segment[0];
        let p1 = segment[1];
        let seg_dx = p1.0 - p0.0;
        let seg_len = segment_length(p0, p1);
        if seg_len <= 0.0 {
            continue;
        }

        let (t0, t1) = if seg_dx.abs() <= 1e-12 {
            if p0.0 >= zone_start_x && p0.0 <= zone_end_x {
                (0.0, 1.0)
            } else {
                continue;
            }
        } else {
            let tx0 = (zone_start_x - p0.0) / seg_dx;
            let tx1 = (zone_end_x - p0.0) / seg_dx;
            let t_lo = tx0.min(tx1).clamp(0.0, 1.0);
            let t_hi = tx0.max(tx1).clamp(0.0, 1.0);
            if t_hi <= t_lo {
                continue;
            }
            (t_lo, t_hi)
        };

        let s0 = cumulative[seg_idx] + seg_len * t0;
        let s1 = cumulative[seg_idx] + seg_len * t1;
        zone_s_min = zone_s_min.min(s0);
        zone_s_max = zone_s_max.max(s1);
    }

    let (window_start_s, window_end_s) =
        if zone_s_min.is_finite() && zone_s_max.is_finite() && zone_s_max > zone_s_min {
            (zone_s_min, zone_s_max)
        } else {
            (0.0, total_len)
        };

    let window_len = (window_end_s - window_start_s).max(1e-12);
    let spacing = window_len / (marker_count as f64 + 1.0);

    (0..marker_count)
        .map(|idx| {
            let s =
                (window_start_s + spacing * (idx as f64 + 1.0)).clamp(window_start_s, window_end_s);
            point_at_arclength(path, &cumulative, s)
        })
        .collect()
}

fn channel_venturi_spec(channel: &ChannelSpec) -> Option<&ChannelVenturiSpec> {
    channel.metadata.as_ref()?.get::<ChannelVenturiSpec>()
}

fn segment_length(p0: Point2D, p1: Point2D) -> f64 {
    (p1.0 - p0.0).hypot(p1.1 - p0.1)
}

fn cumulative_lengths(path: &[Point2D]) -> Vec<f64> {
    let mut cumulative = Vec::with_capacity(path.len());
    cumulative.push(0.0);

    for segment in path.windows(2) {
        let next =
            cumulative.last().copied().unwrap_or(0.0) + segment_length(segment[0], segment[1]);
        cumulative.push(next);
    }

    cumulative
}

fn point_at_arclength(path: &[Point2D], cumulative: &[f64], target_s: f64) -> Point2D {
    if target_s <= 0.0 {
        return path[0];
    }

    let total_len = cumulative.last().copied().unwrap_or(0.0);
    if target_s >= total_len {
        return path[path.len() - 1];
    }

    let mut lo = 0usize;
    while lo + 1 < cumulative.len() && cumulative[lo + 1] < target_s {
        lo += 1;
    }

    let seg_start = path[lo];
    let seg_end = path[lo + 1];
    let s0 = cumulative[lo];
    let s1 = cumulative[lo + 1];
    let seg_len = (s1 - s0).max(1e-12);
    let t = (target_s - s0) / seg_len;

    (
        seg_start.0 + (seg_end.0 - seg_start.0) * t,
        seg_start.1 + (seg_end.1 - seg_start.1) * t,
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::model::ChannelSpec;
    use crate::geometry::{Channel, ChannelType, Node};

    fn simple_role_system() -> ChannelSystem {
        let nodes = vec![
            Node {
                id: 0,
                point: (0.0, 5.0),
                metadata: None,
            },
            Node {
                id: 1,
                point: (1.0, 5.0),
                metadata: None,
            },
            Node {
                id: 2,
                point: (2.0, 5.0),
                metadata: None,
            },
            Node {
                id: 3,
                point: (3.0, 6.0),
                metadata: None,
            },
            Node {
                id: 4,
                point: (3.0, 4.0),
                metadata: None,
            },
            Node {
                id: 5,
                point: (4.0, 5.0),
                metadata: None,
            },
            Node {
                id: 6,
                point: (5.0, 5.0),
                metadata: None,
            },
        ];

        let channels = vec![
            Channel {
                id: 0,
                from_node: 0,
                to_node: 1,
                width: 1.0,
                height: 1.0,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
            Channel {
                id: 1,
                from_node: 1,
                to_node: 2,
                width: 1.0,
                height: 1.0,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
            Channel {
                id: 2,
                from_node: 2,
                to_node: 3,
                width: 1.0,
                height: 1.0,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
            Channel {
                id: 3,
                from_node: 2,
                to_node: 4,
                width: 1.0,
                height: 1.0,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
            Channel {
                id: 4,
                from_node: 3,
                to_node: 5,
                width: 1.0,
                height: 1.0,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
            Channel {
                id: 5,
                from_node: 4,
                to_node: 5,
                width: 1.0,
                height: 1.0,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
            Channel {
                id: 6,
                from_node: 5,
                to_node: 6,
                width: 1.0,
                height: 1.0,
                channel_type: ChannelType::Straight,
                metadata: None,
            },
        ];

        ChannelSystem {
            box_dims: (6.0, 10.0),
            nodes,
            channels,
            box_outline: Vec::new(),
        }
    }

    #[test]
    fn classify_node_roles_assigns_expected_degree_roles() {
        let system = simple_role_system();
        let roles = classify_node_roles(&system);

        assert_eq!(roles.get(&0), Some(&MarkerRole::Inlet));
        assert_eq!(roles.get(&6), Some(&MarkerRole::Outlet));
        assert_eq!(roles.get(&2), Some(&MarkerRole::Split));
        assert_eq!(roles.get(&5), Some(&MarkerRole::Merge));
        assert_eq!(roles.get(&1), Some(&MarkerRole::Internal));
    }

    #[test]
    fn throat_count_uses_ctc_stream_metadata() {
        let mut blueprint = NetworkBlueprint::new("test");

        let c1 = ChannelSpec::new_pipe_rect("c1", "n0", "n1", 1.0, 1.0, 1.0, 1.0, 0.0)
            .with_metadata(ChannelVenturiSpec {
                n_throats: 2,
                is_ctc_stream: true,
                throat_width_m: 40e-6,
                height_m: 50e-6,
                inter_throat_spacing_m: 1.0e-3,
            });

        let c2 = ChannelSpec::new_pipe_rect("c2", "n1", "n2", 1.0, 1.0, 1.0, 1.0, 0.0)
            .with_metadata(ChannelVenturiSpec {
                n_throats: 3,
                is_ctc_stream: false,
                throat_width_m: 40e-6,
                height_m: 50e-6,
                inter_throat_spacing_m: 1.0e-3,
            });

        blueprint.add_channel(c1);
        blueprint.add_channel(c2);

        assert_eq!(throat_count_from_blueprint_metadata(&blueprint), 2);
    }

    #[test]
    fn projected_markers_are_monotonic_and_within_zone() {
        let path = vec![(0.0, 0.0), (100.0, 0.0)];
        let markers = project_markers_along_path(&path, 4, (35.0, 65.0));

        assert_eq!(markers.len(), 4);
        for point in &markers {
            assert!(point.0 > 35.0 && point.0 < 65.0);
        }
        for pair in markers.windows(2) {
            assert!(pair[1].0 > pair[0].0);
        }
    }
}
