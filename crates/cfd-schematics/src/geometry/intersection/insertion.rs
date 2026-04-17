use crate::domain::model::{ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::geometry::Point2D;

use super::detection::detect_crossings;
use super::IntersectionResult;

struct SplitInfo {
    junction_node_id: String,
    t: f64,
    seg: usize,
}

/// Insert junction nodes at all detected channel intersections.
pub fn insert_intersection_nodes(system: &mut NetworkBlueprint) -> IntersectionResult {
    let crossings = detect_crossings(system);
    if crossings.is_empty() {
        return IntersectionResult {
            intersection_count: 0,
            junction_node_ids: Vec::new(),
        };
    }

    let mut junction_ids = Vec::with_capacity(crossings.len());
    let mut new_nodes = Vec::new();
    let initial_node_count = system.nodes.len();

    for (crossing_index, crossing) in crossings.iter().enumerate() {
        let node_id = format!(
            "intersect_{}_{}_{}",
            crossing.channel_a, crossing.channel_b, crossing_index
        );

        system.nodes.push(NodeSpec::new_at(
            node_id.clone(),
            NodeKind::Junction,
            crossing.point,
        ));
        junction_ids.push(initial_node_count + crossing_index);
        new_nodes.push(node_id);
    }

    let mut splits_per_channel: std::collections::HashMap<usize, Vec<SplitInfo>> =
        std::collections::HashMap::new();

    for (crossing_index, crossing) in crossings.iter().enumerate() {
        splits_per_channel
            .entry(crossing.channel_a)
            .or_default()
            .push(SplitInfo {
                junction_node_id: new_nodes[crossing_index].clone(),
                t: crossing.t_a,
                seg: crossing.seg_a,
            });
        splits_per_channel
            .entry(crossing.channel_b)
            .or_default()
            .push(SplitInfo {
                junction_node_id: new_nodes[crossing_index].clone(),
                t: crossing.t_b,
                seg: crossing.seg_b,
            });
    }

    for splits in splits_per_channel.values_mut() {
        splits.sort_by(|left, right| {
            left.seg
                .cmp(&right.seg)
                .then(left.t.partial_cmp(&right.t).expect("structural invariant"))
        });
    }

    let mut next_channel_idx = 0;
    let mut new_channels: Vec<ChannelSpec> = Vec::new();

    for (channel_index, original_channel) in system.channels.iter().enumerate() {
        if let Some(splits) = splits_per_channel.get(&channel_index) {
            let centerline = &original_channel.path;
            if centerline.len() < 2 {
                new_channels.push(original_channel.clone());
                continue;
            }

            let mut waypoints: Vec<(String, Point2D)> = Vec::new();
            waypoints.push((original_channel.from.0.clone(), centerline[0]));

            let mut path_index = 0;
            let mut split_iter = splits.iter().peekable();

            while path_index < centerline.len() - 1 {
                while let Some(split) = split_iter.peek() {
                    if split.seg == path_index {
                        let current_split = split_iter.next().expect("structural invariant");
                        let point = (
                            centerline[path_index]
                                .0
                                .mul_add(1.0 - current_split.t, centerline[path_index + 1].0 * current_split.t),
                            centerline[path_index]
                                .1
                                .mul_add(1.0 - current_split.t, centerline[path_index + 1].1 * current_split.t),
                        );
                        waypoints.push((current_split.junction_node_id.clone(), point));
                    } else {
                        break;
                    }
                }
                path_index += 1;
            }

            waypoints.push((
                original_channel.to.0.clone(),
                *centerline.last().expect("structural invariant"),
            ));

            for pair in waypoints.windows(2) {
                let (from_id, from_point) = &pair[0];
                let (to_id, to_point) = &pair[1];
                next_channel_idx += 1;
                let mut channel_clone = original_channel.clone();
                channel_clone.id = crate::domain::model::EdgeId(format!(
                    "{}_part{}",
                    original_channel.id.0, next_channel_idx
                ));
                channel_clone.from = crate::domain::model::NodeId(from_id.clone());
                channel_clone.to = crate::domain::model::NodeId(to_id.clone());
                channel_clone.path = vec![*from_point, *to_point];
                new_channels.push(channel_clone);
            }
        } else {
            new_channels.push(original_channel.clone());
        }
    }

    system.channels = new_channels;

    IntersectionResult {
        intersection_count: crossings.len(),
        junction_node_ids: junction_ids,
    }
}