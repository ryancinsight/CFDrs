//! `ChannelSystem::to_blueprint` — conversion to solver/render `NetworkBlueprint`.

use crate::error::{GeometryError, GeometryResult};
use crate::domain::therapy_metadata::TherapyZoneMetadata;
use crate::geometry::metadata::{
    ChannelPathMetadata, ChannelVisualRole, JunctionFamily, JunctionGeometryMetadata,
    NodeLayoutMetadata, VenturiGeometryMetadata,
};
use std::collections::HashSet;

use super::channel_system::ChannelSystem;
use super::{centerline_for_channel, polyline_length, ChannelType};

impl ChannelSystem {
    /// Convert this 2D schematic layout into a [`NetworkBlueprint`] for the 1D solver.
    ///
    /// This conversion preserves explicit node positions, channel polylines,
    /// inferred junction geometry, and venturi geometry so the resulting
    /// blueprint can act as the single source of truth for both solver and
    /// schematic rendering paths.
    ///
    /// `scale_m_per_unit` converts schematic coordinate units to SI metres
    /// (e.g. `1e-3` if coordinates are in millimetres, `1e-6` for micrometres).
    ///
    /// # Node classification
    /// - **Inlet** : no incoming channels (source node)
    /// - **Outlet**: no outgoing channels (sink node)
    /// - **Junction**: both incoming and outgoing channels
    ///
    /// # Resistance initialisation
    /// Each `ChannelSpec` is populated with an analytical Hagen-Poiseuille
    /// resistance as a linearisation seed (μ = 0.001 Pa·s, rectangular thin-slit
    /// approximation).  The 1D solver will refine these coefficients with the
    /// actual fluid model on its first Newton step.
    ///
    /// # Errors
    /// Returns [`GeometryError`] when the system has no channels, contains an
    /// out-of-bounds node reference, or produces a network without at least one
    /// inlet and one outlet.
    pub fn to_blueprint(
        &self,
        scale_m_per_unit: f64,
    ) -> GeometryResult<crate::domain::model::NetworkBlueprint> {
        use crate::domain::model::{
            ChannelShape, ChannelSpec, CrossSectionSpec, EdgeId, EdgeKind, NetworkBlueprint,
            NodeId, NodeKind, NodeSpec,
        };

        if self.channels.is_empty() {
            return Err(GeometryError::InvalidChannelPath {
                reason: "Cannot build NetworkBlueprint from an empty ChannelSystem".to_string(),
            });
        }

        // ── topology classification ───────────────────────────────────────
        let max_idx = self.nodes.len();
        for ch in &self.channels {
            if ch.from_node >= max_idx || ch.to_node >= max_idx {
                return Err(GeometryError::ChannelCreationFailed {
                    from_id: ch.from_node,
                    to_id: ch.to_node,
                    reason: format!(
                        "Channel {} references a node index outside the node list (len={})",
                        ch.id, max_idx
                    ),
                });
            }
        }

        let mut has_incoming: HashSet<usize> = HashSet::new();
        let mut has_outgoing: HashSet<usize> = HashSet::new();
        for ch in &self.channels {
            has_outgoing.insert(ch.from_node);
            has_incoming.insert(ch.to_node);
        }

        let classify = |id: usize| -> NodeKind {
            match (has_incoming.contains(&id), has_outgoing.contains(&id)) {
                (false, _) => NodeKind::Inlet,
                (true, false) => NodeKind::Outlet,
                (true, true) => NodeKind::Junction,
            }
        };

        // ── build blueprint ───────────────────────────────────────────────
        let mut blueprint = NetworkBlueprint::new(format!(
            "from_channel_system_{}x{}",
            self.box_dims.0 as u64, self.box_dims.1 as u64
        ));

        let has_inlet = self
            .nodes
            .iter()
            .any(|n| matches!(classify(n.id), NodeKind::Inlet));
        let has_outlet = self
            .nodes
            .iter()
            .any(|n| matches!(classify(n.id), NodeKind::Outlet));
        if !has_inlet {
            return Err(GeometryError::InvalidChannelPath {
                reason: "ChannelSystem has no inlet nodes (every node has incoming channels)"
                    .to_string(),
            });
        }
        if !has_outlet {
            return Err(GeometryError::InvalidChannelPath {
                reason: "ChannelSystem has no outlet nodes (every node has outgoing channels)"
                    .to_string(),
            });
        }

        let node_bp_ids: Vec<String> = self
            .nodes
            .iter()
            .map(|node| {
                node.name
                    .clone()
                    .unwrap_or_else(|| format!("node_{}", node.id))
            })
            .collect();

        for node in &self.nodes {
            let mut node_spec = NodeSpec::new(
                node.name
                    .clone()
                    .unwrap_or_else(|| format!("node_{}", node.id)),
                node.kind.unwrap_or_else(|| classify(node.id)),
            )
                .with_layout(NodeLayoutMetadata {
                    x_mm: node.point.0,
                    y_mm: node.point.1,
                });
            if let Some(junction_geometry) = node.junction_geometry.clone() {
                node_spec = node_spec.with_junction_geometry(junction_geometry);
            }
            node_spec.metadata = node.metadata.clone();
            blueprint.add_node(node_spec);
        }

        // ── channels → ChannelSpec ────────────────────────────────────────
        for channel in &self.channels {
            let centerline = centerline_for_channel(channel, &self.nodes);
            let length_m = channel
                .physical_length_m
                .unwrap_or_else(|| polyline_length(&centerline) * scale_m_per_unit);

            if length_m <= 0.0 {
                return Err(GeometryError::ChannelCreationFailed {
                    from_id: channel.from_node,
                    to_id: channel.to_node,
                    reason: format!("Channel {} has zero or negative path length", channel.id),
                });
            }

            let (avg_width, height) = match &channel.channel_type {
                ChannelType::Frustum {
                    inlet_width,
                    throat_width,
                    outlet_width,
                    ..
                } => (
                    (inlet_width + throat_width + outlet_width) / 3.0,
                    channel.height,
                ),
                _ => (channel.width, channel.height),
            };
            let width_m = channel
                .physical_width_m
                .unwrap_or(avg_width * scale_m_per_unit);
            let height_m = channel
                .physical_height_m
                .unwrap_or(height * scale_m_per_unit);

            // Hagen-Poiseuille thin-slit approximation: R ≈ 12μL / (wh³)
            let mu = 0.001_f64; // Pa·s  (water at 20 °C)
            let resistance =
                (12.0 * mu * length_m / (width_m * height_m.powi(3))).max(f64::EPSILON);

            let cross_section = CrossSectionSpec::Rectangular { width_m, height_m };

            // Detect serpentine channels and estimate bend geometry from path.
            let channel_shape = channel.physical_shape.unwrap_or_else(|| {
                if matches!(channel.channel_type, ChannelType::Serpentine { .. }) {
                    let segments = count_serpentine_segments(&centerline);
                    let bend_radius_m = estimate_bend_radius(&centerline, scale_m_per_unit)
                        .unwrap_or(width_m * 0.5);
                    ChannelShape::Serpentine {
                        segments,
                        bend_radius_m,
                    }
                } else {
                    ChannelShape::Straight
                }
            });

            let mut channel_spec = ChannelSpec {
                id: EdgeId::new(
                    channel
                        .name
                        .clone()
                        .unwrap_or_else(|| format!("ch_{}", channel.id)),
                ),
                kind: EdgeKind::Pipe,
                from: NodeId::new(node_bp_ids[channel.from_node].clone()),
                to: NodeId::new(node_bp_ids[channel.to_node].clone()),
                length_m,
                cross_section,
                channel_shape,
                resistance,
                quad_coeff: 0.0,
                valve_cv: None,
                pump_max_flow: None,
                pump_max_pressure: None,
                path: Some(ChannelPathMetadata {
                    polyline_mm: centerline.clone(),
                    visual_role: channel.visual_role.unwrap_or(ChannelVisualRole::InternalLink),
                }),
                venturi_geometry: channel.venturi_geometry.clone().or_else(|| {
                    venturi_geometry_from_channel_type(
                        &channel.channel_type,
                        height_m,
                        scale_m_per_unit,
                    )
                }),
                metadata: channel.metadata.clone(),
            };
            if let Some(zone) = &channel.therapy_zone {
                channel_spec = channel_spec.with_metadata(TherapyZoneMetadata::new(zone.clone()));
            }
            if let Some(venturi) = channel_spec.venturi_geometry.clone() {
                channel_spec = channel_spec.with_metadata(venturi);
            }
            blueprint.add_channel(channel_spec);
        }

        let inferred_junctions: Vec<(String, JunctionGeometryMetadata)> = blueprint
            .nodes
            .iter()
            .filter(|node| matches!(node.kind, NodeKind::Junction))
            .filter_map(|node| {
                infer_junction_geometry(&blueprint, node.id.as_str())
                    .map(|geometry| (node.id.as_str().to_string(), geometry))
            })
            .collect();
        for (node_id, geometry) in inferred_junctions {
            let Some(node) = blueprint
                .nodes
                .iter_mut()
                .find(|node| node.id.as_str() == node_id)
            else {
                continue;
            };
            if node.junction_geometry.is_none() {
                node.junction_geometry = Some(geometry.clone());
                node.metadata
                    .get_or_insert_with(crate::geometry::metadata::MetadataContainer::new)
                    .insert(geometry);
            }
        }

        Ok(blueprint)
    }
}

fn venturi_geometry_from_channel_type(
    channel_type: &ChannelType,
    height_m: f64,
    scale_m_per_unit: f64,
) -> Option<VenturiGeometryMetadata> {
    let ChannelType::Frustum {
        inlet_width,
        throat_width,
        outlet_width,
        path,
        ..
    } = channel_type
    else {
        return None;
    };
    let length_m = polyline_length(path) * scale_m_per_unit;
    Some(VenturiGeometryMetadata {
        throat_width_m: throat_width * scale_m_per_unit,
        throat_height_m: height_m,
        throat_length_m: length_m,
        inlet_width_m: inlet_width * scale_m_per_unit,
        outlet_width_m: outlet_width * scale_m_per_unit,
        convergent_half_angle_deg: 15.0,
        divergent_half_angle_deg: 15.0,
    })
}

fn infer_junction_geometry(
    blueprint: &crate::domain::model::NetworkBlueprint,
    node_id: &str,
) -> Option<JunctionGeometryMetadata> {
    let outgoing_angles: Vec<f64> = blueprint
        .channels
        .iter()
        .filter(|channel| channel.from.as_str() == node_id)
        .filter_map(|channel| channel.path.as_ref())
        .filter_map(|path| direction_angle_deg(&path.polyline_mm, true))
        .collect();
    let incoming_angles: Vec<f64> = blueprint
        .channels
        .iter()
        .filter(|channel| channel.to.as_str() == node_id)
        .filter_map(|channel| channel.path.as_ref())
        .filter_map(|path| direction_angle_deg(&path.polyline_mm, false))
        .collect();

    let family = match (incoming_angles.len(), outgoing_angles.len()) {
        (1, 2) => JunctionFamily::Bifurcation,
        (1, 3) => JunctionFamily::Trifurcation,
        (2, 1) => JunctionFamily::Merge,
        (3, 1) => JunctionFamily::Merge,
        (1, 1) => JunctionFamily::Merge,
        (_, _) if incoming_angles.len() + outgoing_angles.len() >= 4 => JunctionFamily::Cross,
        _ => return None,
    };

    Some(JunctionGeometryMetadata {
        junction_family: family,
        branch_angles_deg: centered_angles(&outgoing_angles),
        merge_angles_deg: centered_angles(&incoming_angles),
    })
}

fn centered_angles(angles_deg: &[f64]) -> Vec<f64> {
    if angles_deg.is_empty() {
        return Vec::new();
    }
    let mut sorted = angles_deg.to_vec();
    sorted.sort_by(|a, b| a.total_cmp(b));
    let center = sorted[sorted.len() / 2];
    sorted.into_iter().map(|angle| angle - center).collect()
}

fn direction_angle_deg(path: &[(f64, f64)], from_start: bool) -> Option<f64> {
    let (p1, p2) = if from_start {
        path.windows(2)
            .find(|segment| {
                let dx = segment[1].0 - segment[0].0;
                let dy = segment[1].1 - segment[0].1;
                dx.abs() > 1e-9 || dy.abs() > 1e-9
            })
            .map(|segment| (segment[0], segment[1]))?
    } else {
        path.windows(2)
            .rev()
            .find(|segment| {
                let dx = segment[1].0 - segment[0].0;
                let dy = segment[1].1 - segment[0].1;
                dx.abs() > 1e-9 || dy.abs() > 1e-9
            })
            .map(|segment| (segment[1], segment[0]))?
    };
    Some((p2.1 - p1.1).atan2(p2.0 - p1.0).to_degrees())
}

/// Count the number of straight segments in a serpentine path by detecting
/// y-direction reversals.
///
/// A serpentine alternates up/down, so each reversal in y-coordinate direction
/// indicates a new straight segment. Returns at least 1 for a non-empty path.
fn count_serpentine_segments(path: &[(f64, f64)]) -> usize {
    if path.len() < 3 {
        return 1;
    }
    let mut reversals = 0usize;
    let mut prev_dy: f64 = 0.0;
    for i in 1..path.len() {
        let dy = path[i].1 - path[i - 1].1;
        if dy.abs() > 1e-9 {
            if prev_dy != 0.0 && dy.signum() != prev_dy.signum() {
                reversals += 1;
            }
            prev_dy = dy;
        }
    }
    // N reversals → N+1 segments
    (reversals + 1).max(1)
}

/// Estimate bend radius from a serpentine path by measuring the curvature at
/// y-direction reversal points. Returns `None` if the path has no clear bends.
fn estimate_bend_radius(path: &[(f64, f64)], scale: f64) -> Option<f64> {
    if path.len() < 5 {
        return None;
    }
    let mut radii = Vec::new();
    for i in 2..path.len() - 2 {
        let dy_before = path[i].1 - path[i - 1].1;
        let dy_after = path[i + 1].1 - path[i].1;
        // Detect reversal
        if dy_before.abs() > 1e-9
            && dy_after.abs() > 1e-9
            && dy_before.signum() != dy_after.signum()
        {
            // Circumscribed circle radius through three consecutive points
            let (x1, y1) = path[i - 1];
            let (x2, y2) = path[i];
            let (x3, y3) = path[i + 1];
            let a = ((x2 - x1).powi(2) + (y2 - y1).powi(2)).sqrt();
            let b = ((x3 - x2).powi(2) + (y3 - y2).powi(2)).sqrt();
            let c = ((x3 - x1).powi(2) + (y3 - y1).powi(2)).sqrt();
            let s = (a + b + c) / 2.0;
            let area = (s * (s - a) * (s - b) * (s - c)).max(0.0).sqrt();
            if area > 1e-12 {
                let r = (a * b * c) / (4.0 * area);
                radii.push(r * scale); // convert to metres
            }
        }
    }
    if radii.is_empty() {
        return None;
    }
    // Use the median radius for robustness
    radii.sort_by(|a, b| a.total_cmp(b));
    Some(radii[radii.len() / 2])
}
