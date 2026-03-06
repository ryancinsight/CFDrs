use crate::domain::model::{ChannelShape, CrossSectionSpec, NetworkBlueprint};
use crate::error::{VisualizationError, VisualizationResult};
use crate::geometry::metadata::{
    BlueprintRenderHints, ChannelPathMetadata, NodeLayoutMetadata, VenturiGeometryMetadata,
};
use crate::geometry::types::centerline_for_channel;
use crate::geometry::{Channel, ChannelSystem, ChannelType, Node, Point2D};
use crate::visualizations::annotations::SchematicAnnotations;
use crate::visualizations::plotters_backend::PlottersRenderer;
use crate::visualizations::traits::{RenderConfig, SchematicRenderer};
use std::collections::{HashMap, VecDeque};

/// Plot a channel system using the default renderer and configuration
///
/// This function provides backward compatibility while using the new
/// abstracted visualization system internally.
pub fn plot_geometry(system: &ChannelSystem, output_path: &str) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let config = RenderConfig::default();
    renderer.render_system(system, output_path, &config)
}

/// Plot a channel system with custom configuration
///
/// This function allows for more control over the rendering process
/// by accepting a custom configuration.
pub fn plot_geometry_with_config(
    system: &ChannelSystem,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    renderer.render_system(system, output_path, config)
}

/// Plot a channel system with explicit annotation payload.
///
/// This helper leaves the base config unchanged and overlays caller-provided
/// markers/labels through `RenderConfig.annotations` for this render only.
pub fn plot_geometry_with_annotations(
    system: &ChannelSystem,
    output_path: &str,
    config: &RenderConfig,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let mut annotated_config = config.clone();
    annotated_config.annotations = Some(annotations.clone());
    renderer.render_system(system, output_path, &annotated_config)
}

/// Convert a [`NetworkBlueprint`] directly into a renderable [`ChannelSystem`].
///
/// This is the single-source-of-truth schematic path for designs that already
/// exist as blueprints. Explicit layout metadata is used when present; nodes
/// without layout metadata fall back to a deterministic DAG layout.
pub fn channel_system_from_blueprint(
    blueprint: &NetworkBlueprint,
    box_dims_hint: Option<(f64, f64)>,
) -> VisualizationResult<ChannelSystem> {
    let box_dims = box_dims_hint.unwrap_or((127.76, 85.47));
    let auto_layout = auto_layout_positions(blueprint, box_dims);

    let mut node_indices = HashMap::with_capacity(blueprint.nodes.len());
    let mut nodes = Vec::with_capacity(blueprint.nodes.len());
    for (idx, node_spec) in blueprint.nodes.iter().enumerate() {
        let point = node_spec
            .layout
            .map(|layout| (layout.x_mm, layout.y_mm))
            .or_else(|| {
                node_spec
                    .metadata
                    .as_ref()
                    .and_then(|meta| meta.get::<NodeLayoutMetadata>())
                    .map(|layout| (layout.x_mm, layout.y_mm))
            })
            .unwrap_or_else(|| {
                auto_layout
                    .get(node_spec.id.as_str())
                    .copied()
                    .unwrap_or((box_dims.0 * 0.5, box_dims.1 * 0.5))
            });
        nodes.push(Node {
            id: idx,
            name: Some(node_spec.id.as_str().to_string()),
            point,
            kind: Some(node_spec.kind),
            junction_geometry: node_spec.junction_geometry.clone(),
            metadata: node_spec.metadata.clone(),
        });
        node_indices.insert(node_spec.id.as_str().to_string(), idx);
    }

    let mut channels = Vec::with_capacity(blueprint.channels.len());
    for (idx, channel_spec) in blueprint.channels.iter().enumerate() {
        let from_idx = *node_indices.get(channel_spec.from.as_str()).ok_or_else(|| {
            VisualizationError::CoordinateTransformError {
                message: format!(
                    "blueprint channel '{}' references missing from-node '{}'",
                    channel_spec.id.as_str(),
                    channel_spec.from.as_str()
                ),
            }
        })?;
        let to_idx = *node_indices.get(channel_spec.to.as_str()).ok_or_else(|| {
            VisualizationError::CoordinateTransformError {
                message: format!(
                    "blueprint channel '{}' references missing to-node '{}'",
                    channel_spec.id.as_str(),
                    channel_spec.to.as_str()
                ),
            }
        })?;

        let (width_mm, height_mm) = match channel_spec.cross_section {
            CrossSectionSpec::Circular { diameter_m } => (diameter_m * 1e3, diameter_m * 1e3),
            CrossSectionSpec::Rectangular { width_m, height_m } => (width_m * 1e3, height_m * 1e3),
        };
        let path = explicit_or_generated_path(channel_spec, &nodes[from_idx], &nodes[to_idx]);
        let channel_type = channel_type_from_blueprint(channel_spec, &path);
        channels.push(Channel {
            id: idx,
            name: Some(channel_spec.id.as_str().to_string()),
            from_node: from_idx,
            to_node: to_idx,
            width: width_mm,
            height: height_mm,
            channel_type,
            visual_role: channel_spec.path.as_ref().map(|path| path.visual_role),
            physical_length_m: Some(channel_spec.length_m),
            physical_width_m: Some(width_mm * 1e-3),
            physical_height_m: Some(height_mm * 1e-3),
            physical_shape: Some(channel_spec.channel_shape),
            therapy_zone: channel_spec
                .metadata
                .as_ref()
                .and_then(|meta| meta.get::<crate::domain::therapy_metadata::TherapyZoneMetadata>())
                .map(|meta| meta.zone.clone()),
            venturi_geometry: channel_spec.venturi_geometry.clone(),
            metadata: channel_spec.metadata.clone(),
        });
    }

    let box_outline = blueprint_box_outline(blueprint, &nodes, box_dims);

    let system = ChannelSystem {
        box_dims,
        nodes,
        channels,
        box_outline,
    };
    system
        .validate()
        .map_err(|err| VisualizationError::CoordinateTransformError {
            message: format!("blueprint render conversion failed validation: {err}"),
        })?;
    Ok(system)
}

fn blueprint_box_outline(
    blueprint: &NetworkBlueprint,
    nodes: &[Node],
    box_dims: (f64, f64),
) -> Vec<(Point2D, Point2D)> {
    let mut outline = vec![
        ((0.0, 0.0), (box_dims.0, 0.0)),
        ((box_dims.0, 0.0), (box_dims.0, box_dims.1)),
        ((box_dims.0, box_dims.1), (0.0, box_dims.1)),
        ((0.0, box_dims.1), (0.0, 0.0)),
    ];

    for (idx, node_spec) in blueprint.nodes.iter().enumerate() {
        let Some(node) = nodes.get(idx) else {
            continue;
        };

        match node_spec.kind {
            crate::domain::model::NodeKind::Inlet => {
                outline.push(((0.0, node.point.1), node.point));
            }
            crate::domain::model::NodeKind::Outlet => {
                outline.push((node.point, (box_dims.0, node.point.1)));
            }
            crate::domain::model::NodeKind::Reservoir | crate::domain::model::NodeKind::Junction => {
            }
        }
    }

    outline
}

/// Plot a blueprint using blueprint-native layout metadata.
pub fn plot_blueprint(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let system = channel_system_from_blueprint(blueprint, Some((127.76, 85.47)))?;
    plot_geometry_with_config(&system, output_path, config)
}

/// Plot a blueprint using blueprint-native layout metadata and explicit annotations.
pub fn plot_blueprint_with_annotations(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    let system = channel_system_from_blueprint(blueprint, Some((127.76, 85.47)))?;
    plot_geometry_with_annotations(&system, output_path, config, annotations)
}

/// Render a blueprint to PNG with annotations derived entirely from the
/// blueprint's own metadata — no caller-supplied [`SchematicAnnotations`]
/// required.
///
/// Annotations are built by reading [`BlueprintRenderHints`] from the
/// blueprint and running the standard node-role classification + throat
/// projection logic that lives in `cfd-schematics`.  Callers (e.g.
/// `cfd-optim`) need only populate `BlueprintRenderHints` in
/// [`NetworkBlueprint::with_render_hints`] to get fully-annotated figures.
pub fn plot_blueprint_auto_annotated(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let system = channel_system_from_blueprint(blueprint, Some((127.76, 85.47)))?;
    let annotations = build_auto_annotations(blueprint, &system);
    plot_geometry_with_annotations(&system, output_path, config, &annotations)
}

/// Build [`SchematicAnnotations`] from blueprint metadata and geometric
/// analysis.  This consolidates all annotation-building logic so that any
/// renderer calling [`plot_blueprint_auto_annotated`] produces identical
/// figures without domain-specific render code in the caller.
fn build_auto_annotations(
    blueprint: &NetworkBlueprint,
    system: &ChannelSystem,
) -> SchematicAnnotations {
    use crate::visualizations::annotations::{
        center_biased_main_path, classify_node_roles, project_markers_along_path,
        therapy_zone_presence, throat_count_from_blueprint_metadata, AnnotationMarker, MarkerRole,
    };
    use std::collections::BTreeSet;

    let mut annotations = SchematicAnnotations::report_default();
    let mut roles = classify_node_roles(system);
    let mut seen_points: BTreeSet<(i64, i64)> = system
        .nodes
        .iter()
        .map(|n| blueprint_quantize(n.point))
        .collect();

    // ── Therapy-zone role overrides ──────────────────────────────────────────
    let (has_target, has_bypass) = therapy_zone_presence(blueprint);
    if has_target || has_bypass {
        let y_mid = system.box_dims.1 * 0.5;
        let tx_min = system.box_dims.0 * 0.35;
        let tx_max = system.box_dims.0 * 0.65;
        for (node_idx, node) in system.nodes.iter().enumerate() {
            let base = roles
                .get(&node_idx)
                .copied()
                .unwrap_or(MarkerRole::Internal);
            if matches!(
                base,
                MarkerRole::Inlet
                    | MarkerRole::Outlet
                    | MarkerRole::Split
                    | MarkerRole::Merge
            ) {
                continue;
            }
            let near_cl = (node.point.1 - y_mid).abs() <= system.box_dims.1 * 0.18;
            let in_win = node.point.0 >= tx_min && node.point.0 <= tx_max;
            if has_target && near_cl && in_win {
                roles.insert(node_idx, MarkerRole::TherapyTarget);
            } else if has_bypass && !near_cl {
                roles.insert(node_idx, MarkerRole::Bypass);
            }
        }
    }

    // ── Node markers ─────────────────────────────────────────────────────────
    let mut split_idx = 1usize;
    let mut merge_idx = 1usize;
    for (node_idx, node) in system.nodes.iter().enumerate() {
        let role = roles
            .get(&node_idx)
            .copied()
            .unwrap_or(MarkerRole::Internal);
        let marker = match role {
            MarkerRole::Inlet => AnnotationMarker::new(node.point, role).with_label("IN", true),
            MarkerRole::Outlet => AnnotationMarker::new(node.point, role).with_label("OUT", true),
            MarkerRole::Split => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("S{split_idx}"), true);
                split_idx += 1;
                m
            }
            MarkerRole::Merge => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("M{merge_idx}"), true);
                merge_idx += 1;
                m
            }
            _ => AnnotationMarker::new(node.point, role),
        };
        annotations.markers.push(marker);
    }

    // ── Interior centerline vertices ─────────────────────────────────────────
    for point in centerline_vertices(system) {
        let q = blueprint_quantize(point);
        if seen_points.contains(&q) {
            continue;
        }
        seen_points.insert(q);
        annotations
            .markers
            .push(AnnotationMarker::new(point, MarkerRole::Internal));
    }

    // ── Venturi throat markers ────────────────────────────────────────────────
    let hints: Option<&BlueprintRenderHints> = blueprint.render_hints();
    let mut throat_count = throat_count_from_blueprint_metadata(blueprint);
    if throat_count == 0 {
        if let Some(h) = hints {
            throat_count = h.throat_count_hint;
        }
    }
    if throat_count > 0 {
        let main_path = center_biased_main_path(system);
        let zone = (system.box_dims.0 * 0.35, system.box_dims.0 * 0.65);
        for (idx, point) in project_markers_along_path(&main_path, throat_count, zone)
            .iter()
            .enumerate()
        {
            annotations.markers.push(
                AnnotationMarker::new(*point, MarkerRole::VenturiThroat)
                    .with_label(format!("TH{}", idx + 1), true),
            );
        }
    }

    // ── Legend note from render hints ─────────────────────────────────────────
    if let Some(h) = hints {
        let min_w = system
            .channels
            .iter()
            .map(|ch| ch.width)
            .fold(f64::INFINITY, f64::min);
        let max_w = system
            .channels
            .iter()
            .map(|ch| ch.width)
            .fold(f64::NEG_INFINITY, f64::max);
        annotations.legend_note = Some(format!(
            "{}  |  seq: {}  |  layers: {}  |  throats: {}  |  width: {:.2}-{:.2} mm",
            h.treatment_label, h.stage_sequence, h.split_layers, throat_count, min_w, max_w
        ));
        annotations.markers.push(
            AnnotationMarker::new(
                (system.box_dims.0 * 0.50, system.box_dims.1 * 0.82),
                MarkerRole::Internal,
            )
            .with_label(
                format!(
                    "{}  |  {} split layers  |  {} treatment",
                    h.stage_sequence, h.split_layers, h.treatment_label
                ),
                true,
            ),
        );
    }

    annotations
}

fn blueprint_quantize(point: Point2D) -> (i64, i64) {
    (
        (point.0 * 1_000.0).round() as i64,
        (point.1 * 1_000.0).round() as i64,
    )
}

/// Render a [`ChannelSystem`] to an image with auto-generated annotations.
///
/// Node roles (inlet, outlet, split, merge) are inferred from graph
/// degree analysis.  Venturi throats are detected from channel metadata.
/// No [`NetworkBlueprint`] is required — this function operates on the
/// geometry-level data alone.
pub fn plot_geometry_auto_annotated(
    system: &ChannelSystem,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let annotations = build_system_annotations(system);
    plot_geometry_with_annotations(system, output_path, config, &annotations)
}

/// Build [`SchematicAnnotations`] from a [`ChannelSystem`] using graph-degree
/// classification — no blueprint required.
///
/// This mirrors the logic of [`build_auto_annotations`] but operates on the
/// [`ChannelSystem`] directly instead of a [`NetworkBlueprint`]:
///
/// - Node roles via graph degree (inlet/outlet/split/merge)
/// - Therapy-zone overrides from `Channel.therapy_zone`
/// - Interior centerline vertex markers
/// - Venturi throat markers from `Frustum` channel metadata
/// - Legend note with width range and topology summary
fn build_system_annotations(system: &ChannelSystem) -> SchematicAnnotations {
    use crate::visualizations::annotations::{
        classify_node_roles, AnnotationMarker, MarkerRole,
    };
    use std::collections::BTreeSet;

    let mut annotations = SchematicAnnotations::report_default();
    let mut roles = classify_node_roles(system);
    let mut seen_points: BTreeSet<(i64, i64)> = system
        .nodes
        .iter()
        .map(|n| blueprint_quantize(n.point))
        .collect();

    // ── Therapy-zone role overrides from Channel.therapy_zone ─────────────
    let has_target = system.channels.iter().any(|ch| {
        matches!(
            ch.therapy_zone,
            Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)
        )
    });
    let has_bypass = system.channels.iter().any(|ch| {
        matches!(
            ch.therapy_zone,
            Some(crate::domain::therapy_metadata::TherapyZone::HealthyBypass)
        )
    });
    if has_target || has_bypass {
        let y_mid = system.box_dims.1 * 0.5;
        let tx_min = system.box_dims.0 * 0.35;
        let tx_max = system.box_dims.0 * 0.65;
        for (node_idx, node) in system.nodes.iter().enumerate() {
            let base = roles
                .get(&node_idx)
                .copied()
                .unwrap_or(MarkerRole::Internal);
            if matches!(
                base,
                MarkerRole::Inlet
                    | MarkerRole::Outlet
                    | MarkerRole::Split
                    | MarkerRole::Merge
            ) {
                continue;
            }
            let near_cl = (node.point.1 - y_mid).abs() <= system.box_dims.1 * 0.18;
            let in_win = node.point.0 >= tx_min && node.point.0 <= tx_max;
            if has_target && near_cl && in_win {
                roles.insert(node_idx, MarkerRole::TherapyTarget);
            } else if has_bypass && !near_cl {
                roles.insert(node_idx, MarkerRole::Bypass);
            }
        }
    }

    // ── Node markers ─────────────────────────────────────────────────────
    let mut split_idx = 1usize;
    let mut merge_idx = 1usize;
    let mut inlet_idx = 1usize;
    let mut outlet_idx = 1usize;
    for (node_idx, node) in system.nodes.iter().enumerate() {
        let role = roles
            .get(&node_idx)
            .copied()
            .unwrap_or(MarkerRole::Internal);
        let marker = match role {
            MarkerRole::Inlet => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("IN{inlet_idx}"), true);
                inlet_idx += 1;
                m
            }
            MarkerRole::Outlet => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("OUT{outlet_idx}"), true);
                outlet_idx += 1;
                m
            }
            MarkerRole::Split => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("S{split_idx}"), true);
                split_idx += 1;
                m
            }
            MarkerRole::Merge => {
                let m = AnnotationMarker::new(node.point, role)
                    .with_label(format!("M{merge_idx}"), true);
                merge_idx += 1;
                m
            }
            _ => AnnotationMarker::new(node.point, role),
        };
        annotations.markers.push(marker);
    }

    // ── Interior centerline vertices ─────────────────────────────────────
    for point in centerline_vertices(system) {
        let q = blueprint_quantize(point);
        if seen_points.contains(&q) {
            continue;
        }
        seen_points.insert(q);
        annotations
            .markers
            .push(AnnotationMarker::new(point, MarkerRole::Internal));
    }

    // ── Venturi throat markers ────────────────────────────────────────────
    let mut throat_idx = 1usize;
    for channel in &system.channels {
        let is_venturi = match &channel.channel_type {
            crate::geometry::ChannelType::Frustum {
                has_venturi_throat, ..
            } => *has_venturi_throat,
            _ => false,
        };
        if !is_venturi {
            continue;
        }
        if let crate::geometry::ChannelType::Frustum {
            path,
            throat_position,
            ..
        } = &channel.channel_type
        {
            if path.len() >= 2 {
                let idx_f = *throat_position * (path.len() - 1) as f64;
                let lo = idx_f.floor() as usize;
                let hi = (lo + 1).min(path.len() - 1);
                let frac = idx_f - lo as f64;
                let throat_pt = (
                    path[lo].0 + (path[hi].0 - path[lo].0) * frac,
                    path[lo].1 + (path[hi].1 - path[lo].1) * frac,
                );
                annotations.markers.push(
                    AnnotationMarker::new(throat_pt, MarkerRole::VenturiThroat)
                        .with_label(format!("TH{throat_idx}"), true),
                );
                throat_idx += 1;
            }
        }
    }

    // ── Legend note ───────────────────────────────────────────────────────
    let split_count = split_idx - 1;
    let merge_count = merge_idx - 1;
    let throat_count = throat_idx - 1;
    let min_w = system
        .channels
        .iter()
        .map(|ch| ch.width)
        .fold(f64::INFINITY, f64::min);
    let max_w = system
        .channels
        .iter()
        .map(|ch| ch.width)
        .fold(f64::NEG_INFINITY, f64::max);
    if !system.channels.is_empty() {
        annotations.legend_note = Some(format!(
            "splits: {split_count} | merges: {merge_count} | throats: {throat_count} | width: {min_w:.2}\u{2013}{max_w:.2} mm",
        ));
    }

    annotations
}

/// Plot a channel system using a custom renderer
///
/// This function demonstrates the flexibility of the abstracted system
/// by allowing any renderer that implements `SchematicRenderer` to be used.
pub fn plot_geometry_with_renderer<R: SchematicRenderer>(
    system: &ChannelSystem,
    output_path: &str,
    renderer: &R,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    renderer.render_system(system, output_path, config)
}

fn explicit_or_generated_path(
    channel_spec: &crate::domain::model::ChannelSpec,
    from: &Node,
    to: &Node,
) -> Vec<Point2D> {
    if let Some(path) = channel_spec
        .path
        .as_ref()
        .map(|path| path.polyline_mm.clone())
        .or_else(|| {
            channel_spec
                .metadata
                .as_ref()
                .and_then(|meta| meta.get::<ChannelPathMetadata>())
                .map(|path| path.polyline_mm.clone())
        })
    {
        return path;
    }

    match channel_spec.channel_shape {
        ChannelShape::Serpentine {
            segments,
            bend_radius_m,
        } => generated_serpentine_path(from.point, to.point, segments, bend_radius_m * 1e3),
        ChannelShape::Straight => vec![from.point, to.point],
    }
}

fn channel_type_from_blueprint(
    channel_spec: &crate::domain::model::ChannelSpec,
    path: &[Point2D],
) -> ChannelType {
    if let Some(venturi) = channel_spec
        .venturi_geometry
        .as_ref()
        .or_else(|| {
            channel_spec
                .metadata
                .as_ref()
                .and_then(|meta| meta.get::<VenturiGeometryMetadata>())
        })
    {
        let inlet_width = venturi.inlet_width_m * 1e3;
        let throat_width = venturi.throat_width_m * 1e3;
        let outlet_width = venturi.outlet_width_m * 1e3;
        return ChannelType::Frustum {
            path: path.to_vec(),
            widths: frustum_width_profile(path.len(), inlet_width, throat_width, outlet_width),
            inlet_width,
            throat_width,
            outlet_width,
            taper_profile: crate::config::TaperProfile::Smooth,
            throat_position: 0.5,
            has_venturi_throat: true,
        };
    }

    match channel_spec.channel_shape {
        ChannelShape::Serpentine { .. } => ChannelType::Serpentine {
            path: path.to_vec(),
        },
        ChannelShape::Straight => {
            if path.len() > 2 {
                ChannelType::SmoothStraight {
                    path: path.to_vec(),
                }
            } else {
                ChannelType::Straight
            }
        }
    }
}

fn frustum_width_profile(
    point_count: usize,
    inlet_width: f64,
    throat_width: f64,
    outlet_width: f64,
) -> Vec<f64> {
    let n = point_count.max(2);
    let mut widths = Vec::with_capacity(n);
    for idx in 0..n {
        let t = idx as f64 / (n - 1) as f64;
        let width = if t <= 0.5 {
            let local = t / 0.5;
            inlet_width + (throat_width - inlet_width) * local
        } else {
            let local = (t - 0.5) / 0.5;
            throat_width + (outlet_width - throat_width) * local
        };
        widths.push(width.max(1e-6));
    }
    widths
}

fn generated_serpentine_path(
    start: Point2D,
    end: Point2D,
    segments: usize,
    bend_radius_mm: f64,
) -> Vec<Point2D> {
    if segments < 2 {
        return vec![start, end];
    }

    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.abs().max(1e-6);
    let lane_mid = (start.1 + end.1) * 0.5;
    let amplitude = (bend_radius_mm.max(0.8) * 0.9).min((length / 10.0).max(1.2));
    let mut path = Vec::with_capacity(segments + 1);
    path.push(start);
    for step in 1..segments {
        let t = step as f64 / segments as f64;
        let x = start.0 + dx * t;
        let offset = if step % 2 == 0 { -amplitude } else { amplitude };
        let y = lane_mid + dy * t + offset;
        path.push((x, y));
    }
    path.push(end);
    path
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

/// Collect all centerline vertices from a rendered system, including routed
/// polyline bend points. Useful for report annotations that should expose the
/// actual routed path, not just graph-node endpoints.
#[must_use]
pub fn centerline_vertices(system: &ChannelSystem) -> Vec<Point2D> {
    let mut points = Vec::new();
    for channel in &system.channels {
        points.extend(centerline_for_channel(channel, &system.nodes));
    }
    points
}
