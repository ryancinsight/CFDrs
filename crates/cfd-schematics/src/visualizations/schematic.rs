use crate::domain::model::NetworkBlueprint;
use crate::error::VisualizationResult;
use crate::geometry::metadata::{BlueprintRenderHints, ChannelPathMetadata, NodeLayoutMetadata};
use crate::geometry::{ChannelTypeCategory, Point2D};
use crate::visualizations::annotations::SchematicAnnotations;
use crate::visualizations::plotters_backend::PlottersRenderer;
use crate::visualizations::traits::{RenderConfig, SchematicRenderer};
use serde::Serialize;
use std::collections::{HashMap, VecDeque};

pub(crate) struct RenderChannelSystem {
    pub box_outline: Vec<(Point2D, Point2D)>,
    pub channel_paths: Vec<Vec<Point2D>>,
    pub channel_categories: Vec<ChannelTypeCategory>,
}

pub fn plot_geometry(blueprint: &NetworkBlueprint, output_path: &str) -> VisualizationResult<()> {
    plot_blueprint_auto_annotated(blueprint, output_path, &RenderConfig::default())
}

pub fn plot_geometry_with_config(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    renderer.render_system(blueprint, output_path, config)
}

pub fn plot_geometry_with_annotations(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    let renderer = PlottersRenderer;
    let mut annotated_config = config.clone();
    annotated_config.annotations = Some(annotations.clone());
    renderer.render_system(blueprint, output_path, &annotated_config)
}

pub fn plot_blueprint(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    plot_geometry_with_config(blueprint, output_path, config)
}

pub fn plot_blueprint_with_annotations(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
    annotations: &SchematicAnnotations,
) -> VisualizationResult<()> {
    plot_geometry_with_annotations(blueprint, output_path, config, annotations)
}

pub fn plot_blueprint_auto_annotated(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    let annotations = build_auto_annotations(blueprint);
    plot_blueprint_with_annotations(blueprint, output_path, config, &annotations)
}

pub(crate) fn materialize_blueprint_layout(blueprint: &mut NetworkBlueprint) {
    let box_dims = blueprint.box_dims;
    let (node_points, auto_layout_used) = blueprint_node_positions(blueprint, box_dims);
    let channel_paths = resolved_channel_paths(blueprint, &node_points, box_dims);

    for node in &mut blueprint.nodes {
        let Some(&(x_mm, y_mm)) = node_points.get(node.id.as_str()) else {
            continue;
        };
        if auto_layout_used.contains_key(node.id.as_str()) {
            let layout = NodeLayoutMetadata { x_mm, y_mm };
            node.point = (x_mm, y_mm);
            node.layout = Some(layout);
            node.metadata
                .get_or_insert_with(crate::geometry::metadata::MetadataContainer::new)
                .insert(layout);
        }
    }

    for (channel, path) in blueprint.channels.iter_mut().zip(channel_paths) {
        if channel.path.is_empty() && path.len() >= 2 {
            channel.path = path;
        }
    }
}

pub(crate) fn channel_system_from_blueprint(
    blueprint: &NetworkBlueprint,
    box_dims_hint: Option<(f64, f64)>,
    output_path: Option<&str>,
) -> crate::visualizations::schematic::RenderChannelSystem {
    let box_dims = box_dims_hint.unwrap_or(blueprint.box_dims);
    let (node_points, auto_layout_used) = blueprint_node_positions(blueprint, box_dims);
    if !auto_layout_used.is_empty() {
        if let Some(path) = output_path {
            save_auto_layout_json(&auto_layout_used, box_dims, path);
        }
    }
    let mut channel_paths = Vec::with_capacity(blueprint.channels.len());
    let mut channel_categories = Vec::with_capacity(blueprint.channels.len());
    let resolved_paths = resolved_channel_paths(blueprint, &node_points, box_dims);

    for (channel_spec, raw_path) in blueprint.channels.iter().zip(resolved_paths) {
        let category = channel_category_from_blueprint(channel_spec);
        channel_paths.push(render_path_for_display(&raw_path, category, box_dims));
        channel_categories.push(category);
    }

    crate::visualizations::schematic::RenderChannelSystem {
        box_outline: blueprint_box_outline(blueprint, box_dims),
        channel_paths,
        channel_categories,
    }
}

/// Returns `(all_positions, auto_layout_applied)` where `auto_layout_applied` contains
/// only the entries whose coordinates were determined by the topological auto-layout
/// fallback (i.e. the node had no explicit position).  The second map can be serialised
/// to JSON to record — and later replay — the computed layout.
fn blueprint_node_positions(
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

fn resolved_channel_paths(
    blueprint: &NetworkBlueprint,
    node_points: &HashMap<String, Point2D>,
    box_dims: (f64, f64),
) -> Vec<Vec<Point2D>> {
    let parallel_groups = parallel_channel_groups(blueprint);
    let default_point = (box_dims.0 * 0.5, box_dims.1 * 0.5);

    blueprint
        .channels
        .iter()
        .enumerate()
        .map(|(channel_idx, channel_spec)| {
            let from = node_points
                .get(channel_spec.from.as_str())
                .copied()
                .unwrap_or(default_point);
            let to = node_points
                .get(channel_spec.to.as_str())
                .copied()
                .unwrap_or(default_point);
            let key = (
                channel_spec.from.as_str().to_string(),
                channel_spec.to.as_str().to_string(),
            );
            let parallel_slot = parallel_groups.get(&key).and_then(|indices| {
                indices
                    .iter()
                    .position(|idx| *idx == channel_idx)
                    .map(|slot| (slot, indices.len()))
            });
            explicit_or_generated_path(channel_spec, from, to, parallel_slot, box_dims)
        })
        .collect()
}

fn parallel_channel_groups(blueprint: &NetworkBlueprint) -> HashMap<(String, String), Vec<usize>> {
    let mut groups: HashMap<(String, String), Vec<usize>> = HashMap::new();
    for (channel_idx, channel) in blueprint.channels.iter().enumerate() {
        groups
            .entry((
                channel.from.as_str().to_string(),
                channel.to.as_str().to_string(),
            ))
            .or_default()
            .push(channel_idx);
    }
    groups
}

/// Serialisable record saved as `{stem}_layout.json` whenever the topological
/// auto-layout fallback determines node positions.  Loading this file and
/// injecting the coordinates back into `NodeSpec::layout` or
/// `NodeLayoutMetadata` restores the exact same positions without re-running
/// the auto-layout algorithm.
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

fn save_auto_layout_json(
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

fn blueprint_box_outline(
    blueprint: &NetworkBlueprint,
    box_dims: (f64, f64),
) -> Vec<(Point2D, Point2D)> {
    let mut outline = vec![
        ((0.0, 0.0), (box_dims.0, 0.0)),
        ((box_dims.0, 0.0), (box_dims.0, box_dims.1)),
        ((box_dims.0, box_dims.1), (0.0, box_dims.1)),
        ((0.0, box_dims.1), (0.0, 0.0)),
    ];

    for node_spec in &blueprint.nodes {
        match node_spec.kind {
            crate::domain::model::NodeKind::Inlet => {
                outline.push(((0.0, node_spec.point.1), node_spec.point));
            }
            crate::domain::model::NodeKind::Outlet => {
                outline.push((node_spec.point, (box_dims.0, node_spec.point.1)));
            }
            crate::domain::model::NodeKind::Reservoir
            | crate::domain::model::NodeKind::Junction => {}
        }
    }

    outline
}

pub fn plot_geometry_auto_annotated(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    plot_blueprint_auto_annotated(blueprint, output_path, config)
}

fn build_auto_annotations(blueprint: &NetworkBlueprint) -> SchematicAnnotations {
    use crate::visualizations::annotations::{classify_node_roles, AnnotationMarker, MarkerRole};

    let mut annotations = SchematicAnnotations::report_default();
    let mut roles = classify_node_roles(blueprint);

    let has_target =
        blueprint.length_in_zone(crate::domain::therapy_metadata::TherapyZone::CancerTarget) > 0.0;
    let has_bypass =
        blueprint.length_in_zone(crate::domain::therapy_metadata::TherapyZone::HealthyBypass) > 0.0;

    if has_target || has_bypass {
        let box_dims = blueprint.box_dims;
        let y_mid = box_dims.1 * 0.5;
        let tx_min = box_dims.0 * 0.35;
        let tx_max = box_dims.0 * 0.65;
        for (node_idx, node) in blueprint.nodes.iter().enumerate() {
            let base = roles
                .get(&node_idx)
                .copied()
                .unwrap_or(MarkerRole::Internal);
            if matches!(
                base,
                MarkerRole::Inlet | MarkerRole::Outlet | MarkerRole::Split | MarkerRole::Merge
            ) {
                continue;
            }
            let near_cl = (node.point.1 - y_mid).abs() <= box_dims.1 * 0.18;
            let in_win = node.point.0 >= tx_min && node.point.0 <= tx_max;
            if has_target && near_cl && in_win {
                roles.insert(node_idx, MarkerRole::TherapyTarget);
            } else if has_bypass && !near_cl {
                roles.insert(node_idx, MarkerRole::Bypass);
            }
        }
    }

    let mut split_idx = 1usize;
    let mut merge_idx = 1usize;
    let mut inlet_idx = 1usize;
    let mut outlet_idx = 1usize;
    for (node_idx, node) in blueprint.nodes.iter().enumerate() {
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

    let hints: Option<&BlueprintRenderHints> = blueprint.render_hints();
    let mut throat_count =
        crate::visualizations::annotations::throat_count_from_blueprint_metadata(blueprint);
    if throat_count == 0 {
        if let Some(h) = hints {
            throat_count = h.throat_count_hint;
        }
    }

    if throat_count > 0 {
        let actual_points =
            crate::visualizations::annotations::venturi_marker_points_from_blueprint(blueprint);
        let mut th_idx = 1usize;
        if actual_points.is_empty() {
            let main_path = crate::visualizations::annotations::center_biased_main_path(blueprint);
            let box_dims = blueprint.box_dims;
            let zone = (box_dims.0 * 0.35, box_dims.0 * 0.65);
            for point in crate::visualizations::annotations::project_markers_along_path(
                &main_path,
                throat_count,
                zone,
            ) {
                annotations
                    .markers
                    .push(AnnotationMarker::new(point, MarkerRole::VenturiThroat).with_label(format!("TH{th_idx}"), true));
                th_idx += 1;
            }
        } else {
            for point in actual_points {
                annotations
                    .markers
                    .push(AnnotationMarker::new(point, MarkerRole::VenturiThroat).with_label(format!("TH{th_idx}"), true));
                th_idx += 1;
            }
        }
    }

    if let Some(h) = hints {
        let min_w = blueprint
            .channels
            .iter()
            .map(|ch| ch.cross_section.dims().0)
            .fold(f64::INFINITY, f64::min);
        let max_w = blueprint
            .channels
            .iter()
            .map(|ch| ch.cross_section.dims().0)
            .fold(f64::NEG_INFINITY, f64::max);
        annotations.legend_note = Some(format!(
            "{}  |  seq: {}  |  layers: {}  |  throats: {}  |  width: {:.2}-{:.2} mm  |  line thickness ∝ channel width",
            h.treatment_label,
            h.stage_sequence,
            h.split_layers,
            throat_count,
            min_w * 1e3,
            max_w * 1e3
        ));
        annotations.markers.push(
            AnnotationMarker::new(
                (blueprint.box_dims.0 * 0.50, blueprint.box_dims.1 * 0.82),
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

    let volume_note = blueprint.fluid_volume_summary().display_label;
    let y_fraction = 0.88;
    annotations.legend_note = Some(match annotations.legend_note.take() {
        Some(note) => format!("{note} | {volume_note}"),
        None => volume_note.clone(),
    });
    annotations.markers.push(
        crate::visualizations::annotations::AnnotationMarker::new(
            (
                blueprint.box_dims.0 * 0.50,
                blueprint.box_dims.1 * y_fraction,
            ),
            crate::visualizations::annotations::MarkerRole::Internal,
        )
        .with_label(volume_note, true),
    );

    annotations
}

fn explicit_or_generated_path(
    channel_spec: &crate::domain::model::ChannelSpec,
    from: Point2D,
    to: Point2D,
    parallel_slot: Option<(usize, usize)>,
    box_dims: (f64, f64),
) -> Vec<Point2D> {
    if let crate::domain::model::ChannelShape::Serpentine {
        segments,
        bend_radius_m,
        ..
    } = channel_spec.channel_shape
    {
        if !channel_spec.path.is_empty() {
            return if path_has_visible_serpentine_curvature(
                &channel_spec.path,
                bend_radius_m * 1.0e3,
            ) {
                channel_spec.path.clone()
            } else {
                generated_serpentine_path(from, to, segments, bend_radius_m * 1.0e3)
            };
        }
    }

    if !channel_spec.path.is_empty() {
        return channel_spec.path.clone();
    }

    if let Some(path) = channel_spec
        .metadata
        .as_ref()
        .and_then(|meta| meta.get::<ChannelPathMetadata>())
        .map(|path| path.polyline_mm.clone())
    {
        return path;
    }

    if let Some((slot, count)) = parallel_slot.filter(|(_, count)| *count > 1) {
        return generated_parallel_path(from, to, slot, count, box_dims);
    }

    match channel_spec.channel_shape {
        crate::domain::model::ChannelShape::Serpentine {
            segments,
            bend_radius_m,
            ..
        } => generated_serpentine_path(from, to, segments, bend_radius_m * 1.0e3),
        crate::domain::model::ChannelShape::Straight => vec![from, to],
    }
}

fn generated_parallel_path(
    start: Point2D,
    end: Point2D,
    slot: usize,
    count: usize,
    box_dims: (f64, f64),
) -> Vec<Point2D> {
    if count <= 1 {
        return vec![start, end];
    }

    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    if length <= 1e-9 {
        return vec![start, end];
    }

    let nx = -dy / length;
    let ny = dx / length;
    let base_spacing_mm = (length * 0.18).clamp(6.0, box_dims.1 * 0.18);
    let centered_slot = slot as f64 - (count as f64 - 1.0) * 0.5;
    let spread = if count == 2 { 2.0 } else { 1.0 };
    let offset_mm = centered_slot * base_spacing_mm * spread;
    if offset_mm.abs() <= 1e-9 {
        return vec![start, end];
    }

    let shoulder = 0.28;
    let entry = (
        start.0 + dx * shoulder + nx * offset_mm,
        start.1 + dy * shoulder + ny * offset_mm,
    );
    let exit = (
        end.0 - dx * shoulder + nx * offset_mm,
        end.1 - dy * shoulder + ny * offset_mm,
    );
    vec![start, entry, exit, end]
}

fn path_has_visible_serpentine_curvature(path: &[Point2D], bend_radius_mm: f64) -> bool {
    if path.len() < 3 {
        return false;
    }

    let start = path[0];
    let end = *path.last().expect("path has at least one point");
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let chord = dx.hypot(dy);
    if chord <= 1.0e-9 {
        return false;
    }

    let nx = -dy / chord;
    let ny = dx / chord;
    let min_expected_offset_mm = (bend_radius_mm * 0.25).max(0.75);
    let offsets = path
        .iter()
        .map(|point| ((point.0 - start.0) * nx) + ((point.1 - start.1) * ny))
        .filter(|offset| offset.abs() >= min_expected_offset_mm)
        .collect::<Vec<_>>();

    offsets.iter().any(|offset| *offset > 0.0) && offsets.iter().any(|offset| *offset < 0.0)
}

fn channel_category_from_blueprint(
    channel_spec: &crate::domain::model::ChannelSpec,
) -> ChannelTypeCategory {
    if channel_spec.venturi_geometry.is_some() {
        return ChannelTypeCategory::Tapered;
    }
    match channel_spec.channel_shape {
        crate::domain::model::ChannelShape::Serpentine { .. } => ChannelTypeCategory::Curved,
        crate::domain::model::ChannelShape::Straight => ChannelTypeCategory::Straight,
    }
}

/// Generates a square-wave serpentine path between two endpoints.
///
/// Each segment produces one horizontal run offset alternately ±amplitude
/// from the centre-line (start→end), connected by perpendicular 90° steps.
/// This gives a clear square-wave visual distinct from a smooth sine curve.
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
    let length = (dx * dx + dy * dy).sqrt();
    if length < 1e-6 {
        return vec![start, end];
    }

    let px = dx / length;
    let py = dy / length;
    let nx = -py;
    let ny = px;
    let requested_turns = segments.saturating_sub(1).max(2);
    let step = length / requested_turns as f64;
    let shoulder = (step * 0.28).max(length * 0.04);
    let amplitude = (bend_radius_mm * 1.4)
        .max(length / ((requested_turns as f64 + 1.0) * 3.0))
        .min(length * 0.24);

    let mut path = Vec::with_capacity(2 + requested_turns * 5);
    path.push(start);
    for turn_idx in 0..requested_turns {
        let center = step * (turn_idx as f64 + 0.5);
        let sign = if turn_idx % 2 == 0 { 1.0 } else { -1.0 };
        for (along, scale) in [
            (center - shoulder, 0.35),
            (center - shoulder * 0.45, 0.72),
            (center, 1.0),
            (center + shoulder * 0.45, 0.72),
            (center + shoulder, 0.35),
        ] {
            if along <= 1.0e-6 || along >= length - 1.0e-6 {
                continue;
            }
            let offset = sign * amplitude * scale;
            path.push((
                start.0 + px * along + nx * offset,
                start.1 + py * along + ny * offset,
            ));
        }
    }
    path.push(end);
    path
}

fn render_path_for_display(
    path: &[Point2D],
    category: ChannelTypeCategory,
    box_dims: (f64, f64),
) -> Vec<Point2D> {
    if path.len() <= 2 {
        return path.to_vec();
    }

    let tolerance_mm = (box_dims.0.max(box_dims.1) / 1200.0).clamp(0.10, 0.40);
    let simplified = simplify_polyline(path, tolerance_mm);

    if !matches!(category, ChannelTypeCategory::Curved) || simplified.len() < 4 {
        return simplified;
    }

    let smoothed = catmull_rom_resample(&simplified, 5);
    simplify_polyline(&smoothed, tolerance_mm * 0.75)
}

fn simplify_polyline(path: &[Point2D], tolerance_mm: f64) -> Vec<Point2D> {
    if path.len() <= 2 {
        return path.to_vec();
    }

    let mut keep = vec![false; path.len()];
    keep[0] = true;
    keep[path.len() - 1] = true;
    simplify_polyline_range(path, 0, path.len() - 1, tolerance_mm, &mut keep);

    let simplified: Vec<Point2D> = path
        .iter()
        .zip(keep)
        .filter_map(|(point, keep_point)| keep_point.then_some(*point))
        .collect();

    if simplified.len() < 2 {
        vec![path[0], *path.last().expect("path has at least one point")]
    } else {
        simplified
    }
}

fn simplify_polyline_range(
    path: &[Point2D],
    start_idx: usize,
    end_idx: usize,
    tolerance_mm: f64,
    keep: &mut [bool],
) {
    if end_idx <= start_idx + 1 {
        return;
    }

    let start = path[start_idx];
    let end = path[end_idx];
    let mut max_distance = 0.0;
    let mut max_idx = None;

    for (offset, point) in path[start_idx + 1..end_idx].iter().enumerate() {
        let idx = start_idx + 1 + offset;
        let distance = perpendicular_distance_to_segment(*point, start, end);
        if distance > max_distance {
            max_distance = distance;
            max_idx = Some(idx);
        }
    }

    if max_distance > tolerance_mm {
        if let Some(idx) = max_idx {
            keep[idx] = true;
            simplify_polyline_range(path, start_idx, idx, tolerance_mm, keep);
            simplify_polyline_range(path, idx, end_idx, tolerance_mm, keep);
        }
    }
}

fn perpendicular_distance_to_segment(point: Point2D, seg_start: Point2D, seg_end: Point2D) -> f64 {
    let (x0, y0) = point;
    let (x1, y1) = seg_start;
    let (x2, y2) = seg_end;
    let dx = x2 - x1;
    let dy = y2 - y1;

    if dx.abs() < f64::EPSILON && dy.abs() < f64::EPSILON {
        return ((x0 - x1).powi(2) + (y0 - y1).powi(2)).sqrt();
    }

    let t = (((x0 - x1) * dx) + ((y0 - y1) * dy)) / (dx * dx + dy * dy);
    let t = t.clamp(0.0, 1.0);
    let proj_x = x1 + t * dx;
    let proj_y = y1 + t * dy;
    ((x0 - proj_x).powi(2) + (y0 - proj_y).powi(2)).sqrt()
}

fn catmull_rom_resample(path: &[Point2D], samples_per_segment: usize) -> Vec<Point2D> {
    if path.len() < 4 || samples_per_segment == 0 {
        return path.to_vec();
    }

    let mut out = Vec::with_capacity((path.len() - 1) * samples_per_segment + 1);
    out.push(path[0]);

    for idx in 0..path.len() - 1 {
        let p0 = if idx == 0 { path[idx] } else { path[idx - 1] };
        let p1 = path[idx];
        let p2 = path[idx + 1];
        let p3 = if idx + 2 < path.len() {
            path[idx + 2]
        } else {
            path[idx + 1]
        };

        for sample_idx in 1..=samples_per_segment {
            let t = sample_idx as f64 / samples_per_segment as f64;
            out.push(catmull_rom_point(p0, p1, p2, p3, t));
        }
    }

    out
}

fn catmull_rom_point(p0: Point2D, p1: Point2D, p2: Point2D, p3: Point2D, t: f64) -> Point2D {
    let t2 = t * t;
    let t3 = t2 * t;
    let x = 0.5
        * ((2.0 * p1.0)
            + (-p0.0 + p2.0) * t
            + (2.0 * p0.0 - 5.0 * p1.0 + 4.0 * p2.0 - p3.0) * t2
            + (-p0.0 + 3.0 * p1.0 - 3.0 * p2.0 + p3.0) * t3);
    let y = 0.5
        * ((2.0 * p1.1)
            + (-p0.1 + p2.1) * t
            + (2.0 * p0.1 - 5.0 * p1.1 + 4.0 * p2.1 - p3.1) * t2
            + (-p0.1 + 3.0 * p1.1 - 3.0 * p2.1 + p3.1) * t3);
    (x, y)
}

/// Plot a channel system using a custom renderer
///
/// This function demonstrates the flexibility of the abstracted system
/// by allowing any renderer that implements `SchematicRenderer` to be used.
pub fn plot_geometry_with_renderer<R: SchematicRenderer>(
    blueprint: &NetworkBlueprint,
    output_path: &str,
    renderer: &R,
    config: &RenderConfig,
) -> VisualizationResult<()> {
    renderer.render_system(blueprint, output_path, config)
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
        columns
            .entry(node_depth)
            .or_default()
            .push(node.id.as_str());
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
pub fn centerline_vertices(blueprint: &NetworkBlueprint) -> Vec<Point2D> {
    channel_system_from_blueprint(blueprint, Some(blueprint.box_dims), None)
        .channel_paths
        .into_iter()
        .flatten()
        .collect()
}

#[cfg(test)]
mod tests {
    use super::{channel_category_from_blueprint, generated_serpentine_path, render_path_for_display};
    use crate::domain::model::{ChannelShape, ChannelSpec};
    use crate::geometry::ChannelTypeCategory;

    #[test]
    fn generated_serpentine_path_has_mirrored_two_lobe_offsets() {
        let path = generated_serpentine_path((0.0, 0.0), (20.0, 0.0), 2, 2.0);
        let signed_offsets: Vec<f64> = path
            .iter()
            .skip(1)
            .take(path.len().saturating_sub(2))
            .map(|point| point.1)
            .filter(|offset| offset.abs() > 1.0e-6)
            .collect();
        assert!(
            signed_offsets.iter().any(|offset| *offset > 0.0),
            "serpentine path must lobe above the centerline"
        );
        assert!(
            signed_offsets.iter().any(|offset| *offset < 0.0),
            "serpentine path must lobe below the centerline"
        );
    }

    #[test]
    fn serpentine_channels_render_as_curved_paths() {
        let path = generated_serpentine_path((0.0, 0.0), (20.0, 0.0), 5, 2.5);
        let rendered = render_path_for_display(&path, ChannelTypeCategory::Curved, (40.0, 20.0));
        assert!(
            rendered.len() > path.len(),
            "curved serpentine display path should be spline-resampled"
        );
    }

    #[test]
    fn serpentine_channel_category_is_curved() {
        let mut channel = ChannelSpec::new_pipe_rect(
            "serp",
            "a".to_string(),
            "b".to_string(),
            1.0,
            1.0e-3,
            1.0e-3,
            0.0,
            0.0,
        );
        channel.channel_shape = ChannelShape::Serpentine {
            segments: 5,
            bend_radius_m: 2.0e-3,
            wave_type: crate::topology::SerpentineWaveType::Sine,
        };
        assert_eq!(
            channel_category_from_blueprint(&channel),
            ChannelTypeCategory::Curved
        );
    }
}
