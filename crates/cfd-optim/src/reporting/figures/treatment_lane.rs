use std::fmt::Write as _;
use std::path::Path;

use cfd_schematics::domain::model::{ChannelShape, ChannelSpec};
use cfd_schematics::geometry::metadata::ChannelVisualRole;

use super::primitives::{escape_xml, svg_end, svg_start, svg_title};
use super::process::write_placeholder;
use crate::reporting::Milestone12ReportDesign;

#[derive(Clone)]
struct LaneStroke {
    path: Vec<(f64, f64)>,
    color: &'static str,
    width_px: f64,
    venturi_label: Option<String>,
}

#[derive(Clone)]
struct LanePanelData {
    title: String,
    subtitle: String,
    strokes: Vec<LaneStroke>,
    serpentine_channel_count: usize,
    venturi_count: usize,
    treatment_residence_ms: f64,
    cavitation_intensity: f64,
    cavitation_number: f64,
    total_loss_coefficient: f64,
    total_path_length_mm: f64,
    min_bend_radius_mm: Option<f64>,
    max_bend_radius_mm: Option<f64>,
}

pub(super) fn write_treatment_lane_zoom_figure(
    path: &Path,
    option2: &Milestone12ReportDesign,
    ga_best: &Milestone12ReportDesign,
) -> Result<(), Box<dyn std::error::Error>> {
    let option2_panel = panel_data(option2, "Option 2");
    let ga_panel = panel_data(ga_best, "GA");

    if option2_panel.strokes.is_empty() || ga_panel.strokes.is_empty() {
        return write_placeholder(
            path,
            "Treatment-Lane Geometry Zoom - Option 2 vs GA",
            "Treatment-lane path geometry was unavailable for one or both selected designs.",
        );
    }

    let mut svg = String::new();
    let width = 1280.0;
    let height = 760.0;
    svg_start(&mut svg, width, height);
    svg_title(&mut svg, "Treatment-Lane Geometry Zoom - Option 2 vs GA");

    svg.push_str(r##"<rect x="0" y="0" width="1280" height="760" fill="#fbfcfd"/>"##);
    svg.push_str(r##"<text x="640" y="54" text-anchor="middle" font-size="16" fill="#566573">Normalized side-by-side view of the local treatment path only; full-device split-tree context removed</text>"##);

    render_panel(&mut svg, &option2_panel, 60.0, 110.0, 540.0, 560.0);
    render_panel(&mut svg, &ga_panel, 680.0, 110.0, 540.0, 560.0);
    render_delta_overlay(&mut svg, &option2_panel, &ga_panel, 470.0, 640.0, 340.0, 74.0);

    svg.push_str(r##"<line x1="640" y1="120" x2="640" y2="670" stroke="#d5dbdb" stroke-width="2" stroke-dasharray="8 6"/>"##);
    svg.push_str(r##"<line x1="110" y1="706" x2="136" y2="706" stroke="#8C32A0" stroke-width="7" stroke-linecap="round"/>"##);
    svg.push_str(r##"<text x="144" y="711" font-size="13" fill="#34495e">Treatment-channel path</text>"##);
    svg.push_str(r##"<line x1="320" y1="706" x2="346" y2="706" stroke="#E67E22" stroke-width="9" stroke-linecap="round"/>"##);
    svg.push_str(r##"<text x="354" y="711" font-size="13" fill="#34495e">Active venturi throat segment</text>"##);

    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

fn panel_data(design: &Milestone12ReportDesign, panel_label: &str) -> LanePanelData {
    let blueprint = design.candidate.blueprint();
    let mut strokes = Vec::new();
    let mut venturi_index = 1_usize;
    let mut serpentine_channel_count = 0_usize;
    let mut venturi_count = 0_usize;
    let mut total_path_length_mm = 0.0_f64;
    let mut min_bend_radius_mm = f64::INFINITY;
    let mut max_bend_radius_mm = f64::NEG_INFINITY;

    for channel in blueprint.channels.iter().filter(is_treatment_lane_channel) {
        if channel.path.len() < 2 {
            continue;
        }
        let display_path = display_path(channel);
        total_path_length_mm += polyline_length_mm(&display_path);
        if matches!(channel.channel_shape, ChannelShape::Serpentine { .. }) {
            serpentine_channel_count += 1;
        }
        if let Some((local_min_radius_mm, local_max_radius_mm)) = path_bend_radius_range_mm(&display_path) {
            if !matches!(channel.channel_shape, ChannelShape::Serpentine { .. }) {
                serpentine_channel_count += 1;
            }
            min_bend_radius_mm = min_bend_radius_mm.min(local_min_radius_mm);
            max_bend_radius_mm = max_bend_radius_mm.max(local_max_radius_mm);
        } else if let ChannelShape::Serpentine { bend_radius_m, .. } = channel.channel_shape {
            let radius_mm = bend_radius_m * 1.0e3;
            min_bend_radius_mm = min_bend_radius_mm.min(radius_mm);
            max_bend_radius_mm = max_bend_radius_mm.max(radius_mm);
        }
        let is_venturi = channel.visual_role == Some(ChannelVisualRole::VenturiThroat)
            || channel.venturi_geometry.is_some();
        if is_venturi {
            venturi_count += 1;
        }

        let width_mm = channel.cross_section.dims().0 * 1.0e3;
        let width_px = if is_venturi {
            8.0
        } else {
            (2.5 + width_mm / 1.5).clamp(3.0, 7.0)
        };

        strokes.push(LaneStroke {
            path: display_path,
            color: if is_venturi { "#E67E22" } else { "#8C32A0" },
            width_px,
            venturi_label: is_venturi.then(|| {
                let label = format!("TH{venturi_index}");
                venturi_index += 1;
                label
            }),
        });
    }

    LanePanelData {
        title: format!("{panel_label} - {}", design.stage_sequence_label()),
        subtitle: format!(
            "score {:.4} • residence {:.1} ms • cavitation {:.3} • σ {:.3} • K_loss {:.2}",
            design.score,
            design.metrics.mean_residence_time_s.max(0.0) * 1.0e3,
            design.metrics.cavitation_intensity.max(0.0),
            design.metrics.cavitation_number,
            design.metrics.venturi_total_loss_coefficient,
        ),
        strokes,
        serpentine_channel_count,
        venturi_count,
        treatment_residence_ms: design.metrics.mean_residence_time_s.max(0.0) * 1.0e3,
        cavitation_intensity: design.metrics.cavitation_intensity.max(0.0),
        cavitation_number: design.metrics.cavitation_number,
        total_loss_coefficient: design.metrics.venturi_total_loss_coefficient,
        total_path_length_mm,
        min_bend_radius_mm: min_bend_radius_mm.is_finite().then_some(min_bend_radius_mm),
        max_bend_radius_mm: max_bend_radius_mm.is_finite().then_some(max_bend_radius_mm),
    }
}

fn is_treatment_lane_channel(channel: &&ChannelSpec) -> bool {
    matches!(
        channel.visual_role,
        Some(ChannelVisualRole::CenterTreatment | ChannelVisualRole::VenturiThroat)
    )
}

fn display_path(channel: &ChannelSpec) -> Vec<(f64, f64)> {
    match channel.channel_shape {
        ChannelShape::Serpentine {
            segments,
            bend_radius_m,
            ..
        } if channel.path.len() >= 2 => {
            if path_has_visible_serpentine_curvature(&channel.path, bend_radius_m * 1.0e3) {
                channel.path.clone()
            } else {
                generated_serpentine_path(
                    channel.path[0],
                    *channel.path.last().expect("channel.path has at least one point"),
                    segments,
                    bend_radius_m * 1.0e3,
                )
            }
        }
        _ => channel.path.clone(),
    }
}

fn render_panel(
    svg: &mut String,
    panel: &LanePanelData,
    x: f64,
    y: f64,
    width: f64,
    height: f64,
) {
    let _ = write!(
        svg,
        r##"<rect x="{x:.1}" y="{y:.1}" width="{width:.1}" height="{height:.1}" rx="18" fill="#ffffff" stroke="#d5dbdb" stroke-width="2"/>"##
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="20" font-weight="700" fill="#1f2d3d">{}</text>"##,
        x + 22.0,
        y + 32.0,
        escape_xml(&panel.title)
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="13" fill="#5d6d7e">{}</text>"##,
        x + 22.0,
        y + 56.0,
        escape_xml(&panel.subtitle)
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="12" fill="#7b8a8b">Local treatment lane normalized to panel frame</text>"##,
        x + 22.0,
        y + 76.0
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="12" fill="#7b8a8b">Path {:.1} mm • Rmin {} • Rmax {} • throat count {}</text>"##,
        x + 22.0,
        y + 94.0,
        panel.total_path_length_mm,
        format_radius_mm(panel.min_bend_radius_mm),
        format_radius_mm(panel.max_bend_radius_mm),
        panel.venturi_count
    );

    if panel.strokes.is_empty() {
        return;
    }

    let (min_x, max_x, min_y, max_y) = bounds(&panel.strokes);
    let inner_x = x + 24.0;
    let inner_y = y + 108.0;
    let inner_w = width - 48.0;
    let inner_h = height - 162.0;
    let span_x = (max_x - min_x).max(1.0);
    let span_y = (max_y - min_y).max(1.0);
    let scale = (inner_w / span_x).min(inner_h / span_y);
    let x_offset = inner_x + (inner_w - span_x * scale) * 0.5;
    let y_offset = inner_y + (inner_h - span_y * scale) * 0.5;

    let _ = write!(
        svg,
        r##"<rect x="{inner_x:.1}" y="{inner_y:.1}" width="{inner_w:.1}" height="{inner_h:.1}" rx="12" fill="#fcfcfd" stroke="#eef1f3" stroke-width="1.5"/>"##
    );

    for stroke in &panel.strokes {
        let points = stroke
            .path
            .iter()
            .map(|(px, py)| {
                let sx = x_offset + (px - min_x) * scale;
                let sy = y_offset + (max_y - py) * scale;
                format!("{sx:.1},{sy:.1}")
            })
            .collect::<Vec<_>>()
            .join(" ");
        let _ = write!(
            svg,
            r#"<polyline points="{points}" fill="none" stroke="{}" stroke-width="{:.1}" stroke-linecap="round" stroke-linejoin="round"/>"#,
            stroke.color,
            stroke.width_px
        );

        if let Some(label) = &stroke.venturi_label {
            let (cx, cy) = midpoint(&stroke.path);
            let sx = x_offset + (cx - min_x) * scale;
            let sy = y_offset + (max_y - cy) * scale;
            let _ = write!(
                svg,
                r##"<circle cx="{sx:.1}" cy="{sy:.1}" r="8" fill="#fff7ed" stroke="#d35400" stroke-width="2"/>"##
            );
            let _ = write!(
                svg,
                r##"<text x="{:.1}" y="{:.1}" font-size="11" text-anchor="middle" fill="#a04000" font-weight="700">{}</text>"##,
                sx,
                sy + 4.0,
                escape_xml(label)
            );
        }
    }

    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="12" fill="#7b8a8b">Curvature-bearing treatment segments: {}</text>"##,
        x + 22.0,
        y + height - 22.0,
        panel.serpentine_channel_count
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="12" fill="#7b8a8b">Venturi throat segments: {}</text>"##,
        x + width - 190.0,
        y + height - 22.0,
        panel.venturi_count
    );
}

fn bounds(strokes: &[LaneStroke]) -> (f64, f64, f64, f64) {
    let mut min_x = f64::INFINITY;
    let mut max_x = f64::NEG_INFINITY;
    let mut min_y = f64::INFINITY;
    let mut max_y = f64::NEG_INFINITY;

    for stroke in strokes {
        for (x, y) in &stroke.path {
            min_x = min_x.min(*x);
            max_x = max_x.max(*x);
            min_y = min_y.min(*y);
            max_y = max_y.max(*y);
        }
    }

    (min_x, max_x, min_y, max_y)
}

fn midpoint(path: &[(f64, f64)]) -> (f64, f64) {
    path.get(path.len() / 2).copied().unwrap_or((0.0, 0.0))
}

fn render_delta_overlay(
    svg: &mut String,
    option2: &LanePanelData,
    ga: &LanePanelData,
    x: f64,
    y: f64,
    width: f64,
    height: f64,
) {
    let residence_delta_ms = ga.treatment_residence_ms - option2.treatment_residence_ms;
    let cavitation_delta = ga.cavitation_intensity - option2.cavitation_intensity;
    let sigma_delta = ga.cavitation_number - option2.cavitation_number;
    let loss_delta = ga.total_loss_coefficient - option2.total_loss_coefficient;
    let bend_delta = match (option2.min_bend_radius_mm, ga.min_bend_radius_mm) {
        (Some(left), Some(right)) => format!("{:+.2} mm", right - left),
        _ => "n/a".to_string(),
    };
    let throat_delta = ga.venturi_count as isize - option2.venturi_count as isize;
    let serp_delta = ga.serpentine_channel_count as isize - option2.serpentine_channel_count as isize;

    let _ = write!(
        svg,
        r##"<rect x="{x:.1}" y="{y:.1}" width="{width:.1}" height="{height:.1}" rx="14" fill="#f6f8fa" stroke="#cfd8dc" stroke-width="1.5"/>"##
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="14" font-weight="700" fill="#1f2d3d">Geometry Delta (GA - Option 2)</text>"##,
        x + 18.0,
        y + 22.0
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="12" fill="#5d6d7e">Residence {:+.1} ms • cavitation {:+.3} • σ {:+.3} • K_loss {:+.2} • Throats {:+}</text>"##,
        x + 18.0,
        y + 44.0,
        residence_delta_ms,
        cavitation_delta,
        sigma_delta,
        loss_delta,
        throat_delta,
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="12" fill="#7b8a8b">Rmin {} • serpentine channels {:+}. TH labels are ordered inlet to outlet inside each normalized treatment-path frame.</text>"##,
        x + 18.0,
        y + 63.0,
        bend_delta,
        serp_delta
    );
}

fn polyline_length_mm(path: &[(f64, f64)]) -> f64 {
    path.windows(2)
        .map(|segment| {
            let dx = segment[1].0 - segment[0].0;
            let dy = segment[1].1 - segment[0].1;
            dx.hypot(dy)
        })
        .sum()
}

fn generated_serpentine_path(
    start: (f64, f64),
    end: (f64, f64),
    segments: usize,
    bend_radius_mm: f64,
) -> Vec<(f64, f64)> {
    if segments < 2 {
        return vec![start, end];
    }

    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    if length < 1.0e-6 {
        return vec![start, end];
    }

    let px = dx / length;
    let py = dy / length;
    let nx = -py;
    let ny = px;
    let turns = segments.saturating_sub(1).max(2);
    let step = length / turns as f64;
    let shoulder = (step * 0.28).max(length * 0.04);
    let amplitude = (bend_radius_mm * 1.4)
        .max(length / ((turns as f64 + 1.0) * 3.0))
        .min(length * 0.24);

    let mut path = Vec::with_capacity(2 + turns * 5);
    path.push(start);
    for turn_idx in 0..turns {
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

fn path_has_visible_serpentine_curvature(path: &[(f64, f64)], bend_radius_mm: f64) -> bool {
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

fn path_bend_radius_range_mm(path: &[(f64, f64)]) -> Option<(f64, f64)> {
    let radii = path
        .iter()
        .enumerate()
        .filter_map(|(idx, _)| estimate_local_bend_radius_mm(path, idx))
        .filter(|radius| radius.is_finite() && *radius < 1.0e6)
        .collect::<Vec<_>>();
    if radii.is_empty() {
        None
    } else {
        let min_radius = radii.iter().copied().fold(f64::INFINITY, f64::min);
        let max_radius = radii.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        Some((min_radius, max_radius))
    }
}

fn estimate_local_bend_radius_mm(path: &[(f64, f64)], idx: usize) -> Option<f64> {
    if idx == 0 || idx + 1 >= path.len() {
        return None;
    }

    let a = path[idx - 1];
    let b = path[idx];
    let c = path[idx + 1];
    let ab = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
    let bc = ((c.0 - b.0).powi(2) + (c.1 - b.1).powi(2)).sqrt();
    let ac = ((c.0 - a.0).powi(2) + (c.1 - a.1).powi(2)).sqrt();
    if ab <= 1.0e-9 || bc <= 1.0e-9 || ac <= 1.0e-9 {
        return None;
    }

    let twice_area = ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)).abs();
    if twice_area <= 1.0e-9 {
        return None;
    }

    Some((ab * bc * ac) / (2.0 * twice_area))
}

fn format_radius_mm(radius_mm: Option<f64>) -> String {
    radius_mm.map_or_else(|| "n/a".to_string(), |radius| format!("{radius:.2} mm"))
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use cfd_schematics::topology::VenturiPlacementMode;

    use crate::domain::fixtures::{operating_point, stage0_venturi_candidate};
    use crate::reporting::{compute_blueprint_report_metrics, Milestone12ReportDesign};

    use super::write_treatment_lane_zoom_figure;

    #[test]
    fn treatment_lane_zoom_renders_total_loss_annotation() {
        let option2_candidate = stage0_venturi_candidate(
            "lane-option2",
            operating_point(2.0e-6, 30_000.0, 0.18),
            VenturiPlacementMode::StraightSegment,
        );
        let ga_candidate = stage0_venturi_candidate(
            "lane-ga",
            operating_point(2.0e-6, 30_000.0, 0.18),
            VenturiPlacementMode::CurvaturePeakDeanNumber,
        );
        let option2 = Milestone12ReportDesign::new(
            1,
            option2_candidate.clone(),
            compute_blueprint_report_metrics(&option2_candidate).expect("option2 metrics"),
            0.70,
        );
        let ga = Milestone12ReportDesign::new(
            1,
            ga_candidate.clone(),
            compute_blueprint_report_metrics(&ga_candidate).expect("ga metrics"),
            0.75,
        );

        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("m12-treatment-lane-{unique}.svg"));
        write_treatment_lane_zoom_figure(&path, &option2, &ga)
            .expect("treatment lane figure should render");

        let rendered = std::fs::read_to_string(path).expect("rendered svg should exist");
        assert!(rendered.contains("K_loss"));
        assert!(rendered.contains("Geometry Delta (GA - Option 2)"));
    }
}