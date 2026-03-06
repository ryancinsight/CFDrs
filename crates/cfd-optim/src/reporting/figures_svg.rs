//! SVG figure builders for Milestone 12 narrative figures.

use std::fmt::Write as _;
use std::path::Path;

use crate::reporting::ValidationRow;
use crate::RankedDesign;

pub(super) fn write_placeholder(
    path: &Path,
    title: &str,
    message: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut svg = String::new();
    svg_start(&mut svg, 1100.0, 420.0);
    svg_title(&mut svg, title);
    let _ = write!(
        svg,
        r##"<rect x="70" y="120" width="960" height="220" fill="#f5f6fa" stroke="#bdc3c7" stroke-width="2"/>"##
    );
    let _ = write!(
        svg,
        r##"<text x="110" y="240" font-size="24" fill="#2c3e50">{}</text>"##,
        escape_xml(message)
    );
    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

pub(super) fn write_creation_optimization_process_figure(
    path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let mut svg = String::new();
    let w = 1400.0;
    let h = 760.0;
    svg_start(&mut svg, w, h);
    svg_title(
        &mut svg,
        "Milestone 12 Design Creation & Optimization Process",
    );
    let _ = write!(
        svg,
        r##"<text x="70" y="95" font-size="18" fill="#34495e">Primitive selective-routing pipeline for bifurcation/trifurcation treatment trees.</text>"##
    );

    for (x1, y1, x2, y2) in [
        (275.0, 205.0, 365.0, 205.0),
        (570.0, 205.0, 660.0, 205.0),
        (865.0, 205.0, 955.0, 205.0),
        (1160.0, 205.0, 1250.0, 205.0),
        (775.0, 280.0, 775.0, 390.0),
        (680.0, 485.0, 570.0, 485.0),
        (385.0, 485.0, 275.0, 485.0),
    ] {
        arrow(&mut svg, x1, y1, x2, y2);
    }

    process_box(
        &mut svg,
        70.0,
        130.0,
        205.0,
        150.0,
        "#eef6ff",
        "#2e86de",
        "1. Primitive candidate generation",
        &[
            "cfd-schematics::create_geometry",
            "Mirrored Bi/Tri split trees",
            "Branch-diameter biasing for RBC/WBC/CTC routing",
        ],
    );
    process_box(
        &mut svg,
        365.0,
        130.0,
        205.0,
        150.0,
        "#eefcf3",
        "#27ae60",
        "2. ChannelSystem topology SSOT",
        &[
            "Wall ports and true split/merge nodes",
            "Treatment/bypass lanes, serpentines, venturi throats",
            "No report-local topology recreation",
        ],
    );
    process_box(
        &mut svg,
        660.0,
        130.0,
        205.0,
        150.0,
        "#fff7ea",
        "#d68910",
        "3. NetworkBlueprint projection",
        &[
            "Lossless ChannelSystem -> NetworkBlueprint",
            "Downstream solver/export representation only",
            "No independent blueprint layout authoring",
        ],
    );
    process_box(
        &mut svg,
        955.0,
        130.0,
        205.0,
        150.0,
        "#fef0f3",
        "#c0392b",
        "4. Cached cfd-1d evaluation",
        &[
            "MetricsCache avoids repeat solves",
            "Pressure, shear, residence, cavitation",
            "Three-pop separation and safety metrics",
        ],
    );
    process_box(
        &mut svg,
        1250.0,
        130.0,
        95.0,
        150.0,
        "#f4f1ff",
        "#7d3c98",
        "5. Ranked outputs",
        &["Option 1", "Option 2", "GA compare"],
    );
    process_box(
        &mut svg,
        680.0,
        390.0,
        190.0,
        150.0,
        "#fff6fb",
        "#b03a78",
        "Option 2 combined ranking",
        &[
            "CombinedSdtLeukapheresis",
            "Selective routing + venturi cavitation",
            "Leukapheresis and safety gates",
        ],
    );
    process_box(
        &mut svg,
        385.0,
        390.0,
        185.0,
        150.0,
        "#eefcf3",
        "#1e8449",
        "Option 1 acoustic ranking",
        &[
            "SelectiveAcousticTherapy",
            "Center-lane ultrasound/light exposure",
            "RBC-biased bypass outside treatment zone",
        ],
    );
    process_box(
        &mut svg,
        90.0,
        390.0,
        185.0,
        150.0,
        "#f4f6f7",
        "#566573",
        "Report and artifact export",
        &[
            "Schematics from cfd-schematics only",
            "Canonical markdown, narrative, figure manifest",
            "Milestone 12 report figures and JSON",
        ],
    );

    let _ = write!(
        svg,
        r##"<text x="1110" y="415" font-size="18" fill="#34495e" font-weight="600">Option 1</text>
<text x="1110" y="440" font-size="14" fill="#34495e">Selective acoustic shortlist</text>
<text x="1110" y="470" font-size="18" fill="#34495e" font-weight="600">Option 2</text>
<text x="1110" y="495" font-size="14" fill="#34495e">Combined selective venturi shortlist</text>
<text x="1110" y="525" font-size="18" fill="#34495e" font-weight="600">Figure 6</text>
<text x="1110" y="550" font-size="14" fill="#34495e">HydroSDT GA comparison branch</text>"##
    );
    let _ = write!(
        svg,
        r##"<text x="70" y="625" font-size="14" fill="#566573">Selective routing remains the modeled concept: branch-diameter biasing pushes RBCs peripheral, keeps WBC/CTC enrichment in treatment branches, and changes only the treatment mechanism between Option 1 and Option 2.</text>"##
    );
    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

pub(super) fn write_cross_mode_figure(
    path: &Path,
    option1: &[RankedDesign],
    option2: &[RankedDesign],
    ga: &[RankedDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let data = vec![
        ("Opt1 Acoustic", option1.first().map_or(0.0, |d| d.score)),
        ("Opt2 Combined", option2.first().map_or(0.0, |d| d.score)),
        ("GA HydroSDT", ga.first().map_or(0.0, |d| d.score)),
    ];
    write_bar_svg(path, "Cross-Mode Scoring Comparison", &data, 1.0)
}

pub(super) fn write_head_to_head_figure(
    path: &Path,
    option2: &[RankedDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let data: Vec<(String, f64)> = option2
        .iter()
        .take(5)
        .map(|d| (format!("R{}", d.rank), d.score))
        .collect();
    if data.is_empty() {
        return write_placeholder(
            path,
            "Head-to-Head Design Comparison",
            "No Option 2 ranked data available.",
        );
    }
    let data_ref: Vec<(&str, f64)> = data.iter().map(|(k, v)| (k.as_str(), *v)).collect();
    write_bar_svg(path, "Head-to-Head Option 2 Top-5 Scores", &data_ref, 1.0)
}

pub(super) fn write_cavitation_distribution_figure(
    path: &Path,
    option2: &[RankedDesign],
    ga: &[RankedDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut neg = 0usize;
    let mut zero_to_one = 0usize;
    let mut ge_one = 0usize;
    let mut nonfinite = 0usize;

    for d in option2.iter().chain(ga.iter()) {
        let s = d.metrics.cavitation_number;
        if !s.is_finite() {
            nonfinite += 1;
        } else if s < 0.0 {
            neg += 1;
        } else if s < 1.0 {
            zero_to_one += 1;
        } else {
            ge_one += 1;
        }
    }

    let data = vec![
        ("σ<0", neg as f64),
        ("0≤σ<1", zero_to_one as f64),
        ("σ≥1", ge_one as f64),
        ("nonfinite", nonfinite as f64),
    ];
    let max = data.iter().map(|(_, v)| *v).fold(1.0, f64::max);
    write_bar_svg(
        path,
        "Cavitation Number (σ) Distribution",
        &data,
        (max + 1.0).max(1.0),
    )
}

pub(super) fn write_pareto_figure(
    path: &Path,
    option2: &[RankedDesign],
    ga: &[RankedDesign],
    option2_pool_all: &[RankedDesign],
    ga_pool_all: &[RankedDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let selected: Vec<(&RankedDesign, &str)> = option2
        .iter()
        .map(|d| (d, "Option2"))
        .chain(ga.iter().map(|d| (d, "GA")))
        .collect();
    let background: Vec<(&RankedDesign, &str)> = option2_pool_all
        .iter()
        .map(|d| (d, "Pool-O2"))
        .chain(ga_pool_all.iter().map(|d| (d, "Pool-GA")))
        .collect();
    if selected.is_empty() && background.is_empty() {
        return write_placeholder(
            path,
            "Pareto Front — Oncology Objectives",
            "No ranked data available.",
        );
    }

    // Compute data range for auto-scaling axes
    let all_pts: Vec<(f64, f64)> = background
        .iter()
        .chain(selected.iter())
        .map(|(d, _)| {
            (
                d.metrics.cancer_targeted_cavitation,
                d.metrics.therapeutic_window_score,
            )
        })
        .filter(|(x, y)| x.is_finite() && y.is_finite())
        .collect();
    let (x_min, x_max, y_min, y_max) = if all_pts.is_empty() {
        (0.0, 1.0, 0.0, 1.0)
    } else {
        let xmn = all_pts.iter().map(|p| p.0).fold(f64::INFINITY, f64::min);
        let xmx = all_pts.iter().map(|p| p.0).fold(f64::NEG_INFINITY, f64::max);
        let ymn = all_pts.iter().map(|p| p.1).fold(f64::INFINITY, f64::min);
        let ymx = all_pts.iter().map(|p| p.1).fold(f64::NEG_INFINITY, f64::max);
        let x_span = (xmx - xmn).max(0.01);
        let y_span = (ymx - ymn).max(0.01);
        (
            xmn - 0.1 * x_span,
            xmx + 0.1 * x_span,
            ymn - 0.1 * y_span,
            ymx + 0.1 * y_span,
        )
    };

    let mut svg = String::new();
    let w = 1100.0;
    let h = 700.0;
    svg_start(&mut svg, w, h);
    svg_title(&mut svg, "Pareto Front — Oncology Objectives");
    let x0 = 100.0;
    let y0 = 620.0;
    let xw = 900.0;
    let yh = 500.0;
    axis(&mut svg, x0, y0, xw, yh);

    // Tick labels — X axis
    for i in 0..=5 {
        let frac = i as f64 / 5.0;
        let val = x_min + frac * (x_max - x_min);
        let tx = x0 + xw * frac;
        let _ = write!(
            svg,
            r##"<text x="{tx:.1}" y="{:.1}" text-anchor="middle" font-size="12" fill="#7f8c8d">{val:.3}</text>"##,
            y0 + 18.0
        );
    }
    // Tick labels — Y axis
    for i in 0..=5 {
        let frac = i as f64 / 5.0;
        let val = y_min + frac * (y_max - y_min);
        let ty = y0 - yh * frac;
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{ty:.1}" text-anchor="end" font-size="12" fill="#7f8c8d">{val:.3}</text>"##,
            x0 - 8.0
        );
    }

    svg.push_str(r##"<text x="420" y="668" font-size="16" fill="#34495e">tumor_targeted_cavitation_index</text>"##);
    svg.push_str(r##"<text x="18" y="340" transform="rotate(-90 18,340)" font-size="16" fill="#34495e">therapeutic_window_score</text>"##);

    // Background pool: small translucent dots
    for (d, tag) in &background {
        let cx = d.metrics.cancer_targeted_cavitation;
        let cy = d.metrics.therapeutic_window_score;
        if !cx.is_finite() || !cy.is_finite() {
            continue;
        }
        let x = x0 + xw * ((cx - x_min) / (x_max - x_min)).clamp(0.0, 1.0);
        let y = y0 - yh * ((cy - y_min) / (y_max - y_min)).clamp(0.0, 1.0);
        let color = if *tag == "Pool-O2" {
            "#85c1e9"
        } else {
            "#f0b27a"
        };
        let _ = write!(
            svg,
            r#"<circle cx="{x:.2}" cy="{y:.2}" r="4" fill="{color}" fill-opacity="0.35"/>"#
        );
    }

    // Selected top designs: larger, opaque, with stroke
    for (d, tag) in &selected {
        let cx = d.metrics.cancer_targeted_cavitation;
        let cy = d.metrics.therapeutic_window_score;
        if !cx.is_finite() || !cy.is_finite() {
            continue;
        }
        let x = x0 + xw * ((cx - x_min) / (x_max - x_min)).clamp(0.0, 1.0);
        let y = y0 - yh * ((cy - y_min) / (y_max - y_min)).clamp(0.0, 1.0);
        let r = 6.0 + 6.0 * d.score.clamp(0.0, 1.0);
        let (fill, stroke) = if *tag == "Option2" {
            ("#2e86de", "#1a5276")
        } else {
            ("#d35400", "#873600")
        };
        let _ = write!(
            svg,
            r#"<circle cx="{x:.2}" cy="{y:.2}" r="{r:.2}" fill="{fill}" fill-opacity="0.85" stroke="{stroke}" stroke-width="2"/>"#
        );
    }

    // Legend
    let _ = write!(
        svg,
        r#"<circle cx="820" cy="92" r="4" fill="#85c1e9" fill-opacity="0.5"/><text x="830" y="96" font-size="12" fill="#34495e">Option2 pool</text>"#
    );
    let _ = write!(
        svg,
        r#"<circle cx="820" cy="112" r="4" fill="#f0b27a" fill-opacity="0.5"/><text x="830" y="116" font-size="12" fill="#34495e">GA pool</text>"#
    );
    let _ = write!(
        svg,
        r#"<circle cx="820" cy="132" r="6" fill="#2e86de" stroke="#1a5276" stroke-width="2"/><text x="832" y="136" font-size="12" fill="#34495e">Selected Option2</text>"#
    );
    let _ = write!(
        svg,
        r#"<circle cx="820" cy="152" r="6" fill="#d35400" stroke="#873600" stroke-width="2"/><text x="832" y="156" font-size="12" fill="#34495e">Selected GA</text>"#
    );

    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

pub(super) fn write_multifidelity_figure(
    path: &Path,
    rows: &[ValidationRow],
    fast_mode: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    if rows.is_empty() {
        let note = if fast_mode {
            "FAST mode skipped multi-fidelity validation."
        } else {
            "No multi-fidelity validation rows available."
        };
        return write_placeholder(path, "Multi-Fidelity Pressure-Drop Comparison", note);
    }
    let mut data: Vec<(String, f64)> = Vec::new();
    for (i, row) in rows.iter().enumerate() {
        let label_base = if i % 2 == 0 { "A" } else { "B" };
        data.push((format!("{label_base}-1D"), row.dp_1d_bernoulli_pa.max(0.0)));
        data.push((format!("{label_base}-2D"), row.dp_2d_fvm_pa.max(0.0)));
        data.push((format!("{label_base}-3D"), row.dp_3d_fem_pa.max(0.0)));
    }
    let max = data.iter().map(|(_, v)| *v).fold(1.0, f64::max);
    write_bar_svg_owned(
        path,
        "Multi-Fidelity Pressure-Drop Comparison",
        &data,
        max * 1.1,
    )
}

pub(super) fn write_ga_convergence_figure(
    path: &Path,
    best_per_gen: &[f64],
    fast_mode: bool,
) -> Result<(), Box<dyn std::error::Error>> {
    if best_per_gen.is_empty() {
        let note = if fast_mode {
            "FAST mode omits GA convergence history."
        } else {
            "GA convergence history not available."
        };
        return write_placeholder(path, "GA Fitness Convergence", note);
    }
    let mut svg = String::new();
    let w = 1100.0;
    let h = 700.0;
    svg_start(&mut svg, w, h);
    svg_title(&mut svg, "GA Fitness Convergence");
    let x0 = 100.0;
    let y0 = 620.0;
    let xw = 900.0;
    let yh = 500.0;
    axis(&mut svg, x0, y0, xw, yh);

    let n = best_per_gen.len().max(2);
    let mut points = String::new();
    for (i, score) in best_per_gen.iter().enumerate() {
        let x = x0 + xw * (i as f64 / (n - 1) as f64);
        let y = y0 - yh * score.clamp(0.0, 1.0);
        let _ = write!(points, "{x:.2},{y:.2} ");
    }
    let _ = write!(
        svg,
        r##"<polyline fill="none" stroke="#8e44ad" stroke-width="3" points="{points}"/>"##
    );
    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

fn write_bar_svg(
    path: &Path,
    title: &str,
    data: &[(&str, f64)],
    y_max: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    if data.is_empty() {
        return write_placeholder(path, title, "No data available.");
    }
    let mut svg = String::new();
    let w = 1100.0;
    let h = 700.0;
    svg_start(&mut svg, w, h);
    svg_title(&mut svg, title);

    let x0 = 100.0;
    let y0 = 620.0;
    let xw = 900.0;
    let yh = 500.0;
    axis(&mut svg, x0, y0, xw, yh);

    let bw = xw / data.len() as f64 * 0.65;
    let gap = xw / data.len() as f64 * 0.35;
    for (i, (label, value)) in data.iter().enumerate() {
        let left = x0 + (bw + gap) * i as f64 + gap * 0.5;
        let bar_h = if y_max <= 1e-12 {
            0.0
        } else {
            yh * (*value / y_max).clamp(0.0, 1.0)
        };
        let top = y0 - bar_h;
        let color = if i % 2 == 0 { "#2e86de" } else { "#27ae60" };
        let _ = write!(
            svg,
            r#"<rect x="{left:.2}" y="{top:.2}" width="{bw:.2}" height="{bar_h:.2}" fill="{color}" fill-opacity="0.8"/>"#
        );
        let _ = write!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" font-size="13" fill="#2c3e50">{label}</text>"##,
            left,
            y0 + 24.0
        );
        let _ = write!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" font-size="12" fill="#2c3e50">{:.3}</text>"##,
            left,
            top - 8.0,
            value
        );
    }
    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

fn write_bar_svg_owned(
    path: &Path,
    title: &str,
    data: &[(String, f64)],
    y_max: f64,
) -> Result<(), Box<dyn std::error::Error>> {
    let refs: Vec<(&str, f64)> = data.iter().map(|(k, v)| (k.as_str(), *v)).collect();
    write_bar_svg(path, title, &refs, y_max)
}

fn svg_start(svg: &mut String, width: f64, height: f64) {
    let _ = write!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0}" height="{height:.0}" viewBox="0 0 {width:.0} {height:.0}">"#
    );
    svg.push_str(r#"<rect x="0" y="0" width="100%" height="100%" fill="white"/>"#);
    svg.push_str(
        r##"<defs><marker id="arrowhead-flow" markerWidth="10" markerHeight="8" refX="8" refY="4" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,8 L10,4 z" fill="#566573"/></marker></defs>"##,
    );
}

fn svg_end(svg: &mut String) {
    svg.push_str("</svg>");
}

fn svg_title(svg: &mut String, title: &str) {
    let _ = write!(
        svg,
        r##"<text x="70" y="60" font-size="28" font-weight="600" fill="#1f2d3d">{}</text>"##,
        escape_xml(title)
    );
}

fn axis(svg: &mut String, x0: f64, y0: f64, xw: f64, yh: f64) {
    let _ = write!(
        svg,
        r##"<line x1="{x0:.2}" y1="{y0:.2}" x2="{:.2}" y2="{y0:.2}" stroke="#2c3e50" stroke-width="2"/>"##,
        x0 + xw
    );
    let _ = write!(
        svg,
        r##"<line x1="{x0:.2}" y1="{y0:.2}" x2="{x0:.2}" y2="{:.2}" stroke="#2c3e50" stroke-width="2"/>"##,
        y0 - yh
    );
}

fn escape_xml(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
}

fn process_box(
    svg: &mut String,
    x: f64,
    y: f64,
    width: f64,
    height: f64,
    fill: &str,
    stroke: &str,
    title: &str,
    lines: &[&str],
) {
    let _ = write!(
        svg,
        r##"<rect x="{x:.1}" y="{y:.1}" width="{width:.1}" height="{height:.1}" rx="16" ry="16" fill="{fill}" stroke="{stroke}" stroke-width="2.5"/>"##
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" font-size="19" font-weight="600" fill="#1f2d3d">{}</text>"##,
        x + 16.0,
        y + 30.0,
        escape_xml(title)
    );
    for (idx, line) in lines.iter().enumerate() {
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" font-size="14" fill="#34495e">{}</text>"##,
            x + 18.0,
            y + 62.0 + idx as f64 * 24.0,
            escape_xml(line)
        );
    }
}

fn arrow(svg: &mut String, x1: f64, y1: f64, x2: f64, y2: f64) {
    let _ = write!(
        svg,
        r##"<line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" stroke="#566573" stroke-width="3" marker-end="url(#arrowhead-flow)"/>"##
    );
}

#[cfg(test)]
mod tests {
    use super::write_creation_optimization_process_figure;
    use std::time::{SystemTime, UNIX_EPOCH};

    #[test]
    fn creation_optimization_process_figure_is_generated_not_placeholder() {
        let nanos = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("system clock before unix epoch")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("m12_process_{nanos}.svg"));
        write_creation_optimization_process_figure(&path)
            .expect("process figure generation must succeed");
        let svg = std::fs::read_to_string(path).expect("must read generated process svg");
        assert!(svg.contains("Primitive candidate generation"));
        assert!(svg.contains("Cached cfd-1d evaluation"));
        assert!(svg.contains("Option 2 combined ranking"));
        assert!(!svg.contains("placeholder"));
    }
}
