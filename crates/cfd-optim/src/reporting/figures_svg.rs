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

pub(super) fn write_cross_mode_figure(
    path: &Path,
    option1: &[RankedDesign],
    option2: &[RankedDesign],
    rbc: &[RankedDesign],
    ga: &[RankedDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let data = vec![
        ("Option1", option1.first().map_or(0.0, |d| d.score)),
        ("Option2", option2.first().map_or(0.0, |d| d.score)),
        ("RBC", rbc.first().map_or(0.0, |d| d.score)),
        ("GA", ga.first().map_or(0.0, |d| d.score)),
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
    rbc: &[RankedDesign],
    ga: &[RankedDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut neg = 0usize;
    let mut zero_to_one = 0usize;
    let mut ge_one = 0usize;
    let mut nonfinite = 0usize;

    for d in option2.iter().chain(rbc.iter()).chain(ga.iter()) {
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
    rbc: &[RankedDesign],
    ga: &[RankedDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let points: Vec<(&RankedDesign, &str)> = option2
        .iter()
        .map(|d| (d, "Option2"))
        .chain(rbc.iter().map(|d| (d, "RBC")))
        .chain(ga.iter().map(|d| (d, "GA")))
        .collect();
    if points.is_empty() {
        return write_placeholder(
            path,
            "Pareto Front — Oncology Objectives",
            "No ranked data available.",
        );
    }

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
    svg.push_str(r##"<text x="460" y="660" font-size="16" fill="#34495e">cancer_targeted_cavitation</text>"##);
    svg.push_str(r##"<text x="18" y="340" transform="rotate(-90 18,340)" font-size="16" fill="#34495e">therapeutic_window_score</text>"##);

    for (d, tag) in points {
        let x = x0 + xw * d.metrics.cancer_targeted_cavitation.clamp(0.0, 1.0);
        let y = y0 - yh * d.metrics.therapeutic_window_score.clamp(0.0, 1.0);
        let r = 4.0 + 6.0 * d.score.clamp(0.0, 1.0);
        let color = match tag {
            "Option2" => "#2e86de",
            "RBC" => "#27ae60",
            _ => "#d35400",
        };
        let _ = write!(
            svg,
            r#"<circle cx="{x:.2}" cy="{y:.2}" r="{r:.2}" fill="{color}" fill-opacity="0.65"/>"#
        );
    }
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
