//! SVG figure builders for Milestone 12 narrative figures.

use std::fmt::Write as _;
use std::io::Write as IoWrite;
use std::path::Path;

use super::primitives::{axis, svg_end, svg_start, svg_title, write_bar_svg, write_bar_svg_owned};
use super::process::write_placeholder;
use crate::reporting::{Milestone12ReportDesign, ParetoPoint};

/// Adapter that bridges `std::io::Write` → `std::fmt::Write`, enabling
/// `write!()` macros (which use `fmt::Write`) to stream directly to a
/// `BufWriter<File>` without materialising the entire SVG in memory.
struct FmtToIo<W>(W);

impl<W: IoWrite> std::fmt::Write for FmtToIo<W> {
    fn write_str(&mut self, s: &str) -> std::fmt::Result {
        self.0.write_all(s.as_bytes()).map_err(|_| std::fmt::Error)
    }
}

pub(super) fn write_cross_mode_figure(
    path: &Path,
    option1: &[Milestone12ReportDesign],
    option2: &[Milestone12ReportDesign],
    ga: &[Milestone12ReportDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let data = vec![
        ("Opt1 Acoustic", option1.first().map_or(0.0, |d| d.score)),
        ("Opt2 Combined", option2.first().map_or(0.0, |d| d.score)),
        ("GA HydroSDT", ga.first().map_or(0.0, |d| d.score)),
    ];
    write_bar_svg(
        path,
        "Cross-Mode Scoring Comparison",
        &data,
        1.0,
        "Optimization Strategy",
        "Composite Score",
    )
}

pub(super) fn write_head_to_head_figure(
    path: &Path,
    option2: &[Milestone12ReportDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    if option2.is_empty() {
        return write_placeholder(
            path,
            "Head-to-Head Design Comparison",
            "No Option 2 ranked data available.",
        );
    }
    // Single collect into owned (String, f64) — use write_bar_svg_owned to
    // avoid the second Vec<(&str, f64)> allocation.
    let data: Vec<(String, f64)> = option2
        .iter()
        .take(5)
        .map(|d| (format!("R{}", d.rank), d.score))
        .collect();
    write_bar_svg_owned(
        path,
        "Head-to-Head Option 2 Top-5 Scores",
        &data,
        1.0,
        "Design Rank",
        "Composite Score",
    )
}

pub(super) fn write_cavitation_distribution_figure(
    path: &Path,
    option2: &[Milestone12ReportDesign],
    ga: &[Milestone12ReportDesign],
) -> Result<(), Box<dyn std::error::Error>> {
    let data: Vec<(String, f64, &str)> = option2
        .iter()
        .take(5)
        .map(|d| (format!("O2-R{}", d.rank), d.metrics.cavitation_number, "Option2"))
        .chain(
            ga.iter()
                .take(5)
                .map(|d| (format!("GA-R{}", d.rank), d.metrics.cavitation_number, "GA")),
        )
        .filter(|(_, sigma, _)| sigma.is_finite())
        .collect();

    if data.is_empty() {
        return write_placeholder(
            path,
            "Selected-Design Cavitation Number (σ)",
            "No finite cavitation numbers available.",
        );
    }

    let sigma_min = data.iter().map(|(_, sigma, _)| *sigma).fold(f64::INFINITY, f64::min);
    let sigma_max = data
        .iter()
        .map(|(_, sigma, _)| *sigma)
        .fold(f64::NEG_INFINITY, f64::max);
    let y_min = sigma_min.min(0.0) - 0.1 * (sigma_max - sigma_min).max(0.25);
    let y_max = sigma_max.max(1.0) + 0.1 * (sigma_max - sigma_min).max(0.25);

    let mut svg = String::new();
    let w = 1100.0;
    let h = 700.0;
    svg_start(&mut svg, w, h);
    svg_title(&mut svg, "Selected-Design Cavitation Number (σ)");
    let x0 = 100.0;
    let y0 = 620.0;
    let xw = 900.0;
    let yh = 500.0;
    axis(&mut svg, x0, y0, xw, yh);

    for i in 0..=5 {
        let frac = f64::from(i) / 5.0;
        let val = y_min + frac * (y_max - y_min);
        let ty = y0 - yh * frac;
        let _ = write!(
            svg,
            r##"<line x1="{x0:.2}" y1="{ty:.2}" x2="{:.2}" y2="{ty:.2}" stroke="#ecf0f1" stroke-width="1"/>"##,
            x0 + xw
        );
        let _ = write!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" font-size="12" text-anchor="end" fill="#7f8c8d">{:.2}</text>"##,
            x0 - 8.0,
            ty + 4.0,
            val
        );
    }

    for &(threshold, color) in &[(0.0, "#c0392b"), (1.0, "#d4ac0d")] {
        if threshold >= y_min && threshold <= y_max {
            let frac = (threshold - y_min) / (y_max - y_min);
            let ty = y0 - yh * frac;
            let _ = write!(
                svg,
                r##"<line x1="{x0:.2}" y1="{ty:.2}" x2="{:.2}" y2="{ty:.2}" stroke="{color}" stroke-width="2" stroke-dasharray="8 6"/>"##,
                x0 + xw
            );
            let _ = write!(
                svg,
                r##"<text x="{:.2}" y="{:.2}" font-size="12" fill="{color}">σ={threshold:.0}</text>"##,
                x0 + xw - 42.0,
                ty - 6.0
            );
        }
    }

    let step = xw / data.len() as f64;
    for (index, (label, sigma, tag)) in data.iter().enumerate() {
        let cx = x0 + step * (index as f64 + 0.5);
        let cy = y0 - yh * ((*sigma - y_min) / (y_max - y_min)).clamp(0.0, 1.0);
        let (fill, stroke) = if *tag == "Option2" {
            ("#2e86de", "#1a5276")
        } else {
            ("#d35400", "#873600")
        };
        let _ = write!(
            svg,
            r#"<circle cx="{cx:.2}" cy="{cy:.2}" r="8" fill="{fill}" stroke="{stroke}" stroke-width="2"/>"#
        );
        let _ = write!(
            svg,
            r##"<text x="{cx:.2}" y="{:.2}" font-size="12" text-anchor="middle" fill="#2c3e50">{}</text>"##,
            y0 + 24.0,
            super::primitives::escape_xml(label)
        );
        let _ = write!(
            svg,
            r##"<text x="{cx:.2}" y="{:.2}" font-size="12" text-anchor="middle" fill="#2c3e50">{sigma:.3}</text>"##,
            cy - 12.0
        );
    }

    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="668" font-size="16" fill="#34495e">Selected Venturi Designs</text>"##,
        x0 + xw * 0.38
    );
    svg.push_str(r##"<text x="18" y="340" transform="rotate(-90 18,340)" font-size="16" fill="#34495e">Cavitation Number (σ)</text>"##);
    svg.push_str(r##"<circle cx="820" cy="92" r="6" fill="#2e86de" stroke="#1a5276" stroke-width="2"/><text x="832" y="96" font-size="12" fill="#34495e">Option 2</text>"##);
    svg.push_str(r##"<circle cx="820" cy="112" r="6" fill="#d35400" stroke="#873600" stroke-width="2"/><text x="832" y="116" font-size="12" fill="#34495e">GA</text>"##);
    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

pub(super) fn write_pareto_figure(
    path: &Path,
    option2: &[Milestone12ReportDesign],
    ga: &[Milestone12ReportDesign],
    _option2_pool_all: &[ParetoPoint],
    _ga_pool_all: &[ParetoPoint],
) -> Result<(), Box<dyn std::error::Error>> {
    let selected: Vec<(String, f64, f64, f64, &str)> = option2
        .iter()
        .take(5)
        .map(|d| {
            (
                format!("O2-R{}", d.rank),
                d.metrics.cancer_targeted_cavitation,
                d.metrics.rbc_venturi_protection,
                d.score,
                "Option2",
            )
        })
        .chain(ga.iter().take(5).map(|d| {
            (
                format!("GA-R{}", d.rank),
                d.metrics.cancer_targeted_cavitation,
                d.metrics.rbc_venturi_protection,
                d.score,
                "GA",
            )
        }))
        .filter(|(_, x, y, _, _)| x.is_finite() && y.is_finite())
        .collect();

    if selected.is_empty() {
        return write_placeholder(
            path,
            "Selected-Design Oncology Trade-Off Frontier",
            "No ranked data available.",
        );
    }

    let (x_min, x_max, y_min, y_max) = {
        let (mut xmn, mut xmx, mut ymn, mut ymx) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for (_, x, y, _, _) in &selected {
            xmn = xmn.min(*x);
            xmx = xmx.max(*x);
            ymn = ymn.min(*y);
            ymx = ymx.max(*y);
        }
        let x_span = (xmx - xmn).max(0.01);
        let y_span = (ymx - ymn).max(0.01);
        (
            xmn - 0.12 * x_span,
            xmx + 0.12 * x_span,
            ymn - 0.12 * y_span,
            ymx + 0.12 * y_span,
        )
    };

    let frontier: Vec<(f64, f64)> = {
        let mut nondominated = Vec::new();
        'outer: for (_, x, y, _, _) in &selected {
            for (_, ox, oy, _, _) in &selected {
                if (ox > x || oy > y) && ox >= x && oy >= y {
                    continue 'outer;
                }
            }
            nondominated.push((*x, *y));
        }
        nondominated.sort_by(|left, right| left.0.total_cmp(&right.0));
        nondominated.dedup_by(|left, right| {
            (left.0 - right.0).abs() <= 1.0e-9 && (left.1 - right.1).abs() <= 1.0e-9
        });
        nondominated
    };

    let file = std::fs::File::create(path)?;
    let mut svg = FmtToIo(std::io::BufWriter::new(file));
    let w = 1100.0;
    let h = 700.0;
    svg_start(&mut svg, w, h);
    svg_title(&mut svg, "Selected-Design Oncology Trade-Off Frontier");
    let x0 = 100.0;
    let y0 = 620.0;
    let xw = 900.0;
    let yh = 500.0;
    axis(&mut svg, x0, y0, xw, yh);

    // Tick labels — X axis
    for i in 0..=5 {
        let frac = f64::from(i) / 5.0;
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
        let frac = f64::from(i) / 5.0;
        let val = y_min + frac * (y_max - y_min);
        let ty = y0 - yh * frac;
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{ty:.1}" text-anchor="end" font-size="12" fill="#7f8c8d">{val:.3}</text>"##,
            x0 - 8.0
        );
    }

    svg.write_str(r##"<text x="420" y="668" font-size="16" fill="#34495e">Tumor-Targeted Cavitation Index</text>"##)?;
    svg.write_str(r##"<text x="18" y="340" transform="rotate(-90 18,340)" font-size="16" fill="#34495e">Healthy-Cell Protection (RBC Venturi)</text>"##)?;

    if frontier.len() >= 2 {
        let polyline = frontier
            .iter()
            .map(|(cx, cy)| {
                let x = x0 + xw * ((*cx - x_min) / (x_max - x_min)).clamp(0.0, 1.0);
                let y = y0 - yh * ((*cy - y_min) / (y_max - y_min)).clamp(0.0, 1.0);
                format!("{x:.2},{y:.2}")
            })
            .collect::<Vec<_>>()
            .join(" ");
        let _ = write!(
            svg,
            r##"<polyline points="{polyline}" fill="none" stroke="#7d3c98" stroke-width="3" stroke-linejoin="round" stroke-linecap="round"/>"##
        );
    }

    for (label, cx, cy, score, tag) in &selected {
        let x = x0 + xw * ((cx - x_min) / (x_max - x_min)).clamp(0.0, 1.0);
        let y = y0 - yh * ((cy - y_min) / (y_max - y_min)).clamp(0.0, 1.0);
        let r = 6.0 + 6.0 * score.clamp(0.0, 1.0);
        let (fill, stroke) = if *tag == "Option2" {
            ("#2e86de", "#1a5276")
        } else {
            ("#d35400", "#873600")
        };
        let _ = write!(
            svg,
            r#"<circle cx="{x:.2}" cy="{y:.2}" r="{r:.2}" fill="{fill}" fill-opacity="0.85" stroke="{stroke}" stroke-width="2"/>"#
        );
        let _ = write!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" font-size="12" fill="#2c3e50">{}</text>"##,
            x + 10.0,
            y - 10.0,
            super::primitives::escape_xml(label)
        );
    }

    let _ = write!(
        svg,
        r##"<line x1="820" y1="90" x2="846" y2="90" stroke="#7d3c98" stroke-width="3"/><text x="854" y="94" font-size="12" fill="#34495e">Frontier</text>"##
    );
    let _ = write!(
        svg,
        r##"<circle cx="820" cy="112" r="6" fill="#2e86de" stroke="#1a5276" stroke-width="2"/><text x="832" y="116" font-size="12" fill="#34495e">Option 2</text>"##
    );
    let _ = write!(
        svg,
        r##"<circle cx="820" cy="132" r="6" fill="#d35400" stroke="#873600" stroke-width="2"/><text x="832" y="136" font-size="12" fill="#34495e">GA</text>"##
    );

    svg_end(&mut svg);
    svg.0.flush()?;
    Ok(())
}

pub(super) fn write_pediatric_ecv_figure(
    path: &Path,
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
    ga: &Milestone12ReportDesign,
) -> Result<(), Box<dyn std::error::Error>> {
    const PEDIATRIC_LIMIT_ML: f64 = 25.5;

    let mut data = Vec::new();
    if let Some(option1) = option1 {
        data.push((
            "Opt1 Acoustic".to_string(),
            100.0 * option1.metrics.total_ecv_ml / PEDIATRIC_LIMIT_ML,
        ));
    }
    data.push((
        "Opt2 Combined".to_string(),
        100.0 * option2.metrics.total_ecv_ml / PEDIATRIC_LIMIT_ML,
    ));
    data.push((
        "GA HydroSDT".to_string(),
        100.0 * ga.metrics.total_ecv_ml / PEDIATRIC_LIMIT_ML,
    ));
    let y_max = data
        .iter()
        .map(|(_, value)| *value)
        .fold(100.0_f64, f64::max)
        * 1.15;

    write_bar_svg_owned(
        path,
        "Pediatric Circuit Volume Margin (3 kg Reference)",
        &data,
        y_max.max(100.0),
        "Selected Design",
        "% of 25.5 mL Limit",
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
    // Stream SVG directly to file to avoid ~1 MB in-memory string.
    let file = std::fs::File::create(path)?;
    let mut svg = FmtToIo(std::io::BufWriter::new(file));
    let w = 1100.0;
    let h = 700.0;
    svg_start(&mut svg, w, h);
    svg_title(&mut svg, "GA Fitness Convergence");
    let x0 = 100.0;
    let y0 = 620.0;
    let xw = 900.0;
    let yh = 500.0;
    axis(&mut svg, x0, y0, xw, yh);

    // Auto-scale Y axis to data range with 10% margin
    let raw_min = best_per_gen
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::INFINITY, f64::min);
    let raw_max = best_per_gen
        .iter()
        .copied()
        .filter(|v| v.is_finite())
        .fold(f64::NEG_INFINITY, f64::max);
    let (y_lo, y_hi) =
        if !raw_min.is_finite() || !raw_max.is_finite() || (raw_max - raw_min).abs() < 1e-12 {
            let center = if raw_min.is_finite() { raw_min } else { 0.0 };
            (center - 0.05, center + 0.05)
        } else {
            let span = raw_max - raw_min;
            (raw_min - 0.1 * span, raw_max + 0.1 * span)
        };

    // Y-axis tick labels
    for i in 0..=5 {
        let frac = f64::from(i) / 5.0;
        let val = y_lo + frac * (y_hi - y_lo);
        let ty = y0 - yh * frac;
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{ty:.1}" text-anchor="end" font-size="12" fill="#7f8c8d">{val:.4}</text>"##,
            x0 - 8.0
        );
        // Grid line
        let _ = write!(
            svg,
            r##"<line x1="{x0:.1}" y1="{ty:.1}" x2="{:.1}" y2="{ty:.1}" stroke="#ecf0f1" stroke-width="1"/>"##,
            x0 + xw
        );
    }
    // X-axis tick labels (generation numbers)
    let n_gen = best_per_gen.len();
    let x_ticks = 5.min(n_gen.saturating_sub(1)).max(1);
    for i in 0..=x_ticks {
        let gen_idx = if x_ticks == 0 {
            0
        } else {
            i * (n_gen - 1) / x_ticks
        };
        let frac = if n_gen <= 1 {
            0.5
        } else {
            gen_idx as f64 / (n_gen - 1) as f64
        };
        let tx = x0 + xw * frac;
        let _ = write!(
            svg,
            r##"<text x="{tx:.1}" y="{:.1}" text-anchor="middle" font-size="12" fill="#7f8c8d">{gen_idx}</text>"##,
            y0 + 18.0
        );
    }

    // Axis labels
    svg.write_str(
        r##"<text x="500" y="668" font-size="16" fill="#34495e">Generation Number</text>"##,
    )?;
    svg.write_str(r##"<text x="18" y="360" transform="rotate(-90 18,360)" font-size="16" fill="#34495e">Best Fitness Score (HydroSDT Composite)</text>"##)?;

    // Plot convergence curve
    let n = best_per_gen.len().max(2);
    let mut points = String::new();
    for (i, score) in best_per_gen.iter().enumerate() {
        let x = x0 + xw * (i as f64 / (n - 1) as f64);
        let s = if score.is_finite() { *score } else { y_lo };
        let y_frac = ((s - y_lo) / (y_hi - y_lo)).clamp(0.0, 1.0);
        let y = y0 - yh * y_frac;
        let _ = write!(points, "{x:.2},{y:.2} ");
    }
    let _ = write!(
        svg,
        r##"<polyline fill="none" stroke="#8e44ad" stroke-width="3" points="{points}"/>"##
    );

    // Dot markers on each generation
    for (i, score) in best_per_gen.iter().enumerate() {
        let x = x0 + xw * (i as f64 / (n - 1) as f64);
        let s = if score.is_finite() { *score } else { y_lo };
        let y_frac = ((s - y_lo) / (y_hi - y_lo)).clamp(0.0, 1.0);
        let y = y0 - yh * y_frac;
        let _ = write!(
            svg,
            r##"<circle cx="{x:.2}" cy="{y:.2}" r="4" fill="#8e44ad" fill-opacity="0.8"/>"##
        );
    }

    svg_end(&mut svg);
    svg.0.flush()?;
    Ok(())
}
