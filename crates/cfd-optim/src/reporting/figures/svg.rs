//! SVG figure builders for Milestone 12 narrative figures.

use std::fmt::Write as _;
use std::io::Write as IoWrite;
use std::path::Path;

use super::primitives::{axis, svg_end, svg_start, svg_title, write_bar_svg, write_bar_svg_owned};
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
        return Err(
            "cannot render head-to-head Option 2 figure without ranked Option 2 data".into(),
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
        .enumerate()
        .map(|(idx, d)| {
            (
                format!("O2-R{}", idx + 1),
                d.metrics.cavitation_number,
                "Option2",
            )
        })
        .chain(ga.iter().take(5).enumerate().map(|(idx, d)| {
            (
                format!("GA-R{}", idx + 1),
                d.metrics.cavitation_number,
                "GA",
            )
        }))
        .filter(|(_, sigma, _)| sigma.is_finite())
        .collect();

    if data.is_empty() {
        return Err(
            "cannot render cavitation distribution without finite cavitation numbers".into(),
        );
    }

    let sigma_positive_max = data
        .iter()
        .map(|(_, sigma, _)| *sigma)
        .filter(|sigma| *sigma >= 1.0)
        .fold(1.0_f64, f64::max);

    let map_sigma_to_plot = |sigma: f64| -> f64 {
        if sigma < 1.0 {
            sigma.clamp(-2.0, 1.0)
        } else {
            1.0 + sigma
                .log10()
                .clamp(0.0, sigma_positive_max.log10().max(0.1))
        }
    };

    let plot_min = -2.0;
    let plot_max = 1.0 + sigma_positive_max.log10().max(0.1);
    let y_min = plot_min - 0.15 * (plot_max - plot_min).max(1.0);
    let y_max = plot_max + 0.10 * (plot_max - plot_min).max(1.0);

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
        let label = if val <= 1.0 {
            format!("{val:.2}")
        } else {
            format!("1 + log10(σ) = {val:.2}")
        };
        let _ = write!(
            svg,
            r##"<line x1="{x0:.2}" y1="{ty:.2}" x2="{:.2}" y2="{ty:.2}" stroke="#ecf0f1" stroke-width="1"/>"##,
            x0 + xw
        );
        let _ = write!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" font-size="12" text-anchor="end" fill="#7f8c8d">{}</text>"##,
            x0 - 8.0,
            ty + 4.0,
            super::primitives::escape_xml(&label)
        );
    }

    for &(threshold, color) in &[(0.0, "#c0392b"), (1.0, "#d4ac0d")] {
        let plot_threshold = map_sigma_to_plot(threshold);
        if plot_threshold >= y_min && plot_threshold <= y_max {
            let frac = (plot_threshold - y_min) / (y_max - y_min);
            let ty = y0 - yh * frac;
            let _ = write!(
                svg,
                r#"<line x1="{x0:.2}" y1="{ty:.2}" x2="{:.2}" y2="{ty:.2}" stroke="{color}" stroke-width="2" stroke-dasharray="8 6"/>"#,
                x0 + xw
            );
            let _ = write!(
                svg,
                r#"<text x="{:.2}" y="{:.2}" font-size="12" fill="{color}">σ={threshold:.0}</text>"#,
                x0 + xw - 42.0,
                ty - 6.0
            );
        }
    }

    let step = xw / data.len() as f64;
    for (index, (label, sigma, tag)) in data.iter().enumerate() {
        let cx = x0 + step * (index as f64 + 0.5);
        let plot_sigma = map_sigma_to_plot(*sigma);
        let cy = y0 - yh * ((plot_sigma - y_min) / (y_max - y_min)).clamp(0.0, 1.0);
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
            r##"<text x="{cx:.2}" y="{:.2}" font-size="12" text-anchor="middle" fill="#2c3e50">{}</text>"##,
            cy - 12.0,
            if *sigma < 1.0 {
                format!("{sigma:.3}")
            } else {
                format!(">1 ({sigma:.1})")
            }
        );
    }

    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="668" font-size="16" fill="#34495e">Selected Venturi Designs</text>"##,
        x0 + xw * 0.38
    );
    svg.push_str(r##"<text x="18" y="340" transform="rotate(-90 18,340)" font-size="16" fill="#34495e">σ (linear below 1, log-scaled above 1)</text>"##);
    svg.push_str(r##"<text x="110" y="86" font-size="12" fill="#7f8c8d">Raw σ labels are shown above each point; values above 1 are plotted on a log-compressed non-cavitating branch.</text>"##);
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
    option2_pool_all: &[ParetoPoint],
    ga_pool_all: &[ParetoPoint],
) -> Result<(), Box<dyn std::error::Error>> {
    let mut selected: Vec<(String, f64, f64, f64, &str)> = option2
        .iter()
        .take(5)
        .enumerate()
        .map(|(idx, d)| {
            (
                format!("O2-R{}", idx + 1),
                d.metrics.cancer_targeted_cavitation,
                d.metrics.healthy_cell_protection_index,
                d.score,
                "Option2",
            )
        })
        .chain(ga.iter().take(5).enumerate().map(|(idx, d)| {
            (
                format!("GA-R{}", idx + 1),
                d.metrics.cancer_targeted_cavitation,
                d.metrics.healthy_cell_protection_index,
                d.score,
                "GA",
            )
        }))
        .filter(|(_, x, y, _, _)| x.is_finite() && y.is_finite())
        .collect();
    selected.sort_by(|left, right| {
        left.0
            .cmp(&right.0)
            .then_with(|| left.1.total_cmp(&right.1))
            .then_with(|| left.2.total_cmp(&right.2))
            .then_with(|| left.3.total_cmp(&right.3))
            .then_with(|| left.4.cmp(right.4))
    });

    let mut background: Vec<(f64, f64, f64, &str)> = option2_pool_all
        .iter()
        .filter(|p| {
            p.cancer_targeted_cavitation.is_finite() && p.healthy_cell_protection_index.is_finite()
        })
        .map(|p| {
            (
                p.cancer_targeted_cavitation,
                p.healthy_cell_protection_index,
                p.score,
                "Option2",
            )
        })
        .chain(
            ga_pool_all
                .iter()
                .filter(|p| {
                    p.cancer_targeted_cavitation.is_finite()
                        && p.healthy_cell_protection_index.is_finite()
                })
                .map(|p| {
                    (
                        p.cancer_targeted_cavitation,
                        p.healthy_cell_protection_index,
                        p.score,
                        "GA",
                    )
                }),
        )
        .collect();
    background.sort_by(|left, right| {
        left.0
            .total_cmp(&right.0)
            .then_with(|| left.1.total_cmp(&right.1))
            .then_with(|| left.2.total_cmp(&right.2))
            .then_with(|| left.3.cmp(right.3))
    });

    if selected.is_empty() && background.is_empty() {
        return Err(
            "cannot render Pareto frontier without selected or background ranked data".into(),
        );
    }

    let (x_min, x_max, y_min, y_max) = {
        let (mut xmn, mut xmx, mut ymn, mut ymx) = (
            f64::INFINITY,
            f64::NEG_INFINITY,
            f64::INFINITY,
            f64::NEG_INFINITY,
        );
        for (x, y, _, _) in &background {
            xmn = xmn.min(*x);
            xmx = xmx.max(*x);
            ymn = ymn.min(*y);
            ymx = ymx.max(*y);
        }
        for (_, x, y, _, _) in &selected {
            xmn = xmn.min(*x);
            xmx = xmx.max(*x);
            ymn = ymn.min(*y);
            ymx = ymx.max(*y);
        }
        let x_span = if xmx > xmn { xmx - xmn } else { 1.0 };
        let y_span = if ymx > ymn { ymx - ymn } else { 1.0 };
        (
            xmn - 0.12 * x_span,
            xmx + 0.12 * x_span,
            ymn - 0.12 * y_span,
            ymx + 0.12 * y_span,
        )
    };

    let frontier: Vec<(f64, f64)> = {
        let mut nondominated = Vec::new();
        'outer: for (x, y, _, _) in &background {
            for (ox, oy, _, _) in &background {
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
    svg.write_str(r##"<text x="18" y="340" transform="rotate(-90 18,340)" font-size="16" fill="#34495e">Healthy-Cell Protection Index</text>"##)?;

    if !background.is_empty() {
        for (cx, cy, score, tag) in &background {
            let x = x0 + xw * ((*cx - x_min) / (x_max - x_min)).clamp(0.0, 1.0);
            let y = y0 - yh * ((*cy - y_min) / (y_max - y_min)).clamp(0.0, 1.0);
            let r = 2.0 + 1.5 * score.clamp(0.0, 1.0);
            let fill = if *tag == "Option2" {
                "#85c1e9"
            } else {
                "#f5b041"
            };
            let _ = write!(
                svg,
                r#"<circle cx="{x:.2}" cy="{y:.2}" r="{r:.2}" fill="{fill}" fill-opacity="0.20" stroke="none"/>"#
            );
        }
    }

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
        r##"<line x1="820" y1="90" x2="846" y2="90" stroke="#7d3c98" stroke-width="3"/><text x="854" y="94" font-size="12" fill="#34495e">Frontier (full pool)</text>"##
    );
    let _ = write!(
        svg,
        r##"<circle cx="820" cy="102" r="4" fill="#85c1e9" fill-opacity="0.35"/><text x="832" y="106" font-size="12" fill="#34495e">Background pool</text>"##
    );
    let _ = write!(
        svg,
        r##"<circle cx="820" cy="122" r="6" fill="#2e86de" stroke="#1a5276" stroke-width="2"/><text x="832" y="126" font-size="12" fill="#34495e">Option 2 selected</text>"##
    );
    let _ = write!(
        svg,
        r##"<circle cx="820" cy="142" r="6" fill="#d35400" stroke="#873600" stroke-width="2"/><text x="832" y="146" font-size="12" fill="#34495e">GA selected</text>"##
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
        let mode = if fast_mode { "fast" } else { "full" };
        return Err(format!(
            "cannot render GA convergence figure for {mode} run without generation history"
        )
        .into());
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
    let (y_lo, y_hi) = if !raw_min.is_finite() || !raw_max.is_finite() || raw_max == raw_min {
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
        let gen_idx = (i * (n_gen - 1)).checked_div(x_ticks).unwrap_or(0);
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

    let tail_len = best_per_gen.len().min(5);
    let tail_start = best_per_gen.len().saturating_sub(tail_len);
    let tail_delta = if tail_len >= 2 {
        best_per_gen[best_per_gen.len() - 1] - best_per_gen[tail_start]
    } else {
        0.0
    };
    let trajectory_note = if tail_delta > 1.0e-3 {
        format!("Still improving: Δbest(last {tail_len} gen) = +{tail_delta:.4}")
    } else if tail_delta < -1.0e-3 {
        format!("Best fitness regressed over last {tail_len} gen by {tail_delta:.4}")
    } else {
        format!(
            "Near-plateau: |Δbest(last {} gen)| = {:.4}",
            tail_len,
            tail_delta.abs()
        )
    };
    let _ = write!(
        svg,
        r##"<text x="690" y="88" font-size="14" fill="#6c3483">{}</text>"##,
        super::primitives::escape_xml(&trajectory_note)
    );

    svg_end(&mut svg);
    svg.0.flush()?;
    Ok(())
}

/// Per-venturi-placement data for the Dean-venturi figure.
pub(super) struct DeanVenturiPoint {
    pub label: String,
    pub dean_number: f64,
    pub cavitation_number: f64,
    pub throat_velocity_m_s: f64,
    /// Inlet-normalized venturi total loss coefficient.
    pub total_loss_coefficient: f64,
    /// Upstream static pressure at this venturi position [kPa].
    pub upstream_pressure_kpa: f64,
    /// Curvature radius at this bend [mm].
    pub bend_radius_mm: f64,
}

/// Generate a figure showing how venturi placement position along a serpentine
/// channel influences Dean number, cavitation strength, and available upstream
/// pressure at each bend.
///
/// - Blue bars (left Y-axis): Dean number.  Alternates high/low because
///   mirrored serpentine U-turns produce tighter inner bends and wider outer
///   bends.
/// - Orange bars (left Y-axis, shared scale): cavitation strength (1 - sigma).
///   Decreases along the path as upstream pressure is consumed by preceding
///   throats.
/// - Dashed green line: upstream static pressure at each throat position [kPa].
///   Shows the serial pressure decay that limits downstream cavitation.
pub(super) fn write_dean_venturi_placement_figure(
    path: &Path,
    title: &str,
    points: &[DeanVenturiPoint],
) -> Result<(), Box<dyn std::error::Error>> {
    if points.is_empty() {
        return Err(format!("cannot render '{title}' without venturi placement data").into());
    }

    let file = std::fs::File::create(path)?;
    let mut svg = FmtToIo(std::io::BufWriter::new(file));
    let w = 1100.0;
    let h = 780.0;
    svg_start(&mut svg, w, h);
    svg_title(&mut svg, title);

    let x0 = 110.0;
    let y0 = 660.0;
    let xw = 820.0;
    let yh = 510.0;
    axis(&mut svg, x0, y0, xw, yh);

    // ---- Data ranges ----
    let de_max = points.iter().map(|p| p.dean_number).fold(1.0_f64, f64::max);
    let cav_strength: Vec<f64> = points
        .iter()
        .map(|p| (1.0 - p.cavitation_number).max(0.0))
        .collect();
    let cs_max = cav_strength.iter().copied().fold(0.0_f64, f64::max);
    // Shared left-axis max: accommodate both De and cavitation strength.
    // Normalise cavitation strength bars relative to De scale so they
    // are visually comparable.
    let left_max = de_max.max(1.0) * 1.20;

    // Pressure range for the right-axis overlay line.
    let p_min = points
        .iter()
        .map(|p| p.upstream_pressure_kpa)
        .fold(f64::INFINITY, f64::min);
    let p_max = points
        .iter()
        .map(|p| p.upstream_pressure_kpa)
        .fold(f64::NEG_INFINITY, f64::max);
    let p_span = (p_max - p_min).max(1.0);
    let p_lo = p_min - 0.10 * p_span;
    let p_hi = p_max + 0.10 * p_span;
    let format_pressure = |value: f64| {
        if p_span < 5.0 {
            format!("{value:.2}")
        } else {
            format!("{value:.0}")
        }
    };

    // ---- Left Y-axis ticks (De / cavitation strength) ----
    for i in 0..=5 {
        let frac = f64::from(i) / 5.0;
        let val = left_max * frac;
        let ty = y0 - yh * frac;
        let _ = write!(
            svg,
            r##"<line x1="{x0:.1}" y1="{ty:.1}" x2="{:.1}" y2="{ty:.1}" stroke="#ecf0f1" stroke-width="1"/>"##,
            x0 + xw
        );
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" font-size="11" text-anchor="end" fill="#2c3e50">{:.0}</text>"##,
            x0 - 8.0,
            ty + 4.0,
            val
        );
    }

    // ---- Right Y-axis ticks (upstream pressure kPa) ----
    let rx = x0 + xw;
    for i in 0..=5 {
        let frac = f64::from(i) / 5.0;
        let val = p_lo + frac * (p_hi - p_lo);
        let ty = y0 - yh * frac;
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" font-size="11" text-anchor="start" fill="#27ae60">{}</text>"##,
            rx + 8.0,
            ty + 4.0,
            format_pressure(val)
        );
    }
    // Right Y axis line (green for pressure)
    let _ = write!(
        svg,
        r##"<line x1="{rx:.1}" y1="{y0:.1}" x2="{rx:.1}" y2="{:.1}" stroke="#27ae60" stroke-width="2"/>"##,
        y0 - yh
    );

    // ---- Grouped bars ----
    let n = points.len();
    let group_w = xw / n as f64;
    let bar_w = group_w * 0.32;
    let gap = group_w * 0.06;

    // Scale factor to render cavitation strength bars on the same left axis.
    // We scale cs values into De units so they share the axis.
    let cs_to_de = if cs_max > 0.0 { de_max / cs_max } else { 1.0 };

    for (i, pt) in points.iter().enumerate() {
        let group_x = x0 + group_w * i as f64;
        let center_x = group_x + group_w * 0.5;

        // ---- Dean number bar (blue) ----
        let de_h = yh * (pt.dean_number / left_max).clamp(0.0, 1.0);
        let de_top = y0 - de_h;
        let de_x = group_x + gap;
        let _ = write!(
            svg,
            r##"<rect x="{de_x:.1}" y="{de_top:.1}" width="{bar_w:.1}" height="{de_h:.1}" fill="#2e86de" fill-opacity="0.85" rx="2"/>"##
        );
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" font-size="10" text-anchor="middle" fill="#1a5276" font-weight="600">{:.0}</text>"##,
            de_x + bar_w * 0.5,
            de_top - 5.0,
            pt.dean_number
        );

        // ---- Cavitation strength bar (orange), scaled to De axis ----
        let cs = cav_strength[i];
        let cs_scaled = cs * cs_to_de;
        let cs_h = yh * (cs_scaled / left_max).clamp(0.0, 1.0);
        let cs_top = y0 - cs_h;
        let cs_x = de_x + bar_w + gap;
        let _ = write!(
            svg,
            r##"<rect x="{cs_x:.1}" y="{cs_top:.1}" width="{bar_w:.1}" height="{cs_h:.1}" fill="#d35400" fill-opacity="0.85" rx="2"/>"##
        );
        // Show sigma value above the bar.
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" font-size="9" text-anchor="middle" fill="#873600">{}</text>"##,
            cs_x + bar_w * 0.5,
            cs_top - 5.0,
            if pt.cavitation_number.is_finite() {
                format!("{:.2}", pt.cavitation_number)
            } else {
                "inf".to_string()
            }
        );

        // ---- X-axis labels ----
        let _ = write!(
            svg,
            r##"<text x="{center_x:.1}" y="{:.1}" font-size="11" text-anchor="middle" fill="#2c3e50" font-weight="600">{}</text>"##,
            y0 + 18.0,
            super::primitives::escape_xml(&pt.label)
        );
        // Bend radius and velocity sub-labels.
        let _ = write!(
            svg,
            r##"<text x="{center_x:.1}" y="{:.1}" font-size="9" text-anchor="middle" fill="#7f8c8d">R={:.2} mm</text>"##,
            y0 + 32.0,
            pt.bend_radius_mm
        );
        let _ = write!(
            svg,
            r##"<text x="{center_x:.1}" y="{:.1}" font-size="9" text-anchor="middle" fill="#7f8c8d">{:.1} m/s</text>"##,
            y0 + 44.0,
            pt.throat_velocity_m_s
        );
        let _ = write!(
            svg,
            r##"<text x="{center_x:.1}" y="{:.1}" font-size="9" text-anchor="middle" fill="#7f8c8d">K={:.2}</text>"##,
            y0 + 56.0,
            pt.total_loss_coefficient
        );
    }

    // ---- Pressure decay overlay line (green, dashed) ----
    if n >= 2 {
        let mut line_pts = String::new();
        for (i, pt) in points.iter().enumerate() {
            let cx = x0 + group_w * (i as f64 + 0.5);
            let p_frac = ((pt.upstream_pressure_kpa - p_lo) / (p_hi - p_lo)).clamp(0.0, 1.0);
            let cy = y0 - yh * p_frac;
            let _ = write!(line_pts, "{cx:.1},{cy:.1} ");
        }
        let _ = write!(
            svg,
            r##"<polyline fill="none" stroke="#27ae60" stroke-width="2.5" stroke-dasharray="8 5" points="{line_pts}"/>"##
        );
        // Dot markers on the pressure line.
        for (i, pt) in points.iter().enumerate() {
            let cx = x0 + group_w * (i as f64 + 0.5);
            let p_frac = ((pt.upstream_pressure_kpa - p_lo) / (p_hi - p_lo)).clamp(0.0, 1.0);
            let cy = y0 - yh * p_frac;
            let _ = write!(
                svg,
                r##"<circle cx="{cx:.1}" cy="{cy:.1}" r="4" fill="#27ae60"/>"##
            );
            let _ = write!(
                svg,
                r##"<text x="{cx:.1}" y="{:.1}" font-size="9" text-anchor="middle" fill="#1e8449">{} kPa</text>"##,
                cy - 8.0,
                format_pressure(pt.upstream_pressure_kpa)
            );
        }
    }

    // ---- Axis labels ----
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="720" font-size="14" fill="#34495e">Venturi Placement (along serpentine treatment path)</text>"##,
        x0 + xw * 0.22
    );
    let _ = write!(
        svg,
        r##"<text x="16" y="{:.1}" transform="rotate(-90 16,{:.1})" font-size="14" fill="#2c3e50">Dean Number / Cavitation Strength</text>"##,
        y0 - yh * 0.5,
        y0 - yh * 0.5
    );
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="{:.1}" transform="rotate(90 {:.1},{:.1})" font-size="14" fill="#27ae60">Upstream Pressure (kPa)</text>"##,
        rx + 50.0,
        y0 - yh * 0.5,
        rx + 50.0,
        y0 - yh * 0.5
    );

    // ---- Legend ----
    let ly = 86.0;
    let _ = write!(
        svg,
        r##"<rect x="120" y="{ly:.0}" width="12" height="12" fill="#2e86de" fill-opacity="0.85" rx="2"/><text x="138" y="{:.0}" font-size="11" fill="#34495e">Dean number (De)</text>"##,
        ly + 11.0
    );
    let _ = write!(
        svg,
        r##"<rect x="310" y="{ly:.0}" width="12" height="12" fill="#d35400" fill-opacity="0.85" rx="2"/><text x="328" y="{:.0}" font-size="11" fill="#34495e">Cavitation strength (1-σ, label = σ)</text>"##,
        ly + 11.0
    );
    let _ = write!(
        svg,
        r##"<line x1="560" y1="{:.0}" x2="590" y2="{:.0}" stroke="#27ae60" stroke-width="2.5" stroke-dasharray="6 4"/><circle cx="575" cy="{:.0}" r="3" fill="#27ae60"/><text x="598" y="{:.0}" font-size="11" fill="#34495e">Upstream pressure (kPa)</text>"##,
        ly + 6.0,
        ly + 6.0,
        ly + 6.0,
        ly + 11.0
    );

    svg_end(&mut svg);
    svg.0.flush()?;
    Ok(())
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::{write_dean_venturi_placement_figure, DeanVenturiPoint};

    #[test]
    fn dean_venturi_figure_uses_decimal_pressure_labels_for_narrow_ranges() {
        let unique = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock")
            .as_nanos();
        let path = std::env::temp_dir().join(format!("m12-dean-venturi-{unique}.svg"));
        let points = vec![
            DeanVenturiPoint {
                label: "Bend 1 (outer)".to_string(),
                dean_number: 82.0,
                cavitation_number: 0.856,
                throat_velocity_m_s: 16.2,
                total_loss_coefficient: 5307.83,
                upstream_pressure_kpa: 140.00,
                bend_radius_mm: 1.30,
            },
            DeanVenturiPoint {
                label: "Bend 2 (inner)".to_string(),
                dean_number: 84.0,
                cavitation_number: 0.843,
                throat_velocity_m_s: 16.6,
                total_loss_coefficient: 5451.49,
                upstream_pressure_kpa: 140.00,
                bend_radius_mm: 1.30,
            },
            DeanVenturiPoint {
                label: "Bend 3 (outer)".to_string(),
                dean_number: 82.0,
                cavitation_number: 0.858,
                throat_velocity_m_s: 16.2,
                total_loss_coefficient: 5307.83,
                upstream_pressure_kpa: 140.00,
                bend_radius_mm: 1.30,
            },
        ];

        write_dean_venturi_placement_figure(&path, "Dean test", &points)
            .expect("dean venturi figure should render");
        let svg = std::fs::read_to_string(path).expect("rendered svg should exist");

        assert!(svg.contains("139.90"));
        assert!(svg.contains("140.00 kPa"));
    }
}
