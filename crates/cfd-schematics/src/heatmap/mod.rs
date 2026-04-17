//! 96-well plate treatment-zone visualization.
//!
//! Generates an SVG showing a standard SBS 96-well plate (127.76 × 85.47 mm)
//! with the 6×6 SDT treatment zone highlighted and top candidate metrics
//! overlaid as colour bars inside the treatment zone.
//!
//! # SBS plate geometry (ANSI/SLAS standard)
//! - Overall: 127.76 × 85.47 mm
//! - Well pitch: 9.00 mm
//! - First well (A1) centre: x = 14.38 mm, y = 11.24 mm
//! - Grid: 12 columns × 8 rows = 96 wells
//!
//! # Treatment zone
//! - 45 × 45 mm centre-to-centre span centred at (63.88, 42.74) mm
//! - Spans wells B4–G9 (rows B–G, columns 4–9 in 1-indexed notation)
//! - The dashed overlay expands half a pitch beyond the outer wells so the
//!   full 6 × 6 well envelope is visible
//!
//! # Output
//! Pure SVG XML written to a file — no external graphical dependencies.

use std::fmt::Write;
use std::path::Path;

// ── Public types ──────────────────────────────────────────────────────────────

/// Data for a single candidate to overlay on the well-plate diagram.
#[derive(Debug, Clone)]
pub struct CandidateZoneData {
    /// Short display label, e.g. `"CCT-lv3"`.
    pub label: String,
    /// Cancer-targeted cavitation score [0–1].
    pub cancer_cav: f64,
    /// Lysis risk index [0–1].
    pub lysis_risk: f64,
    /// Therapy channel / separation efficiency metric [0–1].
    pub therapy_frac: f64,
}

// ── SVG layout constants ──────────────────────────────────────────────────────

/// Scale factor: 1 mm → 6 px.
const SCALE: f64 = 6.0;

/// Left/right margin [px].
const MARGIN_X: f64 = 62.0;
/// Top margin [px].
const MARGIN_Y: f64 = 70.0;

/// SBS plate width [mm].
const PLATE_W_MM: f64 = 127.76;
/// SBS plate height [mm].
const PLATE_H_MM: f64 = 85.47;

/// Well pitch [mm].
const PITCH: f64 = 9.0;
/// First well (A1) centre X [mm].
const WELL_A1_X: f64 = 14.38;
/// First well (A1) centre Y [mm].
const WELL_A1_Y: f64 = 11.24;
/// Well drawing radius [mm].
const WELL_R: f64 = 3.5;

/// Treatment-zone centre-to-centre span [mm] across 6 wells (5 pitches).
const ZONE_CENTER_SPAN_MM: f64 = 45.0;
/// Treatment zone first column index (0-based) — column 4 in 1-indexed (3 in 0-indexed).
const ZONE_COL_START: usize = 3;
/// Treatment zone first row index (0-based) — row B in A–H labelling (1 in 0-indexed).
const ZONE_ROW_START: usize = 1;
/// Number of wells in each axis of the treatment zone.
const ZONE_WELLS: usize = 6;
/// Full highlighted treatment-zone envelope [mm] including a half-pitch border
/// around the first and last well centres.
const ZONE_ENVELOPE_MM: f64 = PITCH * ZONE_WELLS as f64;

// ── Public function ───────────────────────────────────────────────────────────

/// Write an SVG 96-well plate diagram to `output_path`.
///
/// Up to 5 `top_candidates` are overlaid as colour bars inside the treatment
/// zone.  Bar colour encodes `cancer_cav` (yellow → red gradient).
///
/// # Errors
///
/// Returns a boxed error if the output file cannot be written.
pub fn write_well_plate_diagram_svg(
    top_candidates: &[CandidateZoneData],
    output_path: &Path,
) -> Result<(), Box<dyn std::error::Error>> {
    let svg = build_svg(top_candidates);
    std::fs::write(output_path, svg.as_bytes())?;
    Ok(())
}

// ── Internal helpers ──────────────────────────────────────────────────────────

#[inline]
fn mm_to_px(mm: f64) -> f64 {
    mm * SCALE
}

#[inline]
fn px_x(mm: f64) -> f64 {
    MARGIN_X + mm_to_px(mm)
}

#[inline]
fn px_y(mm: f64) -> f64 {
    MARGIN_Y + mm_to_px(mm)
}

/// Interpolate hex colour: yellow (#FFD700) at 0.0, red (#CC0000) at 1.0.
fn cancer_cav_color(v: f64) -> String {
    let t = v.clamp(0.0, 1.0);
    let r = ((0.8 + 0.2 * (1.0 - t)) * 255.0) as u8;
    let g = ((1.0 - t) * 215.0) as u8;
    format!("#{r:02X}{g:02X}00")
}

/// Build the complete SVG string.
fn build_svg(top_candidates: &[CandidateZoneData]) -> String {
    let plate_px_w = mm_to_px(PLATE_W_MM);
    let plate_px_h = mm_to_px(PLATE_H_MM);
    let svg_w = MARGIN_X * 2.0 + plate_px_w;
    let svg_h = MARGIN_Y * 2.0 + plate_px_h + 110.0; // extra for legend row
    let aspect_ratio = if svg_h > 0.0 { svg_w / svg_h } else { 1.0 };

    let mut s = String::with_capacity(20_000);

    // ── XML header + SVG open ─────────────────────────────────────────────────
    let _ = write!(
        s,
        r##"<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg"
         width="{svg_w:.0}" height="{svg_h:.0}" viewBox="0 0 {svg_w:.0} {svg_h:.0}"
                 preserveAspectRatio="xMidYMin meet" style="width:min(100%, 100vw, calc(100vh * {aspect_ratio:.6}));height:auto;display:block;margin:0 auto;">
  <defs>
    <linearGradient id="cavGrad" x1="0%" y1="0%" x2="100%" y2="0%">
      <stop offset="0%"   style="stop-color:#FFD700"/>
      <stop offset="100%" style="stop-color:#CC0000"/>
    </linearGradient>
  </defs>
  <rect width="{svg_w:.0}" height="{svg_h:.0}" fill="#F5F5F5"/>
"##
    );

    // ── Title ─────────────────────────────────────────────────────────────────
    let _ = write!(
        s,
        r##"  <text x="{cx:.0}" y="30"
        font-family="Arial,sans-serif" font-size="16" font-weight="700"
        fill="#222" text-anchor="middle">SDT Millifluidic Device</text>
  <text x="{cx:.0}" y="49"
        font-family="Arial,sans-serif" font-size="13" font-weight="600"
        fill="#44505a" text-anchor="middle">96-Well Plate Treatment Zone</text>
"##,
        cx = svg_w / 2.0
    );

    // ── Plate background rectangle ────────────────────────────────────────────
    let _ = write!(
        s,
        r##"  <rect x="{MARGIN_X:.1}" y="{MARGIN_Y:.1}" width="{plate_px_w:.1}" height="{plate_px_h:.1}"
        rx="7" ry="7" fill="white" stroke="#AAAAAA" stroke-width="1.5"/>
"##
    );

    // ── Treatment zone highlight ──────────────────────────────────────────────
    // Zone upper-left corner: half a pitch before the first zone well centre
    let zone_mm_x = WELL_A1_X + ZONE_COL_START as f64 * PITCH - PITCH / 2.0;
    let zone_mm_y = WELL_A1_Y + ZONE_ROW_START as f64 * PITCH - PITCH / 2.0;
    let _ = write!(
        s,
        r##"  <rect x="{:.1}" y="{:.1}" width="{:.1}" height="{:.1}"
        rx="4" ry="4" fill="#DCF0FF" stroke="#3377BB"
        stroke-width="2" stroke-dasharray="7,3" opacity="0.6"/>
"##,
        px_x(zone_mm_x),
        px_y(zone_mm_y),
        mm_to_px(ZONE_ENVELOPE_MM),
        mm_to_px(ZONE_ENVELOPE_MM)
    );

    // ── Wells ─────────────────────────────────────────────────────────────────
    let row_labels: [char; 8] = ['A', 'B', 'C', 'D', 'E', 'F', 'G', 'H'];
    #[allow(clippy::needless_range_loop)]
    for row in 0..8_usize {
        for col in 0..12_usize {
            let cx_mm = WELL_A1_X + col as f64 * PITCH;
            let cy_mm = WELL_A1_Y + row as f64 * PITCH;
            let in_zone = (ZONE_COL_START..ZONE_COL_START + ZONE_WELLS).contains(&col)
                && (ZONE_ROW_START..ZONE_ROW_START + ZONE_WELLS).contains(&row);

            let (fill, stroke_col) = if in_zone {
                ("#A8CCF0", "#3377BB")
            } else {
                ("#D5D5D5", "#999999")
            };
            let sw = if in_zone { 1.0_f64 } else { 0.6 };

            let _ = write!(
                s,
                r#"  <circle cx="{:.1}" cy="{:.1}" r="{:.1}"
          fill="{fill}" stroke="{stroke_col}" stroke-width="{sw}"/>
"#,
                px_x(cx_mm),
                px_y(cy_mm),
                mm_to_px(WELL_R)
            );

            // Tiny well label
            let _ = write!(
                s,
                r##"  <text x="{:.1}" y="{:.1}"
          font-family="Arial,sans-serif" font-size="6.5"
          fill="#777" text-anchor="middle" dominant-baseline="central">{}{}</text>
"##,
                px_x(cx_mm),
                px_y(cy_mm),
                row_labels[row],
                col + 1
            );
        }
    }

    // ── Column header numbers (1–12) ──────────────────────────────────────────
    for col in 0..12_usize {
        let cx_mm = WELL_A1_X + col as f64 * PITCH;
        let _ = write!(
            s,
            r##"  <text x="{:.1}" y="{:.0}"
          font-family="Arial,sans-serif" font-size="9"
          fill="#555" text-anchor="middle">{}</text>
"##,
            px_x(cx_mm),
                        MARGIN_Y - 14.0,
            col + 1
        );
    }

    // ── Row header letters (A–H) ──────────────────────────────────────────────
    for (row, &ch) in row_labels.iter().enumerate() {
        let cy_mm = WELL_A1_Y + row as f64 * PITCH;
        let _ = write!(
            s,
            r##"  <text x="{:.0}" y="{:.1}"
          font-family="Arial,sans-serif" font-size="9"
          fill="#555" text-anchor="end" dominant-baseline="central">{ch}</text>
"##,
            MARGIN_X - 6.0,
            px_y(cy_mm)
        );
    }

    // ── Candidate metric bars ─────────────────────────────────────────────────
    let n_cands = top_candidates.len().min(5);
    if n_cands > 0 {
        let zone_px_w = mm_to_px(ZONE_ENVELOPE_MM);
        let zone_px_h = mm_to_px(ZONE_ENVELOPE_MM);
        let bar_gap = 4.0;
        let bar_w = (zone_px_w - bar_gap * (n_cands as f64 + 1.0)) / n_cands as f64;
        let bar_h = zone_px_h - 28.0;
        let bar_top = px_y(zone_mm_y) + 14.0;
        let bar_left0 = px_x(zone_mm_x) + bar_gap;

        for (i, cand) in top_candidates.iter().take(5).enumerate() {
            let bx = bar_left0 + i as f64 * (bar_w + bar_gap);
            let fill = cancer_cav_color(cand.cancer_cav);
            let opacity = (0.55 + 0.40 * cand.cancer_cav).min(0.95);

            // Bar body
            let _ = write!(
                s,
                r#"  <rect x="{bx:.1}" y="{bar_top:.1}" width="{bar_w:.1}" height="{bar_h:.1}"
          fill="{fill}" rx="3" ry="3" opacity="{opacity:.2}"/>
"#
            );

            // Label (top of bar)
            let lbl_x = bx + bar_w / 2.0;
            let _ = write!(
                s,
                "  <text x=\"{:.1}\" y=\"{:.1}\"\n          font-family=\"Arial,sans-serif\" font-size=\"7.5\"\n          fill=\"white\" font-weight=\"bold\" text-anchor=\"middle\">{}</text>\n",
                lbl_x,
                bar_top + 10.0,
                &cand.label
            );

            // Cancer cav score (bottom of bar)
            let _ = write!(
                s,
                "  <text x=\"{:.1}\" y=\"{:.1}\"\n          font-family=\"Arial,sans-serif\" font-size=\"7\"\n          fill=\"white\" text-anchor=\"middle\">{:.2}</text>\n",
                lbl_x,
                bar_top + bar_h - 8.0,
                cand.cancer_cav
            );
        }
    }

    // ── Legend ────────────────────────────────────────────────────────────────
    let leg_y = MARGIN_Y + plate_px_h + 18.0;
    let leg_x = MARGIN_X + 10.0;

    let _ = write!(
        s,
        r##"  <text x="{leg_x:.0}" y="{:.0}"
        font-family="Arial,sans-serif" font-size="9" font-weight="bold" fill="#333">
    Cancer Cavitation Score:</text>
  <rect x="{:.0}" y="{:.0}" width="180" height="11" fill="url(#cavGrad)" rx="2"/>
  <text x="{:.0}" y="{:.0}" font-family="Arial,sans-serif" font-size="8" fill="#555">0.0 (low)</text>
  <text x="{:.0}" y="{:.0}" font-family="Arial,sans-serif" font-size="8"
        fill="#555" text-anchor="end">1.0 (high)</text>
"##,
        leg_y + 12.0,
        leg_x + 180.0,
        leg_y + 18.0,
        leg_x + 180.0,
        leg_y + 38.0,
        leg_x + 365.0,
        leg_y + 38.0
    );

    // Zone legend swatch
    let _ = write!(
        s,
        r##"  <rect x="{:.0}" y="{:.0}" width="13" height="13" fill="#A8CCF0"
        stroke="#3377BB" stroke-width="1.2" stroke-dasharray="4,2" rx="2"/>
  <text x="{:.0}" y="{:.0}" font-family="Arial,sans-serif" font-size="8.5" fill="#444">
    Treatment zone (6×6 wells · {:.0}×{:.0} mm centre span)</text>
"##,
        leg_x + 400.0,
        leg_y + 18.0,
        leg_x + 418.0,
        leg_y + 29.0,
        ZONE_CENTER_SPAN_MM,
        ZONE_CENTER_SPAN_MM,
    );

    // ── Close SVG ─────────────────────────────────────────────────────────────
    s.push_str("</svg>\n");
    s
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    fn sample_candidates() -> Vec<CandidateZoneData> {
        vec![
            CandidateZoneData {
                label: "CCT-lv3".into(),
                cancer_cav: 0.78,
                lysis_risk: 0.001,
                therapy_frac: 0.88,
            },
            CandidateZoneData {
                label: "TBV-cf55".into(),
                cancer_cav: 0.65,
                lysis_risk: 0.002,
                therapy_frac: 0.72,
            },
        ]
    }

    #[test]
    fn svg_starts_with_xml_declaration() {
        let svg = build_svg(&sample_candidates());
        assert!(
            svg.starts_with("<?xml"),
            "expected XML declaration at start"
        );
    }

    #[test]
    fn svg_contains_svg_element() {
        let svg = build_svg(&sample_candidates());
        assert!(svg.contains("<svg "), "expected <svg element");
        assert!(svg.contains("</svg>"), "expected closing </svg>");
        assert!(svg.contains("preserveAspectRatio=\"xMidYMin meet\""));
        assert!(svg.contains("width:min(100%, 100vw, calc(100vh * "));
        assert!(svg.contains("height:auto;display:block;margin:0 auto;"));
    }

    #[test]
    fn svg_uses_split_title_for_plate_heading() {
        let svg = build_svg(&sample_candidates());
        assert!(svg.contains("SDT Millifluidic Device"));
        assert!(svg.contains("96-Well Plate Treatment Zone"));
    }

    #[test]
    fn svg_contains_treatment_zone_rect() {
        let svg = build_svg(&sample_candidates());
        // The dashed treatment zone rect must be present
        assert!(
            svg.contains("stroke-dasharray"),
            "expected dashed border for treatment zone"
        );
    }

    #[test]
    fn svg_treatment_zone_rect_covers_full_six_by_six_envelope() {
        let svg = build_svg(&sample_candidates());
        assert!(
            svg.contains(r#"<rect x="283.3" y="164.4" width="324.0" height="324.0""#),
            "expected treatment zone rect to span the full 6x6 well envelope"
        );
    }

    #[test]
    fn svg_contains_candidate_labels() {
        let svg = build_svg(&sample_candidates());
        assert!(
            svg.contains("CCT-lv3"),
            "expected first candidate label in SVG"
        );
        assert!(
            svg.contains("TBV-cf55"),
            "expected second candidate label in SVG"
        );
    }

    #[test]
    fn svg_empty_candidates_still_valid() {
        let svg = build_svg(&[]);
        assert!(svg.contains("<svg "));
        assert!(svg.contains("</svg>"));
        assert!(!svg.contains("CCT")); // no candidate content
    }

    #[test]
    fn cancer_cav_color_extremes() {
        let low = cancer_cav_color(0.0);
        let high = cancer_cav_color(1.0);
        // Low should be yellowish (#FFD700 area), high should be reddish (#CC0000)
        assert!(low.starts_with('#'), "expected hex color");
        assert!(high.starts_with('#'), "expected hex color");
        assert_ne!(low, high, "colors at 0 and 1 should differ");
    }

    #[test]
    fn write_to_file_creates_file() {
        let tmp = std::env::temp_dir().join("test_well_plate.svg");
        let result = write_well_plate_diagram_svg(&sample_candidates(), &tmp);
        assert!(result.is_ok(), "write should succeed: {:?}", result.err());
        assert!(tmp.exists(), "output file should exist");
        let _ = std::fs::remove_file(&tmp); // cleanup
    }
}
