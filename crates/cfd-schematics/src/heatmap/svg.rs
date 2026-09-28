use super::data::CandidateZoneData;
use super::palette::cancer_cav_color;
use super::plate::{
    MARGIN_X, MARGIN_Y, PITCH, PLATE_H_MM, PLATE_W_MM, WELL_A1_X, WELL_A1_Y, WELL_R,
    ZONE_CENTER_SPAN_MM, ZONE_COL_START, ZONE_ENVELOPE_MM, ZONE_ROW_START, ZONE_WELLS, mm_to_px,
    px_x, px_y,
};
use std::fmt::Write;

/// Build the complete SVG string.
pub(super) fn build_svg(top_candidates: &[CandidateZoneData]) -> String {
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
                cand.label
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
        r##"  <text x="{:.0}" y="{:.0}"
        font-family="Arial,sans-serif" font-size="9" font-weight="bold" fill="#333">
    Cancer Cavitation Score:</text>
  <rect x="{:.0}" y="{:.0}" width="180" height="11" fill="url(#cavGrad)" rx="2"/>
  <text x="{:.0}" y="{:.0}" font-family="Arial,sans-serif" font-size="8" fill="#555">0.0 (low)</text>
  <text x="{:.0}" y="{:.0}" font-family="Arial,sans-serif" font-size="8"
        fill="#555" text-anchor="end">1.0 (high)</text>
"##,
        leg_x,
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
