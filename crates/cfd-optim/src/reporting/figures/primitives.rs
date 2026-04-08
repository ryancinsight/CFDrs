//! Low-level SVG rendering primitives and generic bar-chart builder.

use std::fmt::Write as _;
use std::path::Path;

pub(super) fn write_bar_svg(
    path: &Path,
    title: &str,
    data: &[(&str, f64)],
    y_max: f64,
    x_label: &str,
    y_label: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    if data.is_empty() {
        return super::process::write_placeholder(path, title, "No data available.");
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

    for i in 0..=5 {
        let frac = f64::from(i) / 5.0;
        let val = y_max * frac;
        let ty = y0 - yh * frac;
        let _ = write!(
            svg,
            r##"<line x1="{x0:.2}" y1="{ty:.2}" x2="{:.2}" y2="{ty:.2}" stroke="#ecf0f1" stroke-width="1"/>"##,
            x0 + xw
        );
        let _ = write!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" font-size="12" text-anchor="end" fill="#7f8c8d">{:.1}</text>"##,
            x0 - 8.0,
            ty + 4.0,
            val
        );
    }

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
            r##"<text x="{:.2}" y="{:.2}" font-size="13" text-anchor="middle" fill="#2c3e50">{}</text>"##,
            left + bw * 0.5,
            y0 + 24.0
            ,escape_xml(label)
        );
        let _ = write!(
            svg,
            r##"<text x="{:.2}" y="{:.2}" font-size="12" text-anchor="middle" fill="#2c3e50">{:.3}</text>"##,
            left + bw * 0.5,
            top - 8.0,
            value
        );
    }
    // Axis labels
    let _ = write!(
        svg,
        r##"<text x="{:.1}" y="668" font-size="16" fill="#34495e">{}</text>"##,
        x0 + xw * 0.4,
        escape_xml(x_label)
    );
    let _ = write!(
        svg,
        r##"<text x="18" y="{:.1}" transform="rotate(-90 18,{:.1})" font-size="16" fill="#34495e">{}</text>"##,
        y0 - yh * 0.5,
        y0 - yh * 0.5,
        escape_xml(y_label)
    );
    svg_end(&mut svg);
    std::fs::write(path, svg)?;
    Ok(())
}

pub(super) fn write_bar_svg_owned(
    path: &Path,
    title: &str,
    data: &[(String, f64)],
    y_max: f64,
    x_label: &str,
    y_label: &str,
) -> Result<(), Box<dyn std::error::Error>> {
    let refs: Vec<(&str, f64)> = data.iter().map(|(k, v)| (k.as_str(), *v)).collect();
    write_bar_svg(path, title, &refs, y_max, x_label, y_label)
}

pub(super) fn svg_start(svg: &mut impl std::fmt::Write, width: f64, height: f64) {
    let _ = write!(
        svg,
        r#"<svg xmlns="http://www.w3.org/2000/svg" width="{width:.0}" height="{height:.0}" viewBox="0 0 {width:.0} {height:.0}" preserveAspectRatio="xMidYMin meet" style="max-width:100%;height:auto;">"#
    );
    let _ =
        svg.write_str(r#"<rect x="0" y="0" width="100%" height="100%" fill="white"/>"#);
    let _ = svg.write_str(
        r##"<defs><marker id="arrowhead-flow" markerWidth="10" markerHeight="8" refX="8" refY="4" orient="auto" markerUnits="strokeWidth"><path d="M0,0 L0,8 L10,4 z" fill="#566573"/></marker></defs>"##,
    );
}

pub(super) fn svg_end(svg: &mut impl std::fmt::Write) {
    let _ = svg.write_str("</svg>");
}

pub(super) fn svg_title(svg: &mut impl std::fmt::Write, title: &str) {
    let _ = write!(
        svg,
        r##"<text x="70" y="60" font-size="28" font-weight="600" fill="#1f2d3d">{}</text>"##,
        escape_xml(title)
    );
}

pub(super) fn axis(svg: &mut impl std::fmt::Write, x0: f64, y0: f64, xw: f64, yh: f64) {
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

/// Escape the XML-sensitive subset used by SVG text nodes.
///
/// # Algorithm
/// Performs a single left-to-right pass and expands only `&`, `<`, and `>`.
///
/// # Theorem
/// `escape_xml` preserves the order of all input characters and replaces each
/// escapable character with its XML entity exactly once.
///
/// **Proof sketch**
/// The scan visits each character once in source order. Non-escapable
/// characters are appended unchanged; escapable characters are mapped to one
/// fixed entity string. Because the function never revisits or rewrites output
/// that it has already emitted, it cannot reorder, duplicate, or double-escape
/// characters introduced by earlier replacements.
pub(super) fn escape_xml(text: &str) -> String {
    if !text.as_bytes().iter().any(|byte| matches!(byte, b'&' | b'<' | b'>')) {
        return text.to_string();
    }

    let mut escaped = String::with_capacity(text.len());
    for ch in text.chars() {
        match ch {
            '&' => escaped.push_str("&amp;"),
            '<' => escaped.push_str("&lt;"),
            '>' => escaped.push_str("&gt;"),
            _ => escaped.push(ch),
        }
    }
    escaped
}

/// Estimate average character width for a proportional sans-serif font.
pub(super) fn avg_char_width(font_size: f64) -> f64 {
    font_size * 0.55
}

/// Word-wrap `text` so each line fits within `max_chars` characters.
pub(super) fn wrap_text(text: &str, max_chars: usize) -> Vec<String> {
    let mut lines: Vec<String> = Vec::new();
    let mut current = String::new();
    for word in text.split_whitespace() {
        if current.is_empty() {
            current = word.to_string();
        } else if current.len() + 1 + word.len() <= max_chars {
            current.push(' ');
            current.push_str(word);
        } else {
            lines.push(current);
            current = word.to_string();
        }
    }
    if !current.is_empty() {
        lines.push(current);
    }
    if lines.is_empty() {
        lines.push(String::new());
    }
    lines
}

pub(super) fn process_box(
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
    let title_font = 15.0_f64;
    let body_font = 12.0_f64;
    let pad_x = 14.0_f64;
    let available = (width - 2.0 * pad_x).max(1.0);
    let title_max_chars = (available / avg_char_width(title_font)).floor().max(1.0) as usize;
    let body_max_chars = (available / avg_char_width(body_font)).floor().max(1.0) as usize;

    let _ = write!(
        svg,
        r#"<rect x="{x:.1}" y="{y:.1}" width="{width:.1}" height="{height:.1}" rx="16" ry="16" fill="{fill}" stroke="{stroke}" stroke-width="2.5"/>"#
    );

    let title_lines = wrap_text(title, title_max_chars);
    let title_line_h = title_font * 1.3;
    let mut text_y = y + title_font + 10.0;
    for tl in &title_lines {
        let _ = write!(
            svg,
            r##"<text x="{:.1}" y="{:.1}" font-size="{}" font-weight="600" fill="#1f2d3d">{}</text>"##,
            x + pad_x,
            text_y,
            title_font as u32,
            escape_xml(tl)
        );
        text_y += title_line_h;
    }

    text_y += 4.0;
    let body_line_h = body_font * 1.5;
    for line in lines {
        let wrapped = wrap_text(line, body_max_chars);
        for wl in &wrapped {
            let _ = write!(
                svg,
                r##"<text x="{:.1}" y="{:.1}" font-size="{}" fill="#34495e">{}</text>"##,
                x + pad_x,
                text_y,
                body_font as u32,
                escape_xml(wl)
            );
            text_y += body_line_h;
        }
    }
}

#[cfg(test)]
mod tests {
    use super::{escape_xml, wrap_text};

    #[test]
    fn escape_xml_leaves_plain_text_unchanged() {
        assert_eq!(escape_xml("Alpha 123"), "Alpha 123");
    }

    #[test]
    fn escape_xml_escapes_xml_metacharacters_once() {
        assert_eq!(
            escape_xml("A&B <tag>"),
            "A&amp;B &lt;tag&gt;"
        );
    }

    #[test]
    fn wrap_text_returns_single_empty_line_for_blank_input() {
        assert_eq!(wrap_text("   ", 8), vec![String::new()]);
    }
}

pub(super) fn arrow(svg: &mut String, x1: f64, y1: f64, x2: f64, y2: f64) {
    let _ = write!(
        svg,
        r##"<line x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}" stroke="#566573" stroke-width="3" marker-end="url(#arrowhead-flow)"/>"##
    );
}
