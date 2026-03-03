//! SVG drawing renderer — exports engineering drawings to SVG format.
//!
//! Follows ISO 128 line conventions:
//! - Visible edges: stroke-width 0.5mm, solid, black
//! - Hidden edges: stroke-width 0.25mm, dashed, gray
//! - Silhouette edges: stroke-width 0.5mm, solid, black
//! - Dimension lines: stroke-width 0.25mm, solid, blue
//! - Title block: stroke-width 0.35mm, solid, black

use crate::domain::drawing::dimension::DimensionSpec;
use crate::domain::drawing::sheet::DrawingSheet;
use crate::domain::drawing::title_block::TitleBlock;
use crate::infrastructure::drawing::projection::{OrthographicProjector, ProjectedEdges};
use cfd_mesh::IndexedMesh;

/// Render a `DrawingSheet` with associated meshes to an SVG string.
pub fn render_svg(
    sheet: &DrawingSheet,
    meshes: &[&IndexedMesh<f64>],
) -> String {
    let (sheet_w, sheet_h) = sheet.size.dimensions_mm();
    let margin = sheet.border_margin_mm;

    let mut svg = format!(
        r#"<?xml version="1.0" encoding="UTF-8"?>
<svg xmlns="http://www.w3.org/2000/svg" width="{sheet_w}mm" height="{sheet_h}mm" viewBox="0 0 {sheet_w} {sheet_h}">
<style>
  .visible {{ stroke: black; stroke-width: 0.5; fill: none; }}
  .hidden {{ stroke: #888; stroke-width: 0.25; fill: none; stroke-dasharray: 4,2; }}
  .silhouette {{ stroke: black; stroke-width: 0.5; fill: none; }}
  .dimension {{ stroke: #0066cc; stroke-width: 0.25; fill: none; }}
  .dim-text {{ font-family: sans-serif; font-size: 3.5; fill: #0066cc; text-anchor: middle; }}
  .title-text {{ font-family: sans-serif; font-size: 4; fill: black; }}
  .title-border {{ stroke: black; stroke-width: 0.35; fill: none; }}
  .border {{ stroke: black; stroke-width: 0.5; fill: none; }}
</style>
"#
    );

    // Draw border.
    svg.push_str(&format!(
        r#"<rect class="border" x="{margin}" y="{margin}" width="{}" height="{}"/>
"#,
        sheet_w - 2.0 * margin,
        sheet_h - 2.0 * margin,
    ));

    // Draw each projected view.
    for view in &sheet.views {
        if let Some(mesh) = meshes.get(view.mesh_index) {
            let projector = OrthographicProjector::new(
                view.view_type.direction(),
                view.view_type.up(),
            );
            let edges = projector.project_mesh(mesh);
            let cx = view.position_on_sheet_mm[0];
            let cy = sheet_h - view.position_on_sheet_mm[1]; // SVG Y is top-down
            let scale = view.scale;
            render_edges(&mut svg, &edges, cx, cy, scale);

            // View label.
            if let Some(label) = &view.label {
                svg.push_str(&format!(
                    r#"<text class="dim-text" x="{cx}" y="{}">{label}</text>
"#,
                    cy + 15.0,
                ));
            }
        }
    }

    // Draw dimensions.
    for dim in &sheet.dimensions {
        render_dimension(&mut svg, dim, sheet_h);
    }

    // Draw title block.
    render_title_block(&mut svg, &sheet.title_block, sheet_w, sheet_h, margin);

    svg.push_str("</svg>\n");
    svg
}

fn render_edges(
    svg: &mut String,
    edges: &ProjectedEdges,
    cx: f64,
    cy: f64,
    scale: f64,
) {
    for (a, b) in &edges.visible {
        svg.push_str(&format!(
            r#"<line class="visible" x1="{}" y1="{}" x2="{}" y2="{}"/>
"#,
            cx + a.x * scale, cy - a.y * scale,
            cx + b.x * scale, cy - b.y * scale,
        ));
    }
    for (a, b) in &edges.hidden {
        svg.push_str(&format!(
            r#"<line class="hidden" x1="{}" y1="{}" x2="{}" y2="{}"/>
"#,
            cx + a.x * scale, cy - a.y * scale,
            cx + b.x * scale, cy - b.y * scale,
        ));
    }
    for (a, b) in &edges.silhouette {
        svg.push_str(&format!(
            r#"<line class="silhouette" x1="{}" y1="{}" x2="{}" y2="{}"/>
"#,
            cx + a.x * scale, cy - a.y * scale,
            cx + b.x * scale, cy - b.y * scale,
        ));
    }
}

fn render_dimension(svg: &mut String, dim: &DimensionSpec, _sheet_h: f64) {
    if let DimensionSpec::Linear(lin) = dim {
        let value = lin.text_override.as_deref().unwrap_or("");
        let measured = lin.measured_value();
        let label = if value.is_empty() {
            format!("{measured:.2}")
        } else {
            value.to_owned()
        };
        // Simplified: place text at midpoint.
        let mx = (lin.start[0] + lin.end[0]) / 2.0;
        let my = (lin.start[1] + lin.end[1]) / 2.0;
        svg.push_str(&format!(
            r#"<text class="dim-text" x="{mx}" y="{my}">{label}</text>
"#
        ));
    }
    // Angular, Diameter, Radius rendering is similar; omitted for brevity.
}

fn render_title_block(
    svg: &mut String,
    tb: &TitleBlock,
    sheet_w: f64,
    sheet_h: f64,
    margin: f64,
) {
    let tb_width = 180.0_f64.min(sheet_w - 2.0 * margin);
    let tb_height = 40.0;
    let x = sheet_w - margin - tb_width;
    let y = sheet_h - margin - tb_height;

    // Title block border.
    svg.push_str(&format!(
        r#"<rect class="title-border" x="{x}" y="{y}" width="{tb_width}" height="{tb_height}"/>
"#
    ));

    // Horizontal divider.
    let mid_y = y + tb_height / 2.0;
    svg.push_str(&format!(
        r#"<line class="title-border" x1="{x}" y1="{mid_y}" x2="{}" y2="{mid_y}"/>
"#,
        x + tb_width,
    ));

    // Vertical dividers (3 columns).
    let col1 = x + tb_width / 3.0;
    let col2 = x + 2.0 * tb_width / 3.0;
    svg.push_str(&format!(
        r#"<line class="title-border" x1="{col1}" y1="{y}" x2="{col1}" y2="{}"/>
<line class="title-border" x1="{col2}" y1="{y}" x2="{col2}" y2="{}"/>
"#,
        y + tb_height,
        y + tb_height,
    ));

    // Text entries.
    let text_y1 = y + 8.0;
    let text_y2 = mid_y + 8.0;
    svg.push_str(&format!(
        r#"<text class="title-text" x="{}" y="{text_y1}">{}</text>
<text class="title-text" x="{}" y="{text_y1}">{}</text>
<text class="title-text" x="{}" y="{text_y1}">Dwg: {}</text>
<text class="title-text" x="{}" y="{text_y2}">By: {}</text>
<text class="title-text" x="{}" y="{text_y2}">Rev: {} Date: {}</text>
<text class="title-text" x="{}" y="{text_y2}">Scale: {} Sheet {}/{}</text>
"#,
        x + 4.0, tb.company_name,
        col1 + 4.0, tb.drawing_title,
        col2 + 4.0, tb.drawing_number,
        x + 4.0, tb.author,
        col1 + 4.0, tb.revision, tb.date,
        col2 + 4.0, tb.scale, tb.sheet_number, tb.sheet_total,
    ));
}
