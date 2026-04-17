//! Drawing setup dialog — sheet size, title block, view placement.

use crate::domain::drawing::sheet::{DrawingSheet, SheetSize};
use crate::domain::drawing::view::{ProjectedView, ViewType};

/// Create a standard three-view engineering drawing (front, top, right + isometric).
#[must_use]
pub fn standard_three_view(
    mesh_index: usize,
    sheet_size: SheetSize,
) -> DrawingSheet {
    let mut sheet = DrawingSheet::new(sheet_size);
    let (w, h) = sheet_size.dimensions_mm();

    // Front view — lower-left quadrant.
    sheet.views.push(ProjectedView {
        view_type: ViewType::Front,
        position_on_sheet_mm: [w * 0.25, h * 0.35],
        scale: 1.0,
        mesh_index,
        label: Some("FRONT".to_owned()),
    });

    // Top view — upper-left quadrant.
    sheet.views.push(ProjectedView {
        view_type: ViewType::Top,
        position_on_sheet_mm: [w * 0.25, h * 0.70],
        scale: 1.0,
        mesh_index,
        label: Some("TOP".to_owned()),
    });

    // Right view — lower-right quadrant.
    sheet.views.push(ProjectedView {
        view_type: ViewType::Right,
        position_on_sheet_mm: [w * 0.65, h * 0.35],
        scale: 1.0,
        mesh_index,
        label: Some("RIGHT".to_owned()),
    });

    // Isometric view — upper-right quadrant.
    sheet.views.push(ProjectedView {
        view_type: ViewType::Isometric,
        position_on_sheet_mm: [w * 0.65, h * 0.70],
        scale: 0.8,
        mesh_index,
        label: Some("ISOMETRIC".to_owned()),
    });

    sheet
}
