//! Toolbar element — horizontal bar with action buttons.

use gpui::{
    div, px, Action, ClickEvent, IntoElement, ParentElement, SharedString,
    StatefulInteractiveElement, Styled, InteractiveElement, Window, App,
};

use super::actions::{
    CreateCapsule, CreateCone, CreateCube, CreateCylinder, CreateElbow,
    CreateEllipsoid, CreateFrustum, CreateGeodesicSphere, CreatePipe,
    CreateRoundedCube, CreateSphere, CreateTorus,
    CsgDifference, CsgIntersection, CsgUnion, ExportDxf, ExportDrawing, ExportGlb,
    ExportObj, ExportOpenFoam, ExportPly, ExportStl, FitView, ImportMesh,
    NewProject, OpenProject, Redo, SaveProject, ToggleWireframe, Undo,
};
use super::workspace::gpui_rgba;
use crate::presentation::theme::ThemeColors;

/// Build the main toolbar element.
pub fn render_toolbar(theme: &ThemeColors) -> impl IntoElement {
    div()
        .id("toolbar")
        .flex()
        .flex_row()
        .items_center()
        .gap_1()
        .h(px(36.0))
        .px_2()
        .bg(gpui_rgba(theme.toolbar_bg))
        .border_b_1()
        .border_color(gpui_rgba(theme.border))
        // File group
        .child(toolbar_button("New", NewProject))
        .child(toolbar_button("Open", OpenProject))
        .child(toolbar_button("Save", SaveProject))
        .child(toolbar_separator(theme))
        // Import/Export group
        .child(toolbar_button("Import", ImportMesh))
        .child(toolbar_button("STL", ExportStl))
        .child(toolbar_button("OBJ", ExportObj))
        .child(toolbar_button("PLY", ExportPly))
        .child(toolbar_button("GLB", ExportGlb))
        .child(toolbar_button("DXF", ExportDxf))
        .child(toolbar_button("SVG", ExportDrawing))
        .child(toolbar_button("OF", ExportOpenFoam))
        .child(toolbar_separator(theme))
        // Edit group
        .child(toolbar_button("Undo", Undo))
        .child(toolbar_button("Redo", Redo))
        .child(toolbar_separator(theme))
        // Mesh group
        .child(toolbar_button("Cube", CreateCube))
        .child(toolbar_button("Cyl", CreateCylinder))
        .child(toolbar_button("Sph", CreateSphere))
        .child(toolbar_button("Cone", CreateCone))
        .child(toolbar_button("Torus", CreateTorus))
        .child(toolbar_button("Pipe", CreatePipe))
        .child(toolbar_separator(theme))
        // Extended primitives group
        .child(toolbar_button("Ellips", CreateEllipsoid))
        .child(toolbar_button("Caps", CreateCapsule))
        .child(toolbar_button("Frust", CreateFrustum))
        .child(toolbar_button("Geo", CreateGeodesicSphere))
        .child(toolbar_button("RCube", CreateRoundedCube))
        .child(toolbar_button("Elbow", CreateElbow))
        .child(toolbar_separator(theme))
        // CSG group
        .child(toolbar_button("Union", CsgUnion))
        .child(toolbar_button("Isect", CsgIntersection))
        .child(toolbar_button("Diff", CsgDifference))
        .child(toolbar_separator(theme))
        // View group
        .child(toolbar_button("Fit", FitView))
        .child(toolbar_button("Wire", ToggleWireframe))
}

/// A single toolbar button that dispatches a gpui action on click.
fn toolbar_button(label: &str, action: impl Action) -> impl IntoElement {
    let action_clone = action.boxed_clone();
    div()
        .id(SharedString::from(label.to_owned()))
        .px_2()
        .py_1()
        .rounded_sm()
        .text_sm()
        .text_color(gpui_rgba([220, 220, 230, 255]))
        .cursor_pointer()
        .hover(|s| s.bg(gpui::rgba(0xffffff20)))
        .active(|s| s.bg(gpui::rgba(0xffffff40)))
        .on_click(move |_ev: &ClickEvent, window: &mut Window, cx: &mut App| {
            window.dispatch_action(action_clone.boxed_clone(), cx);
        })
        .child(label.to_owned())
}

/// A thin vertical separator between toolbar groups.
fn toolbar_separator(theme: &ThemeColors) -> impl IntoElement {
    div()
        .w(px(1.0))
        .h(px(20.0))
        .mx_1()
        .bg(gpui_rgba(theme.border))
}
