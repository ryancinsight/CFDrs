//! Panel elements — model tree, property editor, and console output.

use gpui::{
    div, px, InteractiveElement, IntoElement, ParentElement, SharedString,
    StatefulInteractiveElement, Styled,
};

use super::workspace::gpui_rgba;
use crate::presentation::panels::console_panel::{ConsoleState, LogLevel};
use crate::presentation::panels::model_tree_panel;
use crate::presentation::panels::property_editor::{PropertySheet, PropertyValue};
use crate::presentation::theme::ThemeColors;
use crate::presentation::workspace::WorkspaceState;
use crate::domain::document::model_tree::ModelTreeNode;
use crate::domain::scene::graph::SceneEntity;

// ---------------------------------------------------------------------------
// Model tree
// ---------------------------------------------------------------------------

/// Render the model tree panel (left sidebar).
pub fn render_model_tree(state: &WorkspaceState, theme: &ThemeColors) -> impl IntoElement {
    let tree = model_tree_panel::build_model_tree(&state.document.scene);

    let mut panel = div()
        .id("model-tree")
        .flex()
        .flex_col()
        .w(px(200.0))
        .min_w(px(200.0))
        .bg(gpui_rgba(theme.panel_bg))
        .border_r_1()
        .border_color(gpui_rgba(theme.border))
        .overflow_y_scroll()
        .child(panel_header("Model Tree", theme));

    for node in &tree {
        panel = panel.child(render_tree_node(node, theme, 0));
    }
    panel
}

/// Render a single tree node, recursing into children.
fn render_tree_node(
    node: &ModelTreeNode,
    theme: &ThemeColors,
    depth: usize,
) -> impl IntoElement {
    let indent = (depth as f32) * 16.0;
    let label = node.label.clone();
    let id_str = format!("tree-{}-{}", node.scene_index, label);

    let mut row = div()
        .id(SharedString::from(id_str))
        .flex()
        .flex_col()
        .child(
            div()
                .pl(px(indent + 8.0))
                .py(px(2.0))
                .text_xs()
                .text_color(gpui_rgba(theme.text_primary))
                .hover(|s| s.bg(gpui::rgba(0xffffff10)))
                .child(label),
        );

    for child in &node.children {
        row = row.child(render_tree_node(child, theme, depth + 1));
    }
    row
}

// ---------------------------------------------------------------------------
// Properties
// ---------------------------------------------------------------------------

/// Render the properties panel (right sidebar).
pub fn render_properties(state: &WorkspaceState, theme: &ThemeColors) -> impl IntoElement {
    let sheet = if let Some(selected) = state.selection.single() {
        if let Some(node) = state.document.scene.node(selected) {
            match &node.entity {
                SceneEntity::Mesh(handle) => {
                    if let Some(mesh) = state.document.mesh(*handle) {
                        PropertySheet::for_mesh(
                            &node.name,
                            mesh.vertices.len(),
                            mesh.faces.len(),
                        )
                    } else {
                        PropertySheet::empty()
                    }
                }
                _ => PropertySheet::empty(),
            }
        } else {
            PropertySheet::empty()
        }
    } else {
        PropertySheet::empty()
    };

    let mut panel = div()
        .id("properties")
        .flex()
        .flex_col()
        .w(px(220.0))
        .min_w(px(220.0))
        .bg(gpui_rgba(theme.panel_bg))
        .border_l_1()
        .border_color(gpui_rgba(theme.border))
        .child(panel_header("Properties", theme));

    if !sheet.object_name.is_empty() {
        panel = panel.child(
            div()
                .px_2()
                .py_1()
                .text_sm()
                .font_weight(gpui::FontWeight::BOLD)
                .text_color(gpui_rgba(theme.text_primary))
                .child(sheet.object_name.clone()),
        );
    }

    for prop in &sheet.properties {
        panel = panel.child(
            div()
                .flex()
                .flex_row()
                .justify_between()
                .px_2()
                .py(px(2.0))
                .child(
                    div()
                        .text_xs()
                        .text_color(gpui_rgba(theme.text_secondary))
                        .child(prop.name.clone()),
                )
                .child(
                    div()
                        .text_xs()
                        .text_color(gpui_rgba(theme.text_primary))
                        .child(format_value(&prop.value)),
                ),
        );
    }
    panel
}

/// Format a property value for display.
fn format_value(value: &PropertyValue) -> String {
    match value {
        PropertyValue::Text(s) => s.clone(),
        PropertyValue::Float(f) => format!("{f:.4}"),
        PropertyValue::Integer(i) => i.to_string(),
        PropertyValue::Boolean(b) => if *b { "Yes" } else { "No" }.to_owned(),
    }
}

// ---------------------------------------------------------------------------
// Console
// ---------------------------------------------------------------------------

/// Render the console panel (bottom bar).
pub fn render_console(console: &ConsoleState, theme: &ThemeColors) -> impl IntoElement {
    let mut panel = div()
        .id("console")
        .flex()
        .flex_col()
        .h(px(120.0))
        .min_h(px(120.0))
        .bg(gpui_rgba(theme.panel_bg))
        .border_t_1()
        .border_color(gpui_rgba(theme.border))
        .overflow_y_scroll()
        .child(panel_header("Console", theme));

    // Show most recent entries at the top (reversed).
    for entry in console.entries().iter().rev().take(50) {
        let color = match entry.level {
            LogLevel::Info => theme.text_secondary,
            LogLevel::Warning => theme.accent,
            LogLevel::Error => theme.error,
        };
        panel = panel.child(
            div()
                .px_2()
                .py(px(1.0))
                .text_xs()
                .text_color(gpui_rgba(color))
                .child(entry.message.clone()),
        );
    }
    panel
}

// ---------------------------------------------------------------------------
// Shared helpers
// ---------------------------------------------------------------------------

/// Render a panel header row.
fn panel_header(title: &str, theme: &ThemeColors) -> impl IntoElement {
    div()
        .px_2()
        .py_1()
        .text_xs()
        .font_weight(gpui::FontWeight::BOLD)
        .text_color(gpui_rgba(theme.text_secondary))
        .border_b_1()
        .border_color(gpui_rgba(theme.border))
        .child(title.to_owned())
}
