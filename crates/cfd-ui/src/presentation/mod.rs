//! Presentation layer — UI views, panels, toolbars, dialogs.

pub mod actions;
pub mod app;
pub mod dialogs;
pub mod menus;
pub mod panels;
pub mod theme;
pub mod toolbar;
pub mod viewport;
#[cfg(feature = "gpui-window")]
pub mod window;
pub mod workspace;
