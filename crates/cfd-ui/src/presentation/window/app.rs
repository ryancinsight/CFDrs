//! Application bootstrap — creates the gpui window and starts the event loop.

use gpui::{px, size, App, AppContext, Application, WindowBounds, WindowOptions};

use super::actions;
use super::workspace::Workspace;

/// Run the gpui application. This function blocks until the window is closed.
pub fn run() -> anyhow::Result<()> {
    Application::new().run(|app: &mut App| {
        actions::register_keybindings(app);

        let bounds = WindowBounds::Windowed(gpui::Bounds {
            origin: gpui::Point::default(),
            size: size(px(1280.0), px(800.0)),
        });

        let options = WindowOptions {
            window_bounds: Some(bounds),
            focus: true,
            ..Default::default()
        };

        app.open_window(options, |_window, cx| cx.new(Workspace::new))
            .expect("failed to open window");
    });

    Ok(())
}
