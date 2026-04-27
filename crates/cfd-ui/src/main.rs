//! cfd-ui binary entry point.

fn main() -> anyhow::Result<()> {
    tracing_subscriber::fmt::init();
    tracing::info!("cfd-ui starting");

    #[cfg(feature = "gpui-window")]
    {
        cfd_ui::presentation::window::app::run()
    }

    #[cfg(not(feature = "gpui-window"))]
    {
        use cfd_ui::presentation::app::CfdApp;
        let app = CfdApp::new()?;
        tracing::info!(
            "workspace initialized with {} meshes",
            app.state.document.mesh_count(),
        );
        println!("cfd-ui v{}", env!("CARGO_PKG_VERSION"));
        println!(
            "Viewport renderer: {}",
            if app.renderer.is_some() {
                "GPU"
            } else {
                "none"
            }
        );
        println!("Run with --features gpui-window for the graphical interface.");
        Ok(())
    }
}
