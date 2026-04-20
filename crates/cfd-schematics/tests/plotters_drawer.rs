use cfd_schematics::domain::model::NetworkBlueprint;
use cfd_schematics::error::VisualizationError;
use cfd_schematics::visualizations::plotters_backend::{
    PlottersDrawer, PlottersVisualizationEngine,
};
use cfd_schematics::visualizations::traits::{
    Color as CfdColor, GeometricDrawer, RenderConfig, TextStyle, VisualizationEngine,
};
use plotters::prelude::*;

fn synthetic_system() -> NetworkBlueprint {
    NetworkBlueprint {
        name: "test".to_string(),
        box_dims: (10.0, 10.0),
        nodes: Vec::new(),
        channels: Vec::new(),
        box_outline: vec![
            ((0.0, 0.0), (10.0, 0.0)),
            ((10.0, 0.0), (10.0, 10.0)),
            ((10.0, 10.0), (0.0, 10.0)),
            ((0.0, 10.0), (0.0, 0.0)),
        ],
        render_hints: None,
        topology: None,
        lineage: None,
        metadata: None,
        geometry_authored: false,
    }
}

#[test]
fn plotters_drawer_renders_filled_rectangle_and_text() {
    let mut svg = String::new();
    {
        let backend = SVGBackend::with_string(&mut svg, (200, 200));
        let root = backend.into_drawing_area();
        root.fill(&WHITE).expect("background fill must succeed");

        let mut chart = ChartBuilder::on(&root)
            .margin(20)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0.0..10.0, 0.0..10.0)
            .expect("chart must build");

        {
            let mut drawer = PlottersDrawer::new(&root, Some(&mut chart));
            drawer
                .fill_rectangle((1.0, 1.0), (4.0, 5.0), &CfdColor::rgb(10, 20, 30))
                .expect("filled rectangle must render");
            drawer
                .draw_text(
                    (2.0, 3.0),
                    "Sample",
                    &TextStyle::new(CfdColor::rgb(40, 50, 60), 12.0, "sans-serif"),
                )
                .expect("text must render");
        }

        root.present().expect("bitmap must flush");
    }

    let svg_lower = svg.to_ascii_lowercase();
    assert!(svg.contains("Sample"));
    assert!(svg_lower.contains("#0a141e"));
    assert!(svg_lower.contains("#28323c"));
}

#[test]
fn plotters_drawer_requires_chart_context() {
    let mut svg = String::new();
    let result = {
        let backend = SVGBackend::with_string(&mut svg, (200, 200));
        let root = backend.into_drawing_area();
        let mut drawer = PlottersDrawer::new(&root, None);
        drawer.draw_text(
            (2.0, 3.0),
            "Sample",
            &TextStyle::new(CfdColor::rgb(40, 50, 60), 12.0, "sans-serif"),
        )
    };

    match result {
        Err(VisualizationError::CoordinateTransformError { message }) => {
            assert!(message.contains("chart context"));
        }
        other => panic!("expected CoordinateTransformError, got {other:?}"),
    }
}

#[test]
fn plotters_visualization_engine_renders_axes_and_title() {
    let width = 240usize;
    let height = 240usize;
    let mut buffer = vec![255u8; width * height * 3];
    {
        let backend = BitMapBackend::with_buffer(&mut buffer, (width as u32, height as u32));
        let root = backend.into_drawing_area();
        root.fill(&WHITE).expect("background fill must succeed");

        let mut chart = ChartBuilder::on(&root)
            .margin(20)
            .x_label_area_size(30)
            .y_label_area_size(30)
            .build_cartesian_2d(0.0..10.0, 0.0..10.0)
            .expect("chart must build");

        {
            let drawer = PlottersDrawer::new(&root, Some(&mut chart));
            let mut engine = PlottersVisualizationEngine::new(drawer);
            let system = synthetic_system();
            let config = RenderConfig {
                title: "Engine Title".to_string(),
                show_grid: true,
                title_style: TextStyle::new(CfdColor::rgb(15, 30, 45), 18.0, "sans-serif"),
                ..RenderConfig::default()
            };

            engine.add_axes(&system, &config).expect("axes must render");
            engine
                .add_title(&config.title, &config.title_style)
                .expect("title must render");
        }

        root.present().expect("SVG must flush");
    }

    let top_band_non_white = count_non_white_pixels(&buffer, width, height, |_, y| y < 40);
    let edge_band_non_white = count_non_white_pixels(&buffer, width, height, |x, y| {
        x < 50 || x >= width.saturating_sub(50) || y >= height.saturating_sub(50)
    });

    assert!(
        top_band_non_white > 0,
        "title band must contain drawn pixels"
    );
    assert!(
        edge_band_non_white > 0,
        "axes band must contain drawn pixels"
    );
}

fn count_non_white_pixels<F>(buffer: &[u8], width: usize, height: usize, predicate: F) -> usize
where
    F: Fn(usize, usize) -> bool,
{
    assert_eq!(buffer.len(), width * height * 3);
    let mut count = 0usize;
    for (idx, pixel) in buffer.chunks_exact(3).enumerate() {
        let x = idx % width;
        let y = idx / width;
        if predicate(x, y) && pixel[0] != 255 && pixel[1] != 255 && pixel[2] != 255 {
            count += 1;
        }
    }
    count
}
