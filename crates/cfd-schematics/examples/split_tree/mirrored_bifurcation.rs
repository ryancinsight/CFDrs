#![allow(deprecated)] // NetworkBlueprint::new() used intentionally; all nodes use NodeSpec::new_at() with explicit positions.

use cfd_schematics::domain::model::{
    ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec,
};
#[path = "../shared/mod.rs"]
mod shared;

fn main() -> Result<(), Box<dyn std::error::Error>> {

    tracing::info!("Mirrored Bifurcation Demo");
    tracing::info!("---------------------------");

    let box_dims = (100.0, 100.0);
    let r = 5.0;

    fn make_arc(
        center: (f64, f64),
        radius: f64,
        start_angle: f64,
        end_angle: f64,
        points: usize,
    ) -> Vec<(f64, f64)> {
        let mut path = Vec::with_capacity(points);
        for i in 0..points {
            let t = i as f64 / (points - 1) as f64;
            let angle = start_angle + t * (end_angle - start_angle);
            path.push((
                center.0 + radius * angle.cos(),
                center.1 + radius * angle.sin(),
            ));
        }
        path
    }

    let nodes = vec![
        NodeSpec::new_at("0", NodeKind::Inlet, (-5.0, 50.0)),
        NodeSpec::new_at("1", NodeKind::Junction, (25.0, 50.0)),
        NodeSpec::new_at("2", NodeKind::Junction, (75.0, 50.0)),
        NodeSpec::new_at("3", NodeKind::Outlet, (105.0, 50.0)),
    ];

    let channel_diameter = 2.0;
    let mut channels = Vec::new();

    let mut add_ch = |from: &str, to: &str, shape: ChannelShape, path: Option<Vec<(f64, f64)>>| {
        let id = format!("{}", channels.len());
        let mut metadata = cfd_schematics::geometry::metadata::MetadataContainer::new();
        metadata.insert(
            cfd_schematics::geometry::metadata::ChannelGeometryMetadata {
                channel_diameter_mm: channel_diameter,
            },
        );

        let mut ch = ChannelSpec::new_pipe_rect(
            &id,
            from,
            to,
            10.0,
            channel_diameter,
            channel_diameter,
            10.0,
            0.0,
        );
        ch.channel_shape = shape;
        if let Some(p) = path {
            ch.path = p;
        }
        ch.metadata = Some(metadata);
        channels.push(ch);
    };

    use std::f64::consts::PI;

    // Straight channels need explicit paths so the renderer doesn't skip them.
    add_ch(
        "0",
        "1",
        ChannelShape::Straight,
        Some(vec![(-5.0, 50.0), (25.0, 50.0)]),
    );

    let mut top_path = Vec::new();
    top_path.push((25.0, 50.0));
    top_path.extend(make_arc((25.0 + r, 75.0 - r), r, PI, PI / 2.0, 16));
    top_path.extend(make_arc((75.0 - r, 75.0 - r), r, PI / 2.0, 0.0, 16));
    top_path.push((75.0, 50.0));
    add_ch(
        "1",
        "2",
        ChannelShape::Serpentine {
            segments: 1,
            bend_radius_m: r / 1000.0,
        },
        Some(top_path),
    );

    let mut bot_path = Vec::new();
    bot_path.push((25.0, 50.0));
    bot_path.extend(make_arc((25.0 + r, 25.0 + r), r, PI, 3.0 * PI / 2.0, 16));
    bot_path.extend(make_arc(
        (75.0 - r, 25.0 + r),
        r,
        3.0 * PI / 2.0,
        2.0 * PI,
        16,
    ));
    bot_path.push((75.0, 50.0));
    add_ch(
        "1",
        "2",
        ChannelShape::Serpentine {
            segments: 1,
            bend_radius_m: r / 1000.0,
        },
        Some(bot_path),
    );

    add_ch(
        "2",
        "3",
        ChannelShape::Straight,
        Some(vec![(75.0, 50.0), (105.0, 50.0)]),
    );

    let box_outline = vec![
        ((0.0, 0.0), (100.0, 0.0)),
        ((100.0, 0.0), (100.0, 100.0)),
        ((100.0, 100.0), (0.0, 100.0)),
        ((0.0, 100.0), (0.0, 0.0)),
    ];

    let mut system = NetworkBlueprint::new("mirrored_bifurcation");
    system.box_dims = box_dims;
    system.box_outline = box_outline;
    system.nodes = nodes;
    system.channels = channels;

    shared::save_example_output(&system, "mirrored_bifurcation");

    Ok(())
}
