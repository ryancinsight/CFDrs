use cfd_schematics::geometry::{Channel, ChannelSystem, ChannelType, Node};
use cfd_schematics::visualizations::schematic::plot_geometry;
use std::fs;
use std::path::PathBuf;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    let manifest_dir = PathBuf::from(env!("CARGO_MANIFEST_DIR"));
    let out_dir = manifest_dir.join("outputs/mirrored_bifurcation");
    fs::create_dir_all(&out_dir)?;

    println!("Mirrored Bifurcation Demo");
    println!("---------------------------");

    // 100x100 bounding box
    let box_dims = (100.0, 100.0);

    // split at 25% length -> x = 25.0
    // merge at 75% length -> x = 75.0
    // parallel at 25% and 75% width -> y = 25.0, y = 75.0

    // Add nodes for the curved path
    // We will use 5.0mm radius for the corners.
    let r = 5.0;
    
    // Instead of doing manual Arc geometry points which can be tedious, 
    // let's just use the `geometry::generator` to create a `SmoothStraight` 
    // but the generator is for standard bifurcation, so let's stick to manuall building 
    // the layout with SmoothStraight since my `SmoothStraightChannelStrategy` can generate the points.
    // Actually, `SmoothStraight` generates points assuming a wavy transition.
    // What we really want is an `Arc` for the corners. 
    
    // An Arc path is just a Sequence of Point2D. 
    // Let's create a helper to generate a 90-degree quarter-circle arc.
    fn make_arc(center: (f64, f64), radius: f64, start_angle: f64, end_angle: f64, points: usize) -> Vec<(f64, f64)> {
        let mut path = Vec::with_capacity(points);
        for i in 0..points {
            let t = i as f64 / (points - 1) as f64;
            let angle = start_angle + t * (end_angle - start_angle);
            path.push((center.0 + radius * angle.cos(), center.1 + radius * angle.sin()));
        }
        path
    }

    // Replace the simple 8 node design with a design that has straight segments + arc segments.
    // Split at x=25, merge at x=75.
    // Let's place the curve centers at:
    // Top-Left corner: center (25+r, 75-r), start angle 180 (Left), end angle 90 (Up)
    // Top-Right corner: center (75-r, 75-r), start angle 90 (Up), end angle 0 (Right)
    // Bottom-Left corner: center (25+r, 25+r), start angle 180 (Left), end angle 270 (Down)
    // Bottom-Right corner: center (75-r, 25+r), start angle 270 (Down), end angle 0 (Right)

    let nodes = vec![
        Node { id: 0, point: (-5.0, 50.0), metadata: None },    // Inlet
        Node { id: 1, point: (26.0, 50.0), metadata: None },   // Split Jct (Overshoot for CSG intersections)
        Node { id: 2, point: (74.0, 50.0), metadata: None },   // Merge Jct (Overshoot for CSG intersections)
        Node { id: 3, point: (105.0, 50.0), metadata: None },  // Outlet
    ];

    let channel_diameter = 2.0;

    let mut channels = Vec::new();
    let mut add_ch = |from: usize, to: usize, channel_type: ChannelType| {
        let id = channels.len();
        let mut metadata = cfd_schematics::geometry::metadata::MetadataContainer::new();
        metadata.insert(cfd_schematics::geometry::metadata::ChannelGeometryMetadata {
            channel_diameter_mm: channel_diameter,
        });

        channels.push(Channel {
            id,
            from_node: from,
            to_node: to,
            width: channel_diameter,
            height: channel_diameter,
            channel_type,
            metadata: Some(metadata),
        });
    };

    use std::f64::consts::PI;

    // Inlet to Split
    add_ch(0, 1, ChannelType::Straight);
    
    // Top Branch (Split to Merge)
    let mut top_path = Vec::new();
    // Deep intersect inside the void of ch0/ch3
    top_path.push((25.0, 50.0)); 
    // arc from (25.0, 70.0) to (30.0, 75.0)
    top_path.extend(make_arc((25.0 + r, 75.0 - r), r, PI, PI/2.0, 16));
    // arc from (70.0, 75.0) to (75.0, 70.0)
    top_path.extend(make_arc((75.0 - r, 75.0 - r), r, PI/2.0, 0.0, 16));
    top_path.push((75.0, 50.0));
    add_ch(1, 2, ChannelType::Serpentine { path: top_path });
    
    // Bottom Branch (Split to Merge)
    let mut bot_path = Vec::new();
    // Deep intersect inside the void of ch0/ch3
    bot_path.push((25.0, 50.0));
    // arc from (25.0, 30.0) to (30.0, 25.0)
    bot_path.extend(make_arc((25.0 + r, 25.0 + r), r, PI, 3.0*PI/2.0, 16));
    // arc from (70.0, 25.0) to (75.0, 30.0)
    bot_path.extend(make_arc((75.0 - r, 25.0 + r), r, 3.0*PI/2.0, 2.0*PI, 16));
    bot_path.push((75.0, 50.0));
    add_ch(1, 2, ChannelType::Serpentine { path: bot_path });
    
    // Merge to Outlet
    add_ch(2, 3, ChannelType::Straight);
    let box_outline = vec![
        ((0.0, 0.0), (100.0, 0.0)),
        ((100.0, 0.0), (100.0, 100.0)),
        ((100.0, 100.0), (0.0, 100.0)),
        ((0.0, 100.0), (0.0, 0.0)),
    ];

    let system = ChannelSystem {
        box_dims,
        nodes,
        channels,
        box_outline,
    };

    // Export to JSON using the interchange format (which cfd-mesh expects)
    let json_str = system.to_interchange_json()?;
    let json_path = out_dir.join("schematic.json");
    fs::write(&json_path, json_str)?;

    // Export to PNG & SVG
    let png_path = out_dir.join("schematic.png");
    let svg_path = out_dir.join("schematic.svg");
    
    // We need to pass &str to plot_geometry
    plot_geometry(&system, png_path.to_str().unwrap())?;
    plot_geometry(&system, svg_path.to_str().unwrap())?;

    println!("\n✅ Mirrored Bifurcation successfully exported to:");
    println!("  - {}", json_path.display());
    println!("  - {}", png_path.display());
    println!("  - {}", svg_path.display());

    Ok(())
}
