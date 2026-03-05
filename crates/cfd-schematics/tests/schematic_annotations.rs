use cfd_schematics::geometry::{Channel, ChannelSystem, ChannelType, Node};
use cfd_schematics::visualizations::{
    AnnotationMarker, MarkerRole, RenderConfig, SchematicAnnotations,
};
use cfd_schematics::{plot_geometry_with_annotations, plot_geometry_with_config};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn synthetic_system() -> ChannelSystem {
    let nodes = vec![
        Node {
            id: 0,
            point: (0.0, 42.735),
            metadata: None,
        },
        Node {
            id: 1,
            point: (30.0, 42.735),
            metadata: None,
        },
        Node {
            id: 2,
            point: (60.0, 52.0),
            metadata: None,
        },
        Node {
            id: 3,
            point: (60.0, 33.0),
            metadata: None,
        },
        Node {
            id: 4,
            point: (92.0, 42.735),
            metadata: None,
        },
        Node {
            id: 5,
            point: (127.76, 42.735),
            metadata: None,
        },
    ];

    let channels = vec![
        Channel {
            id: 0,
            from_node: 0,
            to_node: 1,
            width: 1.0,
            height: 0.8,
            channel_type: ChannelType::Straight,
            metadata: None,
        },
        Channel {
            id: 1,
            from_node: 1,
            to_node: 2,
            width: 1.0,
            height: 0.8,
            channel_type: ChannelType::Straight,
            metadata: None,
        },
        Channel {
            id: 2,
            from_node: 1,
            to_node: 3,
            width: 1.0,
            height: 0.8,
            channel_type: ChannelType::Straight,
            metadata: None,
        },
        Channel {
            id: 3,
            from_node: 2,
            to_node: 4,
            width: 1.0,
            height: 0.8,
            channel_type: ChannelType::Straight,
            metadata: None,
        },
        Channel {
            id: 4,
            from_node: 3,
            to_node: 4,
            width: 1.0,
            height: 0.8,
            channel_type: ChannelType::Straight,
            metadata: None,
        },
        Channel {
            id: 5,
            from_node: 4,
            to_node: 5,
            width: 1.0,
            height: 0.8,
            channel_type: ChannelType::Straight,
            metadata: None,
        },
    ];

    ChannelSystem {
        box_dims: (127.76, 85.47),
        nodes,
        channels,
        box_outline: vec![
            ((0.0, 0.0), (127.76, 0.0)),
            ((127.76, 0.0), (127.76, 85.47)),
            ((127.76, 85.47), (0.0, 85.47)),
            ((0.0, 85.47), (0.0, 0.0)),
        ],
    }
}

fn unique_svg_path(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}.svg"))
}

#[test]
fn annotated_svg_contains_markers_labels_and_role_colors() {
    let system = synthetic_system();
    let mut annotations = SchematicAnnotations::report_default();
    annotations.markers = vec![
        AnnotationMarker::new((0.0, 42.735), MarkerRole::Inlet).with_label("IN", true),
        AnnotationMarker::new((127.76, 42.735), MarkerRole::Outlet).with_label("OUT", true),
        AnnotationMarker::new((60.0, 42.735), MarkerRole::VenturiThroat).with_label("TH1", true),
    ];

    let config = RenderConfig::well_plate_96_report_annotated();
    let path = unique_svg_path("cfd_schematic_annotations");
    let path_str = path.to_string_lossy();

    plot_geometry_with_annotations(&system, path_str.as_ref(), &config, &annotations)
        .expect("annotated schematic render must succeed");

    let svg = std::fs::read_to_string(&path).expect("must read rendered svg");
    let svg_lower = svg.to_ascii_lowercase();

    assert!(svg.contains("IN"));
    assert!(svg.contains("OUT"));
    assert!(svg.contains("TH1"));
    assert!(svg_lower.contains("circle") || svg_lower.contains("path"));
    assert!(svg_lower.contains("#0b7285"));
    assert!(svg_lower.contains("#c2410c"));
}

#[test]
fn plot_geometry_with_config_still_renders_without_annotations() {
    let system = synthetic_system();
    let config = RenderConfig::well_plate_96();
    let path = unique_svg_path("cfd_schematic_plain");
    let path_str = path.to_string_lossy();

    plot_geometry_with_config(&system, path_str.as_ref(), &config)
        .expect("non-annotated render must still succeed");

    let svg = std::fs::read_to_string(&path).expect("must read rendered svg");
    assert!(svg.contains("<svg"));
}
