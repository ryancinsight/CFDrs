use cfd_schematics::domain::model::{ChannelSpec, NetworkBlueprint, NodeId, NodeKind, NodeSpec};
use cfd_schematics::visualizations::{
    AnnotationMarker, MarkerRole, RenderConfig, SchematicAnnotations,
};
use cfd_schematics::{plot_geometry, plot_geometry_with_annotations, plot_geometry_with_config};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn synthetic_system() -> NetworkBlueprint {
    let nodes = vec![
        NodeSpec {
            id: NodeId::new("0"),
            kind: NodeKind::Inlet,
            point: (0.0, 42.735),
            layout: None,
            junction_geometry: None,
            metadata: None,
        },
        NodeSpec {
            id: NodeId::new("1"),
            kind: NodeKind::Junction,
            point: (30.0, 42.735),
            layout: None,
            junction_geometry: None,
            metadata: None,
        },
        NodeSpec {
            id: NodeId::new("2"),
            kind: NodeKind::Junction,
            point: (60.0, 52.0),
            layout: None,
            junction_geometry: None,
            metadata: None,
        },
        NodeSpec {
            id: NodeId::new("3"),
            kind: NodeKind::Junction,
            point: (60.0, 33.0),
            layout: None,
            junction_geometry: None,
            metadata: None,
        },
        NodeSpec {
            id: NodeId::new("4"),
            kind: NodeKind::Junction,
            point: (92.0, 42.735),
            layout: None,
            junction_geometry: None,
            metadata: None,
        },
        NodeSpec {
            id: NodeId::new("5"),
            kind: NodeKind::Outlet,
            point: (127.76, 42.735),
            layout: None,
            junction_geometry: None,
            metadata: None,
        },
    ];

    let channels = vec![
        ChannelSpec::new_pipe_rect("c0", "0", "1", 30.0, 1.0, 0.8, 10.0, 0.0),
        ChannelSpec::new_pipe_rect("c1", "1", "2", 30.0, 1.0, 0.8, 10.0, 0.0),
        ChannelSpec::new_pipe_rect("c2", "1", "3", 30.0, 1.0, 0.8, 10.0, 0.0),
        ChannelSpec::new_pipe_rect("c3", "2", "4", 30.0, 1.0, 0.8, 10.0, 0.0),
        ChannelSpec::new_pipe_rect("c4", "3", "4", 30.0, 1.0, 0.8, 10.0, 0.0),
        ChannelSpec::new_pipe_rect("c5", "4", "5", 30.0, 1.0, 0.8, 10.0, 0.0),
    ];

    NetworkBlueprint {
        name: "test".to_string(),
        box_dims: (127.76, 85.47),
        nodes,
        channels,
        box_outline: vec![
            ((0.0, 0.0), (127.76, 0.0)),
            ((127.76, 0.0), (127.76, 85.47)),
            ((127.76, 85.47), (0.0, 85.47)),
            ((0.0, 85.47), (0.0, 0.0)),
        ],
        render_hints: None,
        topology: None,
        lineage: None,
        metadata: None,
        geometry_authored: false,
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

#[test]
fn auto_annotations_include_computed_volume_label() {
    let system = synthetic_system();
    let path = unique_svg_path("cfd_schematic_auto_volume");
    let path_str = path.to_string_lossy();

    plot_geometry(&system, path_str.as_ref()).expect("auto-annotated render must succeed");

    let svg = std::fs::read_to_string(&path).expect("must read rendered svg");
    assert!(svg.contains("Volume:"));
    assert!(!svg.contains("Volume: --"));
}
