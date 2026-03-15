use cfd_schematics::domain::model::NodeKind;
use cfd_schematics::interface::presets::symmetric_bifurcation;

#[test]
fn test_symmetric_bifurcation_constraints() {
    let parent_len = 10.0e-3; // 10 mm
    let daughter_len = 15.0e-3; // 15 mm
    let pipe_diam = 4.0e-3; // 4 mm

    let bp = symmetric_bifurcation(
        "therapy_bifurcation",
        parent_len,
        daughter_len,
        pipe_diam,
        pipe_diam,
    );

    // Verify topology: at least 1 inlet, 1 outlet, and junctions.
    // Geometry-authored blueprints include waypoints, so node count > 4.
    let inlets: Vec<_> = bp
        .nodes
        .iter()
        .filter(|n| n.kind == NodeKind::Inlet)
        .collect();
    let outlets: Vec<_> = bp
        .nodes
        .iter()
        .filter(|n| n.kind == NodeKind::Outlet)
        .collect();
    let junctions: Vec<_> = bp
        .nodes
        .iter()
        .filter(|n| n.kind == NodeKind::Junction)
        .collect();

    assert!(
        !inlets.is_empty(),
        "Must have at least 1 inlet, got {}",
        inlets.len()
    );
    assert!(
        !outlets.is_empty(),
        "Must have at least 1 outlet, got {}",
        outlets.len()
    );
    assert!(
        junctions.len() >= 2,
        "Must have at least 2 junctions (diverging + converging), got {}",
        junctions.len()
    );

    // Verify channels exist and have valid geometry.
    assert!(
        bp.channels.len() >= 4,
        "Must have at least 4 channels (parent_in, daughters, parent_out), got {}",
        bp.channels.len()
    );

    // All channels must have positive length and area.
    for ch in &bp.channels {
        assert!(
            ch.length_m > 0.0,
            "Channel '{}' must have positive length",
            ch.id.as_str()
        );
        assert!(
            ch.cross_section.area() > 0.0,
            "Channel '{}' must have positive area",
            ch.id.as_str()
        );
    }
}
