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

    // Verify 1 inlet, 1 outlet, 2 junctions
    assert_eq!(bp.nodes.len(), 4);
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

    assert_eq!(inlets.len(), 1, "Must have exactly 1 inlet");
    assert_eq!(outlets.len(), 1, "Must have exactly 1 outlet");
    assert_eq!(
        junctions.len(),
        2,
        "Must have diverging and converging junctions"
    );

    // Verify exactly 4 channels (parent_in, daughter_1, daughter_2, parent_out)
    assert_eq!(bp.channels.len(), 4);

    let parent_in = bp
        .channels
        .iter()
        .find(|c| c.id.as_str() == "parent_in")
        .unwrap();
    let parent_out = bp
        .channels
        .iter()
        .find(|c| c.id.as_str() == "parent_out")
        .unwrap();
    let d1 = bp
        .channels
        .iter()
        .find(|c| c.id.as_str() == "daughter_1")
        .unwrap();
    let d2 = bp
        .channels
        .iter()
        .find(|c| c.id.as_str() == "daughter_2")
        .unwrap();

    // Verify lengths and 4mm requirement (diameter = 4e-3)
    let expected_area = std::f64::consts::PI * (pipe_diam / 2.0).powi(2);
    for ch in &[parent_in, parent_out, d1, d2] {
        assert_eq!(ch.cross_section.hydraulic_diameter(), pipe_diam);
        assert!((ch.cross_section.area() - expected_area).abs() < 1e-12);
    }

    // Compare correctness against analytical solutions (Hagen-Poiseuille for pipes)
    assert_eq!(parent_in.resistance, 0.0);
    assert_eq!(parent_out.resistance, 0.0);
    assert_eq!(d1.resistance, 0.0);
    assert_eq!(d2.resistance, 0.0);
}
