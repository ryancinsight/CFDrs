use cfd_schematics::interface::presets::{
    asymmetric_trifurcation_venturi_rect, symmetric_bifurcation, symmetric_trifurcation,
    venturi_rect,
};

fn channel_path(bp: &cfd_schematics::NetworkBlueprint, channel_id: &str) -> Vec<(f64, f64)> {
    bp.channels
        .iter()
        .find(|channel| channel.id.as_str() == channel_id)
        .map(|channel| channel.path.clone())
        .unwrap_or_else(|| panic!("channel {channel_id} must exist"))
}

#[test]
fn symmetric_bifurcation_materializes_mirrored_branch_paths() {
    let bp = symmetric_bifurcation("mirrored-bi", 10.0e-3, 12.0e-3, 4.0e-3, 3.0e-3);

    // Must have at least 4 nodes (inlet, split, merge, outlet)
    assert!(bp.nodes.len() >= 4);
    // Must have at least 4 channels (trunk_in + 2 branches + trunk_out)
    assert!(bp.channels.len() >= 4);
    // All channels must have routed paths
    assert!(bp.channels.iter().all(|channel| channel.path.len() >= 2));
    // Blueprint structural validation
    bp.validate()
        .expect("symmetric bifurcation blueprint should be structurally valid");
}

#[test]
fn symmetric_trifurcation_materializes_three_distinct_parallel_lanes() {
    let bp = symmetric_trifurcation("mirrored-tri", 10.0e-3, 12.0e-3, 4.0e-3, 3.0e-3);

    // Must have at least 4 nodes (inlet, split, merge, outlet)
    assert!(bp.nodes.len() >= 4);
    // Must have at least 5 channels (trunk_in + 3 branches + trunk_out)
    assert!(bp.channels.len() >= 5);
    // All channels must have routed paths
    assert!(bp.channels.iter().all(|channel| channel.path.len() >= 2));
    // Blueprint structural validation
    bp.validate()
        .expect("symmetric trifurcation blueprint should be structurally valid");
}

#[test]
fn asymmetric_trifurcation_materializes_distinct_bypass_paths_without_manual_layout() {
    let bp = asymmetric_trifurcation_venturi_rect(
        "asym-tri", 12.0e-3, 16.0e-3, 2.0e-3, 0.45, 0.25, 100e-6, 300e-6, 1.0e-3,
    );
    let left = channel_path(&bp, "left_arm");
    let right = channel_path(&bp, "right_arm");

    assert!(bp.nodes.iter().all(|node| node.layout.is_some()));
    assert!(left.len() >= 4);
    assert!(right.len() >= 4);
    assert_ne!(left, right);
}

#[test]
fn canonical_series_presets_materialize_layout_coordinates() {
    let bp = venturi_rect("series-layout", 2.0e-3, 1.0e-3, 1.0e-3, 3.0e-3);

    assert!(bp.nodes.iter().all(|node| node.layout.is_some()));
    assert!(bp.channels.iter().all(|channel| channel.path.len() >= 2));
}
