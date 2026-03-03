//! Phase 7 — End-to-end pipeline integration tests.
//!
//! These tests exercise the full path from `DesignCandidate` to physics
//! metrics, and from `NetworkBlueprint` to watertight mesh output.

use cfd_optim::{compute_metrics, CrossSectionShape, DesignCandidate, DesignTopology};

// ---------------------------------------------------------------------------
// Helper: build a minimal SingleVenturi candidate with reasonable defaults
// ---------------------------------------------------------------------------

fn single_venturi_candidate() -> DesignCandidate {
    DesignCandidate {
        id: "test-sv".into(),
        topology: DesignTopology::SingleVenturi,
        flow_rate_m3_s: 1e-7,             // 6 mL/min
        inlet_gauge_pa: 5000.0,           // 5 kPa gauge
        throat_diameter_m: 50e-6,         // 50 µm throat
        inlet_diameter_m: 4e-3,           // 4 mm inlet
        throat_length_m: 100e-6,          // 100 µm throat length
        channel_width_m: 2e-3,            // 2 mm channel width
        channel_height_m: 0.5e-3,         // 0.5 mm channel height
        serpentine_segments: 4,
        segment_length_m: 10e-3,          // 10 mm segments
        bend_radius_m: 1e-3,             // 1 mm bend radius
        feed_hematocrit: 0.45,            // whole blood
        trifurcation_center_frac: 1.0 / 3.0,
        cif_pretri_center_frac: 1.0 / 3.0,
        cif_terminal_tri_center_frac: 1.0 / 3.0,
        cif_terminal_bi_treat_frac: 0.68,
        asymmetric_narrow_frac: 0.5,
        trifurcation_left_frac: 1.0 / 3.0,
        cross_section_shape: CrossSectionShape::Rectangular,
    }
}

// ---------------------------------------------------------------------------
// Test 1: SingleVenturi candidate computes metrics successfully
// ---------------------------------------------------------------------------

#[test]
fn single_venturi_candidate_computes_metrics() {
    let candidate = single_venturi_candidate();
    let result = compute_metrics(&candidate);
    assert!(
        result.is_ok(),
        "compute_metrics must succeed: {:?}",
        result.err()
    );

    let metrics = result.unwrap();
    // Cavitation number should be finite for a venturi topology
    assert!(
        metrics.cavitation_number.is_finite(),
        "cavitation_number must be finite, got {}",
        metrics.cavitation_number
    );
    // Total pressure drop should be positive
    assert!(
        metrics.total_pressure_drop_pa > 0.0,
        "total_pressure_drop_pa must be positive, got {}",
        metrics.total_pressure_drop_pa
    );
}

// ---------------------------------------------------------------------------
// Test 2: candidate -> blueprint round-trip
// ---------------------------------------------------------------------------

#[test]
fn candidate_to_blueprint_round_trip() {
    let candidate = single_venturi_candidate();
    let bp = candidate.to_blueprint();

    assert!(!bp.nodes.is_empty(), "Blueprint must have nodes");
    assert!(!bp.channels.is_empty(), "Blueprint must have channels");

    // SingleVenturi should have at least inlet, contraction, throat, outlet
    assert!(
        bp.nodes.len() >= 4,
        "SingleVenturi blueprint must have >= 4 nodes, got {}",
        bp.nodes.len()
    );
    // Three channels: inlet_section, throat_section, diffuser_section
    assert!(
        bp.channels.len() >= 3,
        "SingleVenturi blueprint must have >= 3 channels, got {}",
        bp.channels.len()
    );
}

// ---------------------------------------------------------------------------
// Test 3: BlueprintMeshPipeline produces watertight fluid mesh.
// ---------------------------------------------------------------------------

#[test]
fn blueprint_mesh_pipeline_produces_watertight() {
    use cfd_mesh::application::pipeline::blueprint_mesh::{BlueprintMeshPipeline, PipelineConfig};
    use cfd_schematics::interface::presets::venturi_chain;

    let bp = venturi_chain("test", 0.030, 0.004, 0.002);
    let config = PipelineConfig::default();
    let mut output = BlueprintMeshPipeline::run(&bp, &config)
        .expect("BlueprintMeshPipeline::run must succeed");

    assert!(
        output.fluid_mesh.is_watertight(),
        "Fluid mesh must be watertight"
    );
    assert!(
        output.fluid_mesh.face_count() > 0,
        "Fluid mesh must have faces, got {}",
        output.fluid_mesh.face_count()
    );
}
