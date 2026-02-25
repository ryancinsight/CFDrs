//! Integration tests for the NetworkBlueprint → IndexedMesh pipeline.
#![cfg(feature = "scheme-io")]

use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig, TopologyClass};
use cfd_schematics::interface::presets::{
    bifurcation_rect, serpentine_chain, serpentine_rect, symmetric_bifurcation,
    symmetric_trifurcation, trifurcation_rect, venturi_chain, venturi_rect,
};

fn default_cfg() -> PipelineConfig {
    PipelineConfig::default()
}

fn no_chip_cfg() -> PipelineConfig {
    PipelineConfig {
        include_chip_body: false,
        ..Default::default()
    }
}

// ── Watertight fluid mesh tests ───────────────────────────────────────────────

#[test]
fn venturi_chain_fluid_mesh_watertight() {
    let bp = venturi_chain("v1", 0.030, 0.004, 0.002);
    let mut out = BlueprintMeshPipeline::run(&bp, &no_chip_cfg())
        .expect("venturi_chain pipeline failed");
    assert!(
        out.fluid_mesh.is_watertight(),
        "venturi fluid mesh must be watertight"
    );
    assert!(
        out.fluid_mesh.signed_volume() > 0.0,
        "venturi fluid mesh must have positive volume"
    );
    assert_eq!(out.topology_class, TopologyClass::VenturiChain);
    assert_eq!(out.segment_count, 3);
    // At least one face labelled "inlet" and one "outlet"
    let has_inlet = out.fluid_mesh.boundary_labels.values().any(|l| l == "inlet");
    assert!(has_inlet, "venturi fluid mesh must have inlet label");
}

#[test]
fn bifurcation_fluid_mesh_watertight() {
    let bp = symmetric_bifurcation("b1", 0.010, 0.010, 0.004, 0.003);
    let mut out = BlueprintMeshPipeline::run(&bp, &no_chip_cfg())
        .expect("symmetric_bifurcation pipeline failed");
    assert!(
        out.fluid_mesh.is_watertight(),
        "bifurcation fluid mesh must be watertight"
    );
    assert_eq!(out.topology_class, TopologyClass::Bifurcation);
    let has_inlet = out.fluid_mesh.boundary_labels.values().any(|l| l == "inlet");
    assert!(has_inlet, "bifurcation fluid mesh must have inlet label");
}

#[test]
fn trifurcation_fluid_mesh_watertight() {
    let bp = symmetric_trifurcation("t1", 0.010, 0.008, 0.004, 0.004);
    let mut out = BlueprintMeshPipeline::run(&bp, &no_chip_cfg())
        .expect("symmetric_trifurcation pipeline failed");
    assert!(
        out.fluid_mesh.is_watertight(),
        "trifurcation fluid mesh must be watertight"
    );
    assert_eq!(out.topology_class, TopologyClass::Trifurcation);
}

#[test]
fn serpentine_chain_fluid_mesh_watertight() {
    let bp = serpentine_chain("s1", 3, 0.010, 0.004);
    let mut out = BlueprintMeshPipeline::run(&bp, &no_chip_cfg())
        .expect("serpentine_chain pipeline failed");
    assert!(
        out.fluid_mesh.is_watertight(),
        "serpentine fluid mesh must be watertight"
    );
    assert_eq!(out.segment_count, 3);
}

// ── Chip body tests ───────────────────────────────────────────────────────────

#[test]
fn venturi_chip_body_produced() {
    let bp = venturi_chain("v1", 0.030, 0.004, 0.002);
    let mut out = BlueprintMeshPipeline::run(&bp, &default_cfg())
        .expect("venturi_chain pipeline (with chip body) failed");
    assert!(out.chip_mesh.is_some(), "chip_mesh should be Some");
    let chip = out.chip_mesh.as_mut().unwrap();
    assert!(
        chip.is_watertight(),
        "chip body mesh must be watertight"
    );
}

#[test]
fn chip_body_volume_less_than_substrate() {
    let bp = venturi_chain("v1", 0.030, 0.004, 0.002);
    let mut out = BlueprintMeshPipeline::run(&bp, &default_cfg())
        .expect("venturi pipeline with chip body failed");

    let chip_vol = out.chip_mesh.as_mut().unwrap().signed_volume();
    // Substrate volume using the default chip height from PipelineConfig
    let substrate_vol = 127.76_f64 * 85.47 * default_cfg().chip_height_mm;
    assert!(
        chip_vol < substrate_vol,
        "chip body volume ({chip_vol:.2}) must be less than substrate ({substrate_vol:.2})"
    );
    assert!(chip_vol > 0.0, "chip body volume must be positive");
}

// ── Constraint rejection tests ────────────────────────────────────────────────

#[test]
fn pipeline_rejects_wrong_diameter() {
    // 2 mm < 4 mm ± 0.1 mm → should fail
    let bp = serpentine_chain("x", 3, 0.010, 0.002);
    let result = BlueprintMeshPipeline::run(&bp, &default_cfg());
    assert!(result.is_err(), "2mm diameter should be rejected");
    let msg = result.err().expect("checked above").to_string();
    assert!(
        msg.to_lowercase().contains("hydraulic diameter")
            || msg.to_lowercase().contains("channel error"),
        "error should mention diameter: {msg}"
    );
}

#[test]
fn pipeline_rejects_channels_outside_plate() {
    // Channel longer than plate (127.76 mm) will violate wall clearance
    // total_length = 200 mm + 200 mm + 200 mm > 127.76 mm plate
    let bp = venturi_chain("long", 0.600, 0.004, 0.002); // 600 mm total
    let result = BlueprintMeshPipeline::run(&bp, &default_cfg());
    assert!(
        result.is_err(),
        "channels extending past plate should be rejected"
    );
}

// ── Classification tests ──────────────────────────────────────────────────────

#[test]
fn venturi_rect_classifies_as_venturi() {
    use cfd_mesh::application::pipeline::NetworkTopology;
    let bp = venturi_rect("vr1", 0.004, 0.002, 0.004, 0.005);
    let topo = NetworkTopology::new(&bp);
    assert_eq!(topo.classify(), TopologyClass::VenturiChain);
}

#[test]
fn bifurcation_rect_classifies_correctly() {
    use cfd_mesh::application::pipeline::NetworkTopology;
    let bp = bifurcation_rect("br1", 0.010, 0.010, 0.004, 0.003, 0.004);
    let topo = NetworkTopology::new(&bp);
    assert_eq!(topo.classify(), TopologyClass::Bifurcation);
}

// ── All six therapy designs ───────────────────────────────────────────────────

#[test]
fn all_six_designs_watertight_and_positive_volume() {
    let designs: Vec<(&str, cfd_schematics::NetworkBlueprint)> = vec![
        ("venturi_chain", venturi_chain("d1", 0.030, 0.004, 0.002)),
        (
            "symmetric_bifurcation",
            symmetric_bifurcation("d2", 0.010, 0.010, 0.004, 0.003),
        ),
        (
            "symmetric_trifurcation",
            symmetric_trifurcation("d3", 0.010, 0.008, 0.004, 0.004),
        ),
        ("serpentine_chain", serpentine_chain("d4", 3, 0.010, 0.004)),
        ("venturi_rect", venturi_rect("d5", 0.004, 0.002, 0.004, 0.005)),
        (
            "serpentine_rect",
            serpentine_rect("d6", 3, 0.010, 0.004, 0.004),
        ),
    ];

    for (name, bp) in designs {
        let result = BlueprintMeshPipeline::run(&bp, &no_chip_cfg());
        match result {
            Ok(mut out) => {
                assert!(
                    out.fluid_mesh.is_watertight(),
                    "{name}: fluid_mesh must be watertight"
                );
                assert!(
                    out.fluid_mesh.signed_volume() > 0.0,
                    "{name}: fluid_mesh must have positive volume"
                );
            }
            Err(e) => {
                panic!("{name}: pipeline failed: {e}");
            }
        }
    }
}

// ── trifurcation_rect test ────────────────────────────────────────────────────

#[test]
fn trifurcation_rect_produces_watertight_mesh() {
    let bp = trifurcation_rect("tr1", 0.010, 0.008, 0.004, 0.004, 0.004);
    let mut out = BlueprintMeshPipeline::run(&bp, &no_chip_cfg())
        .expect("trifurcation_rect pipeline failed");
    assert!(
        out.fluid_mesh.is_watertight(),
        "trifurcation_rect fluid mesh must be watertight"
    );
    assert_eq!(out.topology_class, TopologyClass::Trifurcation);
}
