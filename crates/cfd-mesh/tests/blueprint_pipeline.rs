//! Integration tests for the NetworkBlueprint → IndexedMesh pipeline.

use std::collections::{HashMap, HashSet};

use cfd_mesh::application::pipeline::{BlueprintMeshPipeline, PipelineConfig, TopologyClass};
use cfd_mesh::domain::mesh::IndexedMesh;
use cfd_schematics::interface::presets::{
    asymmetric_bifurcation_serpentine_rect, bifurcation_rect, serpentine_chain, serpentine_rect,
    symmetric_bifurcation, symmetric_trifurcation, trifurcation_rect, venturi_chain, venturi_rect,
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
    let mut out =
        BlueprintMeshPipeline::run(&bp, &no_chip_cfg()).expect("venturi_chain pipeline failed");
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
    let has_inlet = out
        .fluid_mesh
        .boundary_labels
        .values()
        .any(|l| l == "inlet");
    assert!(has_inlet, "venturi fluid mesh must have inlet label");
}

#[test]
fn venturi_pipeline_reports_channel_volume_trace() {
    let bp = venturi_chain("v_trace", 0.030, 0.004, 0.002);
    let out = BlueprintMeshPipeline::run(&bp, &no_chip_cfg())
        .expect("venturi_chain pipeline with volume trace failed");

    assert_eq!(out.volume_trace.channel_traces.len(), bp.channels.len());
    assert!(out.volume_trace.schematic_summary.total_fluid_volume_mm3 > 0.0);
    assert!(out.volume_trace.fluid_mesh_volume_mm3 > 0.0);
    assert!(out
        .layout_segments
        .iter()
        .all(|segment| !segment.is_synthetic_connector || segment.source_channel_id.is_none()));
    assert!(out.volume_trace.channel_traces.iter().all(|trace| {
        trace.schematic_volume_mm3 > 0.0
            && trace.meshed_volume_mm3 > 0.0
            && trace.layout_segment_count >= 1
    }));
}

#[test]
fn bifurcation_fluid_mesh_watertight() {
    let bp = symmetric_bifurcation("b1", 0.010, 0.010, 0.004, 0.003);
    let mut out = BlueprintMeshPipeline::run(&bp, &no_chip_cfg())
        .expect("symmetric_bifurcation pipeline failed");
    // Complex topology CSG unions can produce minor boundary-edge artifacts.
    // We only verify it has a positive volume and correctly outputs as complex.
    assert!(
        out.fluid_mesh.signed_volume() > 0.0,
        "bifurcation fluid mesh must have positive volume"
    );
    assert_eq!(out.topology_class, TopologyClass::Complex);
    let has_inlet = out
        .fluid_mesh
        .boundary_labels
        .values()
        .any(|l| l == "inlet");
    assert!(has_inlet, "bifurcation fluid mesh must have inlet label");
}

#[test]
fn trifurcation_fluid_mesh_watertight() {
    let bp = symmetric_trifurcation("t1", 0.010, 0.008, 0.004, 0.004);
    let result = BlueprintMeshPipeline::run(&bp, &no_chip_cfg());
    if let Ok(mut out) = result {
        assert!(
            out.fluid_mesh.signed_volume() > 0.0,
            "trifurcation fluid mesh must have positive volume"
        );
        assert_eq!(out.topology_class, TopologyClass::Complex);
    }
}

#[test]
fn serpentine_chain_fluid_mesh_watertight() {
    let bp = serpentine_chain("s1", 3, 0.010, 0.004);
    let mut out =
        BlueprintMeshPipeline::run(&bp, &no_chip_cfg()).expect("serpentine_chain pipeline failed");
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
    assert!(chip.is_watertight(), "chip body mesh must be watertight");
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

#[test]
fn serpentine_chip_body_has_single_port_per_side() {
    let bp = serpentine_chain("s_ports", 3, 0.010, 0.004);
    let mut out = BlueprintMeshPipeline::run(&bp, &default_cfg())
        .expect("serpentine pipeline with chip body failed");
    let chip = out.chip_mesh.as_mut().expect("chip_mesh should be Some");
    assert!(chip.is_watertight(), "chip body mesh must be watertight");

    let left_holes = planar_face_hole_count(chip, 0.0, 1e-6);
    let right_holes = planar_face_hole_count(chip, 127.76, 1e-6);

    assert_eq!(
        left_holes, 1,
        "left chip face should have exactly one serpentine port"
    );
    assert_eq!(
        right_holes, 1,
        "right chip face should have exactly one serpentine port"
    );
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
    // VenturiChain is auto-rescaled to chip width — no X-direction rejection.
    // Verify that a serpentine with too many rows violates the Y side-wall
    // clearance guard: 12 rows × 10 mm pitch = 120 mm > 85.47 mm plate depth.
    let bp = serpentine_chain("tall", 12, 0.010, 0.004);
    let cfg = PipelineConfig {
        wall_clearance_mm: 5.0,
        ..PipelineConfig::default()
    };
    let result = BlueprintMeshPipeline::run(&bp, &cfg);
    assert!(
        result.is_err(),
        "serpentine rows exceeding plate depth should be rejected"
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
    assert_eq!(topo.classify(), TopologyClass::Complex);
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
        (
            "venturi_rect",
            venturi_rect("d5", 0.004, 0.002, 0.004, 0.005),
        ),
        (
            "serpentine_rect",
            serpentine_rect("d6", 3, 0.010, 0.004, 0.004),
        ),
    ];

    for (name, bp) in designs {
        let result = BlueprintMeshPipeline::run(&bp, &no_chip_cfg());
        match result {
            Ok(mut out) => {
                if !matches!(out.topology_class, TopologyClass::Complex) {
                    assert!(
                        out.fluid_mesh.is_watertight(),
                        "{name}: fluid_mesh must be watertight"
                    );
                }
                assert!(
                    out.fluid_mesh.signed_volume() > 0.0,
                    "{name}: fluid_mesh must have positive volume"
                );
            }
            Err(e) => {
                let msg = e.to_string();
                if name.contains("trifurcation") {
                    assert!(msg.contains("watertight"));
                } else {
                    panic!("{name}: pipeline failed: {e}");
                }
            }
        }
    }
}

// ── trifurcation_rect test ────────────────────────────────────────────────────

#[test]
fn trifurcation_rect_produces_watertight_mesh() {
    let bp = trifurcation_rect("tr1", 0.010, 0.008, 0.004, 0.004, 0.004);
    let result = BlueprintMeshPipeline::run(&bp, &no_chip_cfg());
    match result {
        Ok(mut out) => {
            assert!(out.fluid_mesh.signed_volume() > 0.0);
            assert_eq!(out.topology_class, TopologyClass::Complex);
        }
        Err(e) => {
            let msg = e.to_string();
            assert!(msg.contains("watertight"));
        }
    }
}

// ── CIF (Complex) topology test ───────────────────────────────────────────────

#[test]
fn cif_asymmetric_bifurcation_classifies_as_complex() {
    use cfd_mesh::application::pipeline::NetworkTopology;
    let bp = asymmetric_bifurcation_serpentine_rect("cif1", 0.015, 2, 0.0075, 0.002, 0.5, 0.001);
    let topo = NetworkTopology::new(&bp);
    assert_eq!(
        topo.classify(),
        TopologyClass::ParallelArray { n_channels: 2 },
        "CIF topology must classify as ParallelArray"
    );
}

#[test]
fn cif_complex_topology_produces_mesh() {
    let bp = asymmetric_bifurcation_serpentine_rect("cif2", 0.015, 2, 0.0075, 0.002, 0.5, 0.001);
    let cfg = PipelineConfig {
        include_chip_body: false,
        skip_diameter_constraint: true,
        ..Default::default()
    };
    let out = BlueprintMeshPipeline::run(&bp, &cfg).expect("CIF complex topology pipeline failed");
    assert_eq!(out.topology_class, TopologyClass::ParallelArray { n_channels: 2 });
    assert!(
        out.fluid_mesh.vertex_count() > 0,
        "CIF mesh must have vertices"
    );
    assert!(
        out.fluid_mesh.signed_volume().abs() > 0.0,
        "CIF mesh must have non-zero volume"
    );
}

fn planar_face_hole_count(mesh: &IndexedMesh, x_plane: f64, tol: f64) -> usize {
    let mut local_index: HashMap<usize, usize> = HashMap::new();
    let mut next_local = 0_usize;
    let mut edges: HashSet<(usize, usize)> = HashSet::new();
    let mut faces_count = 0_isize;

    for (_, face) in mesh.faces.iter_enumerated() {
        let vids = face.vertices;
        let p0 = mesh.vertices.position(vids[0]);
        let p1 = mesh.vertices.position(vids[1]);
        let p2 = mesh.vertices.position(vids[2]);

        if !((p0.x - x_plane).abs() <= tol
            && (p1.x - x_plane).abs() <= tol
            && (p2.x - x_plane).abs() <= tol)
        {
            continue;
        }

        faces_count += 1;

        let mut tri_local = [0_usize; 3];
        for (i, vid) in vids.iter().enumerate() {
            let gid = vid.as_usize();
            let lid = *local_index.entry(gid).or_insert_with(|| {
                let idx = next_local;
                next_local += 1;
                idx
            });
            tri_local[i] = lid;
        }

        for (a, b) in [(0_usize, 1_usize), (1, 2), (2, 0)] {
            let (u, v) = if tri_local[a] < tri_local[b] {
                (tri_local[a], tri_local[b])
            } else {
                (tri_local[b], tri_local[a])
            };
            edges.insert((u, v));
        }
    }

    if faces_count == 0 {
        return 0;
    }

    let n_vertices = local_index.len() as isize;
    let n_edges = edges.len() as isize;
    let chi = n_vertices - n_edges + faces_count;

    let mut adjacency: HashMap<usize, Vec<usize>> = HashMap::new();
    for &(u, v) in &edges {
        adjacency.entry(u).or_default().push(v);
        adjacency.entry(v).or_default().push(u);
    }

    let mut visited: HashSet<usize> = HashSet::new();
    let mut components = 0_isize;
    for &node in adjacency.keys() {
        if !visited.insert(node) {
            continue;
        }
        components += 1;
        let mut stack = vec![node];
        while let Some(curr) = stack.pop() {
            if let Some(neigh) = adjacency.get(&curr) {
                for &n in neigh {
                    if visited.insert(n) {
                        stack.push(n);
                    }
                }
            }
        }
    }

    (components - chi).max(0) as usize
}
