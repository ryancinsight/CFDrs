use cfd_optim::{
    save_schematic_svg, CrossSectionShape, DesignCandidate, DesignTopology, PrimitiveSplitSequence,
    TreatmentZoneMode,
};
use std::path::PathBuf;
use std::time::{SystemTime, UNIX_EPOCH};

fn unique_svg_path(prefix: &str) -> PathBuf {
    let nanos = SystemTime::now()
        .duration_since(UNIX_EPOCH)
        .expect("system clock before unix epoch")
        .as_nanos();
    std::env::temp_dir().join(format!("{prefix}_{nanos}.svg"))
}

fn candidate_with_topology(id: &str, topology: DesignTopology) -> DesignCandidate {
    DesignCandidate {
        id: id.to_string(),
        topology,
        flow_rate_m3_s: 150.0 / 6.0e7,
        inlet_gauge_pa: 300_000.0,
        throat_diameter_m: 45e-6,
        inlet_diameter_m: 540e-6,
        throat_length_m: 250e-6,
        channel_width_m: 5e-3,
        channel_height_m: 1.5e-3,
        serpentine_segments: 4,
        segment_length_m: 8e-3,
        bend_radius_m: 3e-3,
        feed_hematocrit: 0.30,
        trifurcation_center_frac: 1.0 / 3.0,
        cif_pretri_center_frac: 0.54,
        cif_terminal_tri_center_frac: 0.50,
        cif_terminal_bi_treat_frac: 0.84,
        asymmetric_narrow_frac: 0.5,
        trifurcation_left_frac: 1.0 / 3.0,
        cross_section_shape: CrossSectionShape::Rectangular,
        treatment_zone_mode: TreatmentZoneMode::VenturiThroats,
        centerline_venturi_throat_count: 2,
    }
}

#[test]
fn save_schematic_svg_adds_throat_markers_for_venturi_design() {
    let candidate = candidate_with_topology("venturi_test", DesignTopology::SingleVenturi);
    let out = unique_svg_path("cfd_optim_venturi");

    save_schematic_svg(&candidate, &out).expect("venturi schematic export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("IN"));
    assert!(svg.contains("OUT"));
    assert!(svg.contains("TH1"));
}

#[test]
fn save_schematic_svg_keeps_node_markers_without_throats_for_non_venturi_design() {
    let mut candidate = candidate_with_topology("nonventuri_test", DesignTopology::SerpentineGrid);
    candidate.treatment_zone_mode = TreatmentZoneMode::UltrasoundOnly;
    let out = unique_svg_path("cfd_optim_nonventuri");

    save_schematic_svg(&candidate, &out).expect("non-venturi schematic export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("IN"));
    assert!(svg.contains("OUT"));
    assert!(!svg.contains("TH1"));
}

#[test]
fn save_schematic_svg_renders_selective_blueprint_topology_from_blueprint_ssot() {
    let candidate = candidate_with_topology(
        "selective_cif",
        DesignTopology::PrimitiveSelectiveTree {
            sequence: PrimitiveSplitSequence::TriBi,
        },
    );
    let out = unique_svg_path("cfd_optim_selective_cif");

    save_schematic_svg(&candidate, &out).expect("selective schematic export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("TH1"));
    assert!(svg.contains("S1"));
    assert!(svg.contains("layers"));
}

#[test]
fn save_schematic_svg_renders_full_tree_split_markers_for_dtcv() {
    let candidate = candidate_with_topology(
        "selective_dtcv",
        DesignTopology::PrimitiveSelectiveTree {
            sequence: PrimitiveSplitSequence::TriTri,
        },
    );
    let out = unique_svg_path("cfd_optim_selective_dtcv");

    save_schematic_svg(&candidate, &out).expect("dtcv schematic export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("TH1"));
    assert!(svg.contains("S4"));
}

#[test]
fn save_schematic_svg_renders_acoustic_dtcv_without_throats() {
    let mut candidate = candidate_with_topology(
        "selective_dtcv_acoustic",
        DesignTopology::PrimitiveSelectiveTree {
            sequence: PrimitiveSplitSequence::TriTri,
        },
    );
    candidate.treatment_zone_mode = TreatmentZoneMode::UltrasoundOnly;
    let out = unique_svg_path("cfd_optim_selective_dtcv_acoustic");

    save_schematic_svg(&candidate, &out).expect("acoustic dtcv schematic export must succeed");

    let svg = std::fs::read_to_string(&out).expect("must read rendered svg");
    assert!(svg.contains("S4"));
    assert!(!svg.contains("TH1"));
}
