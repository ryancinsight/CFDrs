//! Canonical preset constructors for [`BlueprintTopologySpec`].
//!
//! Each constructor returns a declarative spec that can be fed to
//! [`BlueprintTopologyFactory::build()`] to generate the corresponding
//! [`NetworkBlueprint`].  The GA composes arbitrary topologies by chaining
//! [`BlueprintTopologyMutation`] operations on these seeds.
//!
//! ## Three-stage optimisation pipeline
//!
//! | Stage | Preset family | What varies |
//! |-------|---------------|------------|
//! | 1 — Residence/Separation | `asymmetric_split_tree_spec` | Split kinds, depths, per-branch widths |
//! | 2 — Venturi Cavitation | + `with_venturi_placements()` | Throat geometry, channel targeting |
//! | 3 — GA Refinement | + `with_serpentine()` / `with_dean_venturi()` | Serpentine insertion, Dean placement |

mod helpers;
mod milestone12;
mod modifiers;
mod parallel;
pub mod sequence;
mod series;
mod tree;

pub use milestone12::{
    build_milestone12_blueprint, build_milestone12_topology_spec, enumerate_milestone12_topologies,
    milestone12_default_stage_layouts, milestone12_primitive_selective_tree_spec,
    promote_milestone12_option1_to_option2, Milestone12PrimitiveSelectiveSpec,
    Milestone12StageBranchSpec, Milestone12StageLayout, Milestone12TopologyRequest,
};
pub use modifiers::{
    with_branch_serpentine, with_dean_venturi_placement, with_venturi, with_venturi_placements,
};
pub use parallel::{parallel_microchannel_array_spec, parallel_path_spec};
pub use sequence::{PrimitiveSplitSequence, ALL_SELECTIVE_SEQUENCES, TRI_FIRST_SEQUENCES};
pub use series::{
    constriction_expansion_series_spec, serial_double_venturi_series_spec, series_path_spec,
    serpentine_bend_venturi_series_spec, serpentine_series_spec, single_venturi_series_spec,
    spiral_serpentine_series_spec, venturi_serpentine_series_spec,
};
pub use tree::{asymmetric_split_tree_spec, symmetric_n_furcation_spec};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::model::{SplitKind, VenturiConfig, VenturiPlacementMode};

    #[test]
    fn symmetric_n_furcation_spec_produces_equal_branches() {
        let spec = symmetric_n_furcation_spec("test_bi", 2, 2, 6e-3, 1e-3, 5e-3);
        assert_eq!(spec.split_stages.len(), 2);
        for stage in &spec.split_stages {
            assert_eq!(stage.split_kind, SplitKind::NFurcation(2));
            assert_eq!(stage.branches.len(), 2);
        }
        let leaf_count = spec.terminal_branch_count();
        assert_eq!(leaf_count, 4, "Bi→Bi = 4 leaves");
    }

    #[test]
    fn symmetric_n_furcation_spec_produces_odd_branches() {
        let spec = symmetric_n_furcation_spec("test_tri", 1, 3, 6e-3, 1e-3, 5e-3);
        let stage = &spec.split_stages[0];
        assert_eq!(stage.branches.len(), 3);
        assert!(!stage.branches[0].treatment_path);
        assert!(stage.branches[1].treatment_path);
        assert!(!stage.branches[2].treatment_path);

        let center_w = stage.branches[1].route.width_m;
        let left_w = stage.branches[0].route.width_m;
        assert!(
            (center_w - 2e-3).abs() < 1e-9,
            "center = 33.3% of 6mm = 2mm, got {center_w}"
        );
        assert!(
            (left_w - 2e-3).abs() < 1e-9,
            "peripheral = 33.3% of 6mm = 2mm, got {left_w}"
        );
    }

    #[test]
    fn with_venturi_targets_treatment_channels() {
        let base = symmetric_n_furcation_spec("test_v", 1, 3, 6e-3, 1e-3, 5e-3);
        let spec = with_venturi_placements(
            base,
            50e-6,
            1e-3,
            250e-6,
            1,
            VenturiPlacementMode::StraightSegment,
        );
        assert!(spec.has_venturi());
        assert!(spec.uses_venturi_treatment());
        assert_eq!(spec.venturi_placements.len(), 1);
    }

    #[test]
    fn with_venturi_supports_explicit_target_list() {
        let base = parallel_path_spec(
            "parallel_v",
            2.0e-3,
            2.0e-3,
            8.0e-3,
            8.0e-3,
            vec![
                crate::topology::ParallelChannelSpec {
                    channel_id: "lane_a".to_string(),
                    route: crate::topology::ChannelRouteSpec {
                        length_m: 10.0e-3,
                        width_m: 1.0e-3,
                        height_m: 1.0e-3,
                        serpentine: None,
                        therapy_zone: crate::domain::therapy_metadata::TherapyZone::CancerTarget,
                    },
                },
                crate::topology::ParallelChannelSpec {
                    channel_id: "lane_b".to_string(),
                    route: crate::topology::ChannelRouteSpec {
                        length_m: 10.0e-3,
                        width_m: 1.2e-3,
                        height_m: 1.0e-3,
                        serpentine: None,
                        therapy_zone: crate::domain::therapy_metadata::TherapyZone::CancerTarget,
                    },
                },
            ],
            crate::topology::TreatmentActuationMode::UltrasoundOnly,
        );
        let spec = with_venturi(
            base,
            VenturiConfig {
                target_channel_ids: vec!["lane_b".to_string()],
                serial_throat_count: 1,
                throat_geometry: crate::topology::ThroatGeometrySpec {
                    throat_width_m: 60e-6,
                    throat_height_m: 1e-3,
                    throat_length_m: 200e-6,
                    inlet_width_m: 0.0,
                    outlet_width_m: 0.0,
                    convergent_half_angle_deg: 0.0,
                    divergent_half_angle_deg: 0.0,
                },
                placement_mode: VenturiPlacementMode::StraightSegment,
            },
        );
        assert_eq!(spec.venturi_placements.len(), 1);
        assert_eq!(spec.venturi_placements[0].target_channel_id, "lane_b");
    }

    #[test]
    fn with_branch_serpentine_modifies_treatment_branch() {
        let base = symmetric_n_furcation_spec("test_serp", 1, 3, 6e-3, 1e-3, 5e-3);
        let spec = with_branch_serpentine(base, "stage_0", "center", 6, 4e-3, 10e-3)
            .expect("branch should be found");
        assert!(spec.has_serpentine());
    }
}
