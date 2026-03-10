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
mod modifiers;
mod parallel;
mod series;
mod tree;

pub use modifiers::{with_branch_serpentine, with_dean_venturi_placement, with_venturi_placements};
pub use parallel::{parallel_microchannel_array_spec, parallel_path_spec};
pub use series::{
    constriction_expansion_series_spec, serial_double_venturi_series_spec,
    serpentine_bend_venturi_series_spec, serpentine_series_spec, series_path_spec,
    single_venturi_series_spec, spiral_serpentine_series_spec, venturi_serpentine_series_spec,
};
pub use tree::{asymmetric_split_tree_spec, asymmetric_trifurcation_spec, symmetric_bifurcation_spec};

#[cfg(test)]
mod tests {
    use super::*;
    use crate::topology::model::{SplitKind, VenturiPlacementMode};

    #[test]
    fn symmetric_bifurcation_spec_produces_equal_branches() {
        let spec = symmetric_bifurcation_spec("test_bi", 2, 6e-3, 1e-3, 5e-3);
        assert_eq!(spec.split_stages.len(), 2);
        for stage in &spec.split_stages {
            assert_eq!(stage.split_kind, SplitKind::Bifurcation);
            assert_eq!(stage.branches.len(), 2);
        }
        let leaf_count = spec.terminal_branch_count();
        assert_eq!(leaf_count, 4, "Bi→Bi = 4 leaves");
    }

    #[test]
    fn asymmetric_trifurcation_spec_center_arm_is_treatment() {
        let spec = asymmetric_trifurcation_spec("test_tri", 1, 6e-3, 0.5, 1e-3, 5e-3);
        let stage = &spec.split_stages[0];
        assert_eq!(stage.branches.len(), 3);
        assert!(stage.branches[0].treatment_path);
        assert!(!stage.branches[1].treatment_path);
        assert!(!stage.branches[2].treatment_path);

        let center_w = stage.branches[0].route.width_m;
        let left_w = stage.branches[1].route.width_m;
        assert!(
            (center_w - 3e-3).abs() < 1e-9,
            "center = 50% of 6mm = 3mm, got {center_w}"
        );
        assert!(
            (left_w - 1.5e-3).abs() < 1e-9,
            "peripheral = 25% of 6mm = 1.5mm, got {left_w}"
        );
    }

    #[test]
    fn with_venturi_targets_treatment_channels() {
        let base = asymmetric_trifurcation_spec("test_v", 1, 6e-3, 0.4, 1e-3, 5e-3);
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
    fn with_branch_serpentine_modifies_treatment_branch() {
        let base = asymmetric_trifurcation_spec("test_serp", 1, 6e-3, 0.4, 1e-3, 5e-3);
        let spec = with_branch_serpentine(base, "stage_0", "center", 6, 4e-3, 10e-3)
            .expect("branch should be found");
        assert!(spec.has_serpentine());
    }
}
