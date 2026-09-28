//! Canonical preset constructors for
//! [`BlueprintTopologySpec`](crate::topology::BlueprintTopologySpec).
//!
//! Each constructor returns a declarative spec that can be fed to
//! [`BlueprintTopologyFactory::build()`](crate::topology::BlueprintTopologyFactory::build)
//! to generate the corresponding
//! [`NetworkBlueprint`](crate::domain::model::NetworkBlueprint). The GA
//! composes arbitrary topologies by chaining
//! [`BlueprintTopologyMutation`](crate::topology::BlueprintTopologyMutation)
//! operations on these seeds.
//!
//! ## Three-stage optimisation pipeline
//!
//! | Stage | Preset family | What varies |
//! |-------|---------------|------------|
//! | 1 — Residence/Separation | `asymmetric_split_tree_spec` | Split kinds, depths, per-branch widths |
//! | 2 — Venturi Cavitation | + `with_venturi_placements()` | Throat geometry, channel targeting |
//! | 3 — GA Refinement | + `with_serpentine()` / `with_dean_venturi()` | Serpentine insertion, Dean placement |

mod milestone12;
mod modifiers;
mod parallel;
mod plate_presets;
pub mod sequence;
mod series;
mod tree;

pub use milestone12::{
    Milestone12PrimitiveSelectiveSpec, Milestone12StageBranchSpec, Milestone12StageLayout,
    Milestone12TopologyRequest, build_milestone12_blueprint, build_milestone12_topology_spec,
    enumerate_milestone12_topologies, milestone12_default_stage_layouts,
    milestone12_primitive_selective_tree_spec, promote_milestone12_option1_to_option2,
};
pub use modifiers::{
    with_branch_serpentine, with_dean_venturi_placement, with_venturi, with_venturi_placements,
};
pub use parallel::{parallel_microchannel_array_spec, parallel_path_spec};
pub use sequence::{ALL_SELECTIVE_SEQUENCES, PrimitiveSplitSequence, TRI_FIRST_SEQUENCES};
pub use series::{
    constriction_expansion_series_spec, serial_double_venturi_series_spec, series_path_spec,
    serpentine_bend_venturi_series_spec, serpentine_series_spec, single_venturi_series_spec,
    spiral_serpentine_series_spec, venturi_serpentine_series_spec,
};
pub use tree::{asymmetric_split_tree_spec, symmetric_n_furcation_spec};

#[cfg(test)]
mod tests;
