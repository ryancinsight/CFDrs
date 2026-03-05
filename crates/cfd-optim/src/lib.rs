//! # cfd-optim — Millifluidic SDT Design Optimiser
//!
//! Generates and ranks the **top-5 millifluidic designs** for Sonodynamic
//! Therapy (SDT) using either:
//!
//! - **Hydrodynamic cavitation** through venturi throats — collapses microbubbles
//!   to activate sonosensitisers in targeted tissue wells.
//! - **Uniform light / ultrasound exposure** across a 6 × 6 well-plate centre
//!   zone — ensures every well receives equal treatment dose.
//!
//! All designs are constrained to fit within a standard **96-well plate footprint**
//! (ANSI/SLAS 1-2004, 127.76 × 85.47 mm) and comply with the FDA haemolysis
//! guideline for blood-contacting devices (≤ 150 Pa sustained wall shear).
//!
//! ## Quick start
//!
//! ```rust,ignore
//! use cfd_optim::{SdtOptimizer, OptimMode, SdtWeights};
//!
//! // SDT cavitation mode — top 5 designs
//! let top5 = SdtOptimizer::new(OptimMode::SdtCavitation, SdtWeights::default())
//!     .top_5()
//!     .expect("optimisation failed");
//!
//! for d in &top5 {
//!     println!("{}", d.summary());
//! }
//! ```
//!
//! ## Module structure
//!
//! | Module | Purpose |
//! |--------|---------|
//! [`constraints`]    | 96-well plate geometry, blood properties, FDA limits |
//! [`design`]         | Topology families, candidate struct, parameter sweep |
//! [`metrics`]        | Physics-based metric computation (cfd-1d + Casson blood) |
//! [`scoring`]        | Multi-objective scoring functions |
//! [`orchestration`]  | Parametric sweep + ranking engine (`SdtOptimizer`) |
//! [`evo`]            | Genetic algorithm (`GeneticOptimizer`) |
//! [`delivery`]       | JSON/SVG export and mesh pipeline |
//! [`reporting`]      | Deterministic ranking helpers and canonical markdown writer |
//! [`error`]          | Error types |

#![warn(clippy::all)]
#![warn(clippy::pedantic)]
// Numerical cast patterns common in optimisation scoring.
#![allow(clippy::cast_precision_loss)]
#![allow(clippy::cast_possible_truncation)]
#![allow(clippy::cast_sign_loss)]
// Physics variable names (x, y, k, n, …) are conventional.
#![allow(clippy::similar_names)]
#![allow(clippy::many_single_char_names)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::module_name_repetitions)]
#![allow(clippy::missing_errors_doc)]
#![allow(clippy::missing_panics_doc)]
#![allow(clippy::must_use_candidate)]
#![allow(clippy::redundant_closure_for_method_calls)]
#![allow(clippy::doc_markdown)]
#![allow(clippy::unreadable_literal)]
#![allow(clippy::manual_let_else)]
#![allow(clippy::unnecessary_wraps)]
#![allow(clippy::match_same_arms)]
#![allow(clippy::useless_conversion)]
#![allow(clippy::inline_always)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::needless_pass_by_value)]
#![allow(clippy::items_after_statements)]
#![allow(clippy::format_push_string)]
#![allow(clippy::field_reassign_with_default)]
#![allow(clippy::new_without_default)]
#![allow(clippy::trivially_copy_pass_by_ref)]
#![allow(clippy::float_cmp)]
#![allow(clippy::return_self_not_must_use)]
#![allow(clippy::cast_possible_wrap)]
#![allow(clippy::too_many_arguments)]
#![allow(clippy::unused_self)]
#![allow(clippy::ptr_arg)]
#![allow(clippy::struct_excessive_bools)]

// ── SDT-specific modules ─────────────────────────────────────────────────────

pub mod analysis;
pub mod constraints;
pub mod delivery;
pub mod design;
pub mod error;
pub mod evo;
pub mod metrics;
pub mod metrics_cache;
pub mod orchestration;
pub mod reporting;
pub mod scoring;

// ── Top-level re-exports ─────────────────────────────────────────────────────

pub use delivery::{
    save_all_modes_json, save_annotated_cif_svg, save_comparison_svg, save_schematic_svg,
    save_top5_json,
};
pub use design::{
    build_candidate_space, sample_random_candidates, CrossSectionShape, DesignCandidate,
    DesignTopology, TreatmentZoneMode,
};
pub use error::OptimError;
pub use evo::{candidate_to_genome, decode_genome, GeneticOptimizer, MillifluidicGenome};
pub use metrics::{compute_metrics, giersiepen_hi, ChannelHemolysis, SdtMetrics};
pub use metrics_cache::{MetricsCache, MetricsCacheStats, METRICS_CACHE_VERSION};
pub use orchestration::{
    OptimStats, RankedDesign, RobustScoreStats, RobustSweepConfig, SdtOptimizer,
};
pub use scoring::{score_candidate, score_description, OptimMode, SdtWeights};

#[cfg(feature = "mesh-export")]
pub use delivery::{DesignArtifacts, DesignPipeline};

// ── Legacy API (backward compatibility) ──────────────────────────────────────
//
// The types below provide backward-compatible API surface so that any
// external code that references them continues to compile without changes.

use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ScenarioResult {
    pub id: String,
    pub scheme_json_path: String,
    pub flow_uniformity: f64,
    pub pressure_drop_pa: f64,
    pub max_shear_pa: f64,
    pub throughput_ml_min: f64,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct ObjectiveWeights {
    pub flow_uniformity_weight: f64,
    pub pressure_drop_weight: f64,
    pub max_shear_weight: f64,
    pub throughput_weight: f64,
}

impl Default for ObjectiveWeights {
    fn default() -> Self {
        Self {
            flow_uniformity_weight: 0.40,
            pressure_drop_weight: 0.20,
            max_shear_weight: 0.20,
            throughput_weight: 0.20,
        }
    }
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct RankedScenario {
    pub scenario: ScenarioResult,
    pub score: f64,
}

// Legacy error variants are handled by the re-exported `OptimError` from `error.rs`.

pub fn optimize_scenarios(
    scenarios: &[ScenarioResult],
    weights: ObjectiveWeights,
) -> Result<Vec<RankedScenario>, OptimError> {
    if scenarios.is_empty() {
        return Err(OptimError::EmptyInput);
    }

    let max_dp = scenarios
        .iter()
        .map(|s| s.pressure_drop_pa)
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let max_shear = scenarios
        .iter()
        .map(|s| s.max_shear_pa)
        .fold(0.0_f64, f64::max)
        .max(1.0);
    let max_throughput = scenarios
        .iter()
        .map(|s| s.throughput_ml_min)
        .fold(0.0_f64, f64::max)
        .max(1.0);

    let mut ranked: Vec<RankedScenario> = scenarios
        .iter()
        .cloned()
        .map(|scenario| {
            let dp_term = 1.0 - (scenario.pressure_drop_pa / max_dp).clamp(0.0, 1.0);
            let shear_term = 1.0 - (scenario.max_shear_pa / max_shear).clamp(0.0, 1.0);
            let q_term = (scenario.throughput_ml_min / max_throughput).clamp(0.0, 1.0);
            let uniformity_term = scenario.flow_uniformity.clamp(0.0, 1.0);

            let score = weights.flow_uniformity_weight * uniformity_term
                + weights.pressure_drop_weight * dp_term
                + weights.max_shear_weight * shear_term
                + weights.throughput_weight * q_term;

            RankedScenario { scenario, score }
        })
        .collect();

    ranked.sort_by(|a, b| b.score.total_cmp(&a.score));
    Ok(ranked)
}

pub fn top_candidates(
    ranked: &[RankedScenario],
    top_k: usize,
) -> Result<Vec<RankedScenario>, OptimError> {
    if top_k == 0 {
        return Err(OptimError::InvalidTopK(top_k));
    }

    Ok(ranked.iter().take(top_k).cloned().collect())
}

#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshHandoffRecord {
    pub scenario_id: String,
    pub scheme_json_path: String,
    pub mesh_output_path: String,
}

#[cfg(feature = "mesh-export")]
pub fn generate_mesh_handoff(
    ranked: &[RankedScenario],
    top_k: usize,
    mesh_output_dir: &std::path::Path,
) -> Result<Vec<MeshHandoffRecord>, OptimError> {
    let _ = mesh_output_dir; // suppress unused-variable warning

    let selected = top_candidates(ranked, top_k)?;
    let records = Vec::with_capacity(selected.len());

    for item in selected {
        let scenario = item.scenario;
        let design_dir = mesh_output_dir.join(&scenario.id);
        std::fs::create_dir_all(&design_dir).map_err(|err| OptimError::MeshExportFailed {
            scenario_id: scenario.id.clone(),
            message: err.to_string(),
        })?;

        // Parse the interchange JSON to reconstruct the channel system, then
        // build a blueprint from it.  NetworkBlueprint does not implement
        // serde::Deserialize, so we parse as InterchangeChannelSystem instead
        // and then fall back to an error (full round-trip is not yet supported).
        let json_str = std::fs::read_to_string(&scenario.scheme_json_path).map_err(|err| {
            OptimError::MeshExportFailed {
                scenario_id: scenario.id.clone(),
                message: format!("Failed to read scheme JSON: {err}"),
            }
        })?;
        let interchange: cfd_schematics::geometry::InterchangeChannelSystem =
            serde_json::from_str(&json_str).map_err(|err| OptimError::MeshExportFailed {
                scenario_id: scenario.id.clone(),
                message: format!("Failed to parse interchange JSON: {err}"),
            })?;

        // InterchangeChannelSystem → NetworkBlueprint round-trip is not yet
        // implemented upstream. For now, return an informative error.
        return Err(OptimError::MeshExportFailed {
            scenario_id: scenario.id.clone(),
            message: format!(
                "Legacy generate_mesh_handoff does not support InterchangeChannelSystem → \
                 NetworkBlueprint conversion (schema v{}). Use DesignPipeline::export_design \
                 instead.",
                interchange.schema_version,
            ),
        });
    }

    Ok(records)
}

#[cfg(feature = "mesh-export")]
#[allow(dead_code)]
fn reassign_regions_from_labels(mesh: &mut cfd_mesh::domain::mesh::IndexedMesh) {
    use cfd_mesh::domain::core::index::{FaceId, RegionId};
    use std::collections::HashMap;

    let region_for: HashMap<FaceId, RegionId> = mesh
        .boundary_labels
        .iter()
        .map(|(&fid, label)| {
            let region = match label.as_str() {
                "inlet" => RegionId::new(0),
                "outlet" => RegionId::new(1),
                _ => RegionId::new(2),
            };
            (fid, region)
        })
        .collect();

    for (fid, face) in mesh.faces.iter_mut_enumerated() {
        face.region = *region_for.get(&fid).unwrap_or(&RegionId::new(2));
    }
}
