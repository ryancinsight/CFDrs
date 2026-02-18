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
//! ```rust,no_run
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
//! [`constraints`] | 96-well plate geometry, blood properties, FDA limits |
//! [`design`]      | Topology families, candidate struct, parameter sweep |
//! [`metrics`]     | Physics-based metric computation (cfd-1d + Casson blood) |
//! [`scoring`]     | Multi-objective scoring functions |
//! [`optimizer`]   | Parametric sweep + ranking engine |
//! [`error`]       | Error types |

// ── SDT-specific modules ─────────────────────────────────────────────────────

pub mod constraints;
pub mod design;
pub mod error;
pub mod metrics;
pub mod optimizer;
pub mod scoring;

// ── Top-level re-exports ─────────────────────────────────────────────────────

pub use design::{build_candidate_space, DesignCandidate, DesignTopology};
pub use error::OptimError;
pub use metrics::{compute_metrics, giersiepen_hi, SdtMetrics};
pub use optimizer::{OptimStats, RankedDesign, SdtOptimizer};
pub use scoring::{score_candidate, score_description, OptimMode, SdtWeights};

// ── Legacy API (backward compatibility) ──────────────────────────────────────
//
// The types below are preserved from the original generic stub so that any
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
    use blue2mesh::{export::CfdExporter, MeshGenerator};

    let selected = top_candidates(ranked, top_k)?;
    let mut records = Vec::with_capacity(selected.len());

    for item in selected {
        let scenario = item.scenario;
        let mesh = MeshGenerator::new()
            .generate_from_scheme(&scenario.scheme_json_path)
            .map_err(|err| OptimError::MeshExportFailed {
                scenario_id: scenario.id.clone(),
                message: err.to_string(),
            })?;

        let output = mesh_output_dir.join(format!("{}_mesh.vtk", scenario.id));
        CfdExporter::new()
            .export_vtk(&mesh, &output)
            .map_err(|err| OptimError::MeshExportFailed {
                scenario_id: scenario.id.clone(),
                message: err.to_string(),
            })?;

        records.push(MeshHandoffRecord {
            scenario_id: scenario.id,
            scheme_json_path: scenario.scheme_json_path,
            mesh_output_path: output.to_string_lossy().to_string(),
        });
    }

    Ok(records)
}
