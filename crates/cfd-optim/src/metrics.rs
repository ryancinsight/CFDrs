//! Physics-based metric computation for SDT millifluidic design candidates.
//!
//! Each [`DesignCandidate`] is evaluated against the 1-D physics models from
//! `cfd-1d` (venturi + serpentine) and the Casson blood rheology from
//! `cfd-core`.  The resulting [`SdtMetrics`] struct gathers all values needed
//! for multi-objective scoring.
//!
//! # Cavitation physics
//! The cavitation number is defined as:
//! ```text
//! σ = (P_inlet − P_vapor) / (½ ρ V_throat²)
//! ```
//! Cavitation onset occurs when **σ < 1**.  In practical millifluidic devices
//! driven at 1–3 bar gauge by syringe pumps, σ < 1 is achievable only with
//! very small throat diameters (50–150 μm) at the higher flow rates (5–10 mL/min).
//!
//! # Haemolysis model
//! The combined haemolysis index (HI) per device pass is computed via the
//! Giersiepen (1990) model summed over each flow segment:
//! ```text
//! ΔHb/Hb = C · t^α · τ_wall^β
//! ```
//! where *t* is the transit time through a segment and *τ_wall* is the local
//! wall shear stress.  A per-pass HI < 0.001 (0.1 %) is the acceptability
//! threshold used here.

use cfd_1d::{FlowConditions, SerpentineModel, VenturiModel};
use cfd_core::physics::fluid::blood::CassonBlood;
use serde::{Deserialize, Serialize};

use crate::constraints::*;
use crate::design::{DesignCandidate, DesignTopology};
use crate::error::OptimError;

// ── SDT metrics struct ───────────────────────────────────────────────────────

/// All physics-derived metrics for one [`DesignCandidate`].
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SdtMetrics {
    // ── Cavitation ──
    /// Cavitation number σ at the primary venturi throat.
    /// `< 1` implies active cavitation inception; lower values = more intense.
    /// `f64::INFINITY` when the topology has no venturi.
    pub cavitation_number: f64,

    /// Dimensionless cavitation potential: `max(0, 1 − σ)`.
    /// `0` = no cavitation; approaching `1` = very strong cavitation.
    pub cavitation_potential: f64,

    /// Wall shear rate at the venturi throat [1/s].
    pub throat_shear_rate_inv_s: f64,

    /// Estimated wall shear stress at the venturi throat [Pa].
    /// May greatly exceed the FDA 150 Pa limit; transit time is < 1 μs–1 ms.
    pub throat_shear_pa: f64,

    /// Whether the throat shear stress exceeds FDA 150 Pa.
    /// Expected `true` for effective SDT; cumulative haemolysis is tracked
    /// separately via [`hemolysis_index_per_pass`].
    pub throat_exceeds_fda: bool,

    // ── Main-channel safety ──
    /// Maximum wall shear stress in sustained-contact (non-throat) channels [Pa].
    /// Must be < 150 Pa for FDA compliance.
    pub max_main_channel_shear_pa: f64,

    /// FDA compliance of main distribution / serpentine channels.
    /// `true` iff `max_main_channel_shear_pa < 150 Pa`.
    pub fda_main_compliant: bool,

    // ── Cumulative haemolysis ──
    /// Giersiepen (1990) haemolysis index per device pass (dimensionless fraction).
    /// Sum of throat + main-channel contributions.
    pub hemolysis_index_per_pass: f64,

    // ── Distribution / exposure ──
    /// Flow uniformity across parallel outlet branches: `1 − CV(Q_outlets)`.
    /// `1.0` = perfectly uniform.  Symmetric designs reach 1.0 by construction.
    pub flow_uniformity: f64,

    /// Fraction of the 36 treatment wells covered by channel passes.
    /// `1.0` = all 36 wells are in the channel path.
    pub well_coverage_fraction: f64,

    /// Mean fluid residence time in the treatment zone [s].
    pub mean_residence_time_s: f64,

    // ── System-level ──
    /// Total pressure drop across the device [Pa].
    pub total_pressure_drop_pa: f64,

    /// Total channel path length [mm].
    pub total_path_length_mm: f64,

    /// `true` if the total ΔP < available gauge pressure (syringe pump can
    /// drive the flow without stalling).
    pub pressure_feasible: bool,
}

// ── Public entry point ───────────────────────────────────────────────────────

/// Evaluate all physics-based metrics for a single design candidate.
///
/// Uses [`CassonBlood::normal_blood`] (Casson non-Newtonian blood at 37 °C)
/// and the `cfd-1d` models for venturi and serpentine stages.
///
/// # Errors
/// Returns [`OptimError::PhysicsError`] if a 1-D model fails (e.g. zero
/// throat diameter passed to a venturi-bearing topology).
pub fn compute_metrics(candidate: &DesignCandidate) -> Result<SdtMetrics, OptimError> {
    let blood = CassonBlood::<f64>::normal_blood();

    // ── Venturi stage ────────────────────────────────────────────────────────
    let (
        cavitation_number,
        cavitation_potential,
        throat_shear_rate,
        throat_shear_pa,
        dp_venturi,
        path_venturi_mm,
        v_throat,
    ) = if candidate.topology.has_venturi() {
        eval_venturi_stage(candidate, &blood)?
    } else {
        (f64::INFINITY, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0)
    };

    let throat_exceeds_fda = throat_shear_pa > FDA_MAX_WALL_SHEAR_PA;

    // ── Distribution / serpentine stage ─────────────────────────────────────
    let DistributionResult {
        flow_uniformity,
        well_coverage_fraction,
        max_main_shear_pa,
        dp_distribution_pa,
        path_distribution_mm,
        residence_time_s,
    } = eval_distribution_stage(candidate, &blood)?;

    let fda_main_compliant = max_main_shear_pa <= FDA_MAX_WALL_SHEAR_PA;

    // ── Cumulative haemolysis index ─────────────────────────────────────────
    let hemolysis_index_per_pass = {
        // Throat contribution (brief high-shear transit)
        let hi_throat = if candidate.topology.has_venturi() && v_throat > 0.0 {
            let t_throat = candidate.throat_length_m / v_throat; // transit time [s]
            giersiepen_hi(throat_shear_pa, t_throat)
        } else {
            0.0
        };

        // Main-channel contribution (long low-shear exposure)
        let hi_main = if residence_time_s > 0.0 {
            giersiepen_hi(max_main_shear_pa, residence_time_s)
        } else {
            0.0
        };

        hi_throat + hi_main
    };

    // ── Totals ───────────────────────────────────────────────────────────────
    let total_dp = dp_venturi + dp_distribution_pa;
    let total_path = path_venturi_mm + path_distribution_mm;
    let pressure_feasible = total_dp <= candidate.inlet_gauge_pa;

    Ok(SdtMetrics {
        cavitation_number,
        cavitation_potential,
        throat_shear_rate_inv_s: throat_shear_rate,
        throat_shear_pa,
        throat_exceeds_fda,
        max_main_channel_shear_pa: max_main_shear_pa,
        fda_main_compliant,
        hemolysis_index_per_pass,
        flow_uniformity,
        well_coverage_fraction,
        mean_residence_time_s: residence_time_s,
        total_pressure_drop_pa: total_dp,
        total_path_length_mm: total_path,
        pressure_feasible,
    })
}

// ── Venturi stage evaluator ──────────────────────────────────────────────────

/// Returns `(σ, cav_potential, shear_rate_1/s, shear_pa, dp_pa, path_mm, v_throat)`.
fn eval_venturi_stage(
    candidate: &DesignCandidate,
    blood: &CassonBlood<f64>,
) -> Result<(f64, f64, f64, f64, f64, f64, f64), OptimError> {
    let d_inlet = candidate.inlet_diameter_m;
    let d_throat = candidate.throat_diameter_m;
    let l_throat = candidate.throat_length_m.max(d_throat * 2.0); // ensure non-trivial length

    if d_throat <= 0.0 {
        return Err(OptimError::InvalidParameter(format!(
            "candidate '{}': throat_diameter_m must be > 0 for topology '{}'",
            candidate.id,
            candidate.topology.name()
        )));
    }

    // Per-venturi flow rate and inlet velocity
    let q_v = candidate.per_venturi_flow();
    let a_inlet = std::f64::consts::FRAC_PI_4 * d_inlet * d_inlet;
    let v_inlet = q_v / a_inlet;

    let model = VenturiModel::<f64>::millifluidic(d_inlet, d_throat, l_throat);
    let mut cond = FlowConditions::<f64>::new(v_inlet);
    cond.pressure = candidate.inlet_pressure_pa();

    let analysis = model.analyze(blood, &cond).map_err(|e| OptimError::PhysicsError {
        id: candidate.id.clone(),
        reason: e.to_string(),
    })?;

    let v_throat = analysis.throat_velocity;
    let shear_rate = analysis.throat_shear_rate;
    let mu_eff = blood.apparent_viscosity(shear_rate);
    let shear_pa = mu_eff * shear_rate;

    // Cavitation number: σ = (P_inlet − P_vapor) / (½ρV²)
    let sigma = if v_throat > 0.0 {
        let dynamic_pressure = 0.5 * BLOOD_DENSITY_KG_M3 * v_throat * v_throat;
        (candidate.inlet_pressure_pa() - BLOOD_VAPOR_PRESSURE_PA) / dynamic_pressure
    } else {
        f64::INFINITY
    };

    let cav_potential = if sigma < SIGMA_CRIT && sigma >= 0.0 {
        (1.0 - sigma / SIGMA_CRIT).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Approximate venturi path length: inlet cone (4×d) + throat (l_throat) + diffuser (6×d)
    let path_mm = (d_inlet * 10.0 + l_throat) * candidate.topology.venturi_count() as f64 * 1000.0;

    // Total pressure drop across ALL parallel venturis (same ΔP, different Q → same ΔP)
    let dp = analysis.dp_total;

    Ok((sigma, cav_potential, shear_rate, shear_pa, dp, path_mm, v_throat))
}

// ── Distribution stage evaluator ────────────────────────────────────────────

struct DistributionResult {
    flow_uniformity: f64,
    well_coverage_fraction: f64,
    max_main_shear_pa: f64,
    dp_distribution_pa: f64,
    path_distribution_mm: f64,
    residence_time_s: f64,
}

fn eval_distribution_stage(
    candidate: &DesignCandidate,
    blood: &CassonBlood<f64>,
) -> Result<DistributionResult, OptimError> {
    match candidate.topology {
        // ── Serpentine-based topologies ──────────────────────────────────
        DesignTopology::SerpentineGrid | DesignTopology::VenturiSerpentine => {
            eval_serpentine(candidate, blood)
        }

        // ── Bifurcation / trifurcation ───────────────────────────────────
        DesignTopology::BifurcationVenturi | DesignTopology::TrifurcationVenturi => {
            eval_branching(candidate, blood)
        }

        // ── Single venturi with no distribution stage ────────────────────
        DesignTopology::SingleVenturi => {
            // Wall shear in the short outlet channel downstream of the venturi
            // (approximated as main channel shear at inlet velocity)
            let q = candidate.flow_rate_m3_s;
            let w = candidate.channel_width_m;
            let h = candidate.channel_height_m;
            let v = q / (w * h);
            let shear_rate = 6.0 * v / h;
            let mu_eff = blood.apparent_viscosity(shear_rate);
            let shear_pa = mu_eff * shear_rate;

            Ok(DistributionResult {
                flow_uniformity: 1.0,
                // Single outlet → 1 cavitation site, coverage ≈ 1 well / 36
                well_coverage_fraction: DesignTopology::SingleVenturi.nominal_well_coverage(),
                max_main_shear_pa: shear_pa,
                dp_distribution_pa: 0.0, // negligible outlet stub
                path_distribution_mm: 0.0,
                residence_time_s: 0.0,
            })
        }
    }
}

/// Evaluate the serpentine distribution stage.
fn eval_serpentine(
    candidate: &DesignCandidate,
    blood: &CassonBlood<f64>,
) -> Result<DistributionResult, OptimError> {
    let w = candidate.channel_width_m;
    let h = candidate.channel_height_m;
    let seg_len = candidate.segment_length_m;
    let n_segs = candidate.serpentine_segments;
    let bend_r = candidate.bend_radius_m;

    if n_segs == 0 {
        return Err(OptimError::InvalidParameter(format!(
            "candidate '{}': serpentine_segments must be >= 1",
            candidate.id
        )));
    }

    // For VenturiSerpentine the full flow enters the serpentine downstream of
    // the single venturi; for SerpentineGrid the full flow enters directly.
    let q_ch = candidate.flow_rate_m3_s;
    let v_ch = q_ch / (w * h);

    let model =
        SerpentineModel::<f64>::millifluidic_rectangular(w, h, seg_len, n_segs, bend_r);
    let cond = FlowConditions::<f64>::new(v_ch);

    let analysis = model.analyze(blood, &cond).map_err(|e| OptimError::PhysicsError {
        id: candidate.id.clone(),
        reason: e.to_string(),
    })?;

    let shear_rate = analysis.wall_shear_rate;
    let mu_eff = blood.apparent_viscosity(shear_rate);
    let shear_pa = mu_eff * shear_rate;

    // Path length: straights + semicircular bends
    let n_bends = (n_segs.saturating_sub(1)) as f64;
    let path_m =
        (n_segs as f64) * seg_len + n_bends * std::f64::consts::PI * bend_r;
    let path_mm = path_m * 1000.0;

    // Well coverage: normalise path length against the minimum for full coverage
    // (6 row-sweeps × 45 mm = 270 mm).
    let coverage = (path_m / FULL_GRID_SERPENTINE_LENGTH_M).min(1.0);

    // Residence time = channel volume / flow rate
    let vol = path_m * w * h;
    let residence = vol / q_ch;

    // Single serpentine → perfect flow uniformity (single stream)
    Ok(DistributionResult {
        flow_uniformity: 1.0,
        well_coverage_fraction: coverage,
        max_main_shear_pa: shear_pa,
        dp_distribution_pa: analysis.dp_total,
        path_distribution_mm: path_mm,
        residence_time_s: residence,
    })
}

/// Evaluate the distribution channels in branching topologies
/// (bifurcation + trifurcation).
///
/// Uses fully-developed rectangular-channel shear estimate for the main
/// trunk and daughter channels.
fn eval_branching(
    candidate: &DesignCandidate,
    blood: &CassonBlood<f64>,
) -> Result<DistributionResult, OptimError> {
    let n_branches = candidate.topology.outlet_count() as f64;
    let w = candidate.channel_width_m;
    let h = candidate.channel_height_m;

    // Main trunk: carries full flow Q
    let q_trunk = candidate.flow_rate_m3_s;
    let shear_trunk = compute_rect_wall_shear(q_trunk, w, h, blood);

    // Daughter branches: each carries Q / n_branches
    let q_branch = q_trunk / n_branches;
    let shear_branch = compute_rect_wall_shear(q_branch, w, h, blood);

    let max_shear = shear_trunk.max(shear_branch);

    // Approximate ΔP of the distribution network using Hagen-Poiseuille scaling:
    //   ΔP ≈ 12 μ L Q / (h³ w)   [rectangular channel]
    // Use average viscosity from Casson model at the dominant shear rate.
    let mu_trunk = blood.apparent_viscosity(6.0 * (q_trunk / (w * h)) / h);
    // Branch length ≈ half-pitch from junction to well row
    let l_branch = TREATMENT_HEIGHT_MM * 0.5e-3; // ~22.5 mm
    let dh = 2.0 * w * h / (w + h); // hydraulic diameter
    let re_branch = BLOOD_DENSITY_KG_M3 * (q_branch / (w * h)) * dh / mu_trunk;
    let f = if re_branch > 0.0 { 64.0 / re_branch } else { 64.0 };
    let v_branch = q_branch / (w * h);
    let dp_branch = f * (l_branch / dh) * 0.5 * BLOOD_DENSITY_KG_M3 * v_branch * v_branch;

    // Path length = inlet trunk + branch channels
    let trunk_len_mm = TREATMENT_HEIGHT_MM * 0.5; // 22.5 mm to bifurcation point
    let branch_len_mm = l_branch * 1000.0 * n_branches;
    let path_mm = trunk_len_mm + branch_len_mm;

    // Residence time in distribution channels
    let vol_trunk = (trunk_len_mm * 1e-3) * w * h;
    let vol_branches = (l_branch * n_branches) * w * h;
    let residence = (vol_trunk + vol_branches) / q_trunk;

    // Symmetric designs: perfect flow uniformity
    Ok(DistributionResult {
        flow_uniformity: 1.0,
        well_coverage_fraction: candidate.topology.nominal_well_coverage(),
        max_main_shear_pa: max_shear,
        dp_distribution_pa: dp_branch,
        path_distribution_mm: path_mm,
        residence_time_s: residence,
    })
}

// ── Physics helpers ──────────────────────────────────────────────────────────

/// Giersiepen (1990) haemolysis index contribution:
/// `HI = C · t^α · τ^β`
///
/// Returns 0 for non-positive inputs (no damage).
#[inline]
pub fn giersiepen_hi(tau_pa: f64, t_s: f64) -> f64 {
    if tau_pa <= 0.0 || t_s <= 0.0 {
        return 0.0;
    }
    GIERSIEPEN_C * t_s.powf(GIERSIEPEN_ALPHA) * tau_pa.powf(GIERSIEPEN_BETA)
}

/// Wall shear stress in a rectangular channel (fully developed Poiseuille flow).
///
/// Approximation: τ_wall ≈ μ_eff × 6V/h  (valid when w >> h).
fn compute_rect_wall_shear(q: f64, w: f64, h: f64, blood: &CassonBlood<f64>) -> f64 {
    let v = q / (w * h);
    let gamma = 6.0 * v / h;
    let mu = blood.apparent_viscosity(gamma);
    mu * gamma
}
