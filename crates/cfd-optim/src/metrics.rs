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


use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_schematics::CrossSectionSpec;
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

    // ── Cell separation (2-population: cancer vs RBC) ──
    /// Inertial focusing separation efficiency for MCF-7 cancer cells vs RBCs.
    ///
    /// Defined as `|x̃_cancer − x̃_rbc|` ∈ [0, 1] where `x̃` is the normalised
    /// lateral equilibrium position (0 = center, 1 = wall).
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub cell_separation_efficiency: f64,

    /// Fraction of MCF-7 cancer cells collected in the center channel.
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub cancer_center_fraction: f64,

    /// Fraction of RBCs collected in the peripheral (wall) channels.
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub rbc_peripheral_fraction: f64,

    // ── Three-population separation (WBC+cancer→center, RBC→periphery) ──
    /// Three-population separation efficiency.
    ///
    /// Defined as `x̃_rbc − max(x̃_cancer, x̃_wbc)`, quantifying how well
    /// RBCs are pushed to the wall relative to both WBCs and cancer cells.
    /// Higher values mean both WBC and cancer stay central while RBCs migrate
    /// to the wall.
    ///
    /// `0.0` for topologies without a three-population separation stage.
    pub three_pop_sep_efficiency: f64,

    /// Fraction of WBCs collected in the center channel.
    ///
    /// `0.0` for topologies without a three-population separation stage.
    pub wbc_center_fraction: f64,

    /// Fraction of RBCs collected in the peripheral channels in the
    /// three-population separation topology.
    ///
    /// `0.0` for topologies without a three-population separation stage.
    pub rbc_peripheral_fraction_three_pop: f64,

    /// Normalised lateral equilibrium position of WBCs (0=center, 1=wall).
    ///
    /// `0.5` (dispersed) for topologies without a separation stage.
    pub wbc_equilibrium_pos: f64,

    /// Normalised lateral equilibrium position of cancer cells (0=center, 1=wall).
    ///
    /// `0.5` (dispersed) for topologies without a separation stage.
    pub cancer_equilibrium_pos: f64,

    /// Normalised lateral equilibrium position of RBCs (0=center, 1=wall).
    ///
    /// `0.5` (dispersed) for topologies without a separation stage.
    pub rbc_equilibrium_pos: f64,
}

// ── Public entry point ───────────────────────────────────────────────────────

/// Evaluate all physics-based metrics for a single design candidate.
///
/// All topologies in the design space are **symmetric** (equal pressure at all
/// parallel branches by construction), so flow splits equally among parallel
/// paths.  This allows direct analytical computation without a network solver:
///
/// 1. Trunk / distribution channels: rectangular duct Poiseuille (Shah-London).
/// 2. Venturi sections: [`VenturiModel::analyze`] with per-venturi flow.
/// 3. Serpentine sections: [`SerpentineModel::analyze`] with full flow.
/// 4. Cell separation: `CellSeparationModel` / `lateral_equilibrium`.
///
/// # Errors
/// Returns [`OptimError::PhysicsError`] if a 1-D model returns an error for
/// physically unrealisable parameters (e.g. throat velocity = 0).
pub fn compute_metrics(candidate: &DesignCandidate) -> Result<SdtMetrics, OptimError> {
    use cfd_1d::{FlowConditions, SerpentineModel};

    let blood = CassonBlood::<f64>::normal_blood();

    // ── Geometry from cfd-schematics blueprint ────────────────────────────────
    // Build the NetworkBlueprint via cfd-schematics presets so that channel
    // cross-section geometry is defined in the design library, not recomputed here.
    let bp = candidate.to_blueprint();

    // Main channel cross-section (inlet approach / serpentine / trunk).
    // Look for "inlet_section", "segment_1", or "parent" — whichever is first.
    let main_cs = bp.channels.iter()
        .find(|c| {
            let id = c.id.as_str();
            id == "inlet_section" || id == "parent" || id.starts_with("segment_")
        })
        .map(|c| c.cross_section)
        .unwrap_or(CrossSectionSpec::Rectangular {
            width_m:  candidate.channel_width_m,
            height_m: candidate.channel_height_m,
        });

    let (w, h) = match main_cs {
        CrossSectionSpec::Rectangular { width_m, height_m } => (width_m, height_m),
        CrossSectionSpec::Circular { diameter_m } => (diameter_m, diameter_m),
    };

    // Throat cross-section for venturi topologies.
    let throat_cs = bp.channels.iter()
        .find(|c| c.id.as_str() == "throat_section")
        .map(|c| c.cross_section);

    let q   = candidate.flow_rate_m3_s;
    let n_pv  = candidate.topology.parallel_venturi_count();
    let n_out = candidate.topology.outlet_count().max(1);
    let q_per_venturi = if n_pv > 0 { q / n_pv as f64 } else { 0.0 };
    let q_per_outlet  = q / n_out as f64;

    let mut total_dp          = 0.0_f64;
    let mut max_main_shear_pa = 0.0_f64;
    let mut total_hi          = 0.0_f64;
    let mut venturi_sigma     = f64::INFINITY;
    let mut venturi_shear_rate = 0.0_f64;
    let mut venturi_shear_pa  = 0.0_f64;
    let mut total_path_len_m  = 0.0_f64;
    let mut residence_time_s  = 0.0_f64;

    // ── Rectangular-duct helper (Shah-London Poiseuille) ─────────────────────
    // Returns (ΔP [Pa], wall_shear_stress [Pa]) for a straight rectangular channel
    // carrying flow `q_ch` with dimensions w_ch × h_ch and length `len`.
    // Wall shear rate at narrow wall: γ ≈ 6 Q / (w h²).
    // ΔP = 12 μ L Q / (w h³) × f_SL   where f_SL = 1 − 0.63/(w/h).
    let rect_metrics = |q_ch: f64, w_ch: f64, h_ch: f64, len: f64| -> (f64, f64) {
        let v       = q_ch / (w_ch * h_ch);
        let gamma   = 6.0 * v / h_ch;
        let mu      = blood.apparent_viscosity(gamma);
        let aspect  = (w_ch / h_ch).max(1.0);
        let f_sl    = 1.0 - 0.63 / aspect;
        let dp      = (12.0 * mu * len * q_ch / (w_ch * h_ch.powi(3)) * f_sl).max(0.0);
        let shear   = (mu * gamma).max(0.0);
        (dp, shear)
    };

    // ── Accumulate one straight-channel segment ───────────────────────────────
    let mut add_rect = |q_ch: f64, w_ch: f64, h_ch: f64, len: f64| {
        let (dp, shear) = rect_metrics(q_ch, w_ch, h_ch, len);
        total_dp          += dp;
        max_main_shear_pa  = max_main_shear_pa.max(shear);
        let t = len * w_ch * h_ch / q_ch.max(1e-12);
        total_hi          += giersiepen_hi(shear, t);
        total_path_len_m  += len;
        residence_time_s  += t;
    };

    // ── Per-topology flow path accumulation ──────────────────────────────────
    match candidate.topology {
        DesignTopology::SingleVenturi | DesignTopology::SerialDoubleVenturi => {
            // Inlet → venturi(s) → outlet.  No distribution tree.
        }
        DesignTopology::BifurcationVenturi => {
            // inlet → trunk(Q) → [Bi] → 2 branches(Q/2) → venturi in each
            let trunk_len  = TREATMENT_HEIGHT_MM * 0.25e-3;
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            add_rect(q,       w, h, trunk_len);   // trunk carries full Q
            add_rect(q / 2.0, w, h, branch_len);  // branch: Q/2 each
        }
        DesignTopology::TrifurcationVenturi => {
            // inlet → trunk(Q) → [Tri] → 3 branches(Q/3) → venturi in each
            let trunk_len  = TREATMENT_HEIGHT_MM * 0.5e-3;
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            add_rect(q,       w, h, trunk_len);
            add_rect(q / 3.0, w, h, branch_len);  // 3 branches, each Q/3
        }
        DesignTopology::VenturiSerpentine | DesignTopology::SerpentineGrid => {
            // Venturi and serpentine handled separately below; no extra trunk.
        }
        DesignTopology::CellSeparationVenturi | DesignTopology::WbcCancerSeparationVenturi => {
            // inlet → trunk(Q) → split → venturi_center(Q/2) + peri_bypass(Q/2)
            // The peripheral bypass carries Q/2 through the same cross-section, so
            // its wall shear is exactly half that of the Q trunk already tracked by
            // add_rect below — the max() is always dominated by the full-Q trunk.
            let trunk_len = TREATMENT_HEIGHT_MM * 0.5e-3;
            add_rect(q, w, h, trunk_len);
        }
        DesignTopology::DoubleBifurcationVenturi => {
            // inlet → trunk(Q) → [Bi] → branch1(Q/2) → [Bi] → 4×(Q/4) → venturi
            let trunk_len   = TREATMENT_HEIGHT_MM * 0.20e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            add_rect(q,       w, h, trunk_len);
            add_rect(q / 2.0, w, h, branch1_len);
            add_rect(q / 4.0, w, h, branch2_len);
        }
        DesignTopology::TripleBifurcationVenturi => {
            // inlet → trunk → [Bi] → [Bi] → [Bi] → 8×(Q/8) → venturi
            let trunk_len   = TREATMENT_HEIGHT_MM * 0.15e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            let branch3_len = TREATMENT_HEIGHT_MM * 0.08e-3;
            add_rect(q,       w, h, trunk_len);
            add_rect(q / 2.0, w, h, branch1_len);
            add_rect(q / 4.0, w, h, branch2_len);
            add_rect(q / 8.0, w, h, branch3_len);
        }
        DesignTopology::DoubleTrifurcationVenturi => {
            // inlet → trunk → [Tri] → [Tri] → 9×(Q/9) → venturi
            let trunk_len   = TREATMENT_HEIGHT_MM * 0.20e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            add_rect(q,       w, h, trunk_len);
            add_rect(q / 3.0, w, h, branch1_len);
            add_rect(q / 9.0, w, h, branch2_len);
        }
        DesignTopology::BifurcationTrifurcationVenturi => {
            // inlet → trunk → [Bi] → [Tri] → 6×(Q/6) → venturi
            let trunk_len   = TREATMENT_HEIGHT_MM * 0.25e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            add_rect(q,       w, h, trunk_len);
            add_rect(q / 2.0, w, h, branch1_len);
            add_rect(q / 6.0, w, h, branch2_len);
        }
        DesignTopology::BifurcationSerpentine => {
            // inlet → trunk(Q) → [Bi] → 2 serpentine arms (Q/2 each)
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            add_rect(q, w, h, trunk_len);
        }
        DesignTopology::TrifurcationSerpentine => {
            // inlet → trunk(Q) → [Tri] → 3 serpentine arms (Q/3 each)
            let trunk_len = TREATMENT_HEIGHT_MM * 0.33e-3;
            add_rect(q, w, h, trunk_len);
        }
        DesignTopology::AdaptiveTree { levels, split_types } => {
            // Variable-depth tree: trunk(Q) + one level_len segment per tree level.
            // Each successive level carries Q / (cumulative fan product).
            if levels > 0 {
                // Distribute total available height evenly across (levels + 1) segments
                let level_len = TREATMENT_HEIGHT_MM * 1e-3 * 0.5 / (levels as f64 + 1.0);
                let mut divisor = 1usize;
                add_rect(q, w, h, level_len); // trunk carries full Q
                for i in 0..levels as usize {
                    let fan = if (split_types >> i) & 1 == 0 { 2 } else { 3 };
                    divisor *= fan;
                    add_rect(q / divisor as f64, w, h, level_len);
                }
            }
            // levels == 0: single straight channel — distribution handled by venturi block
        }
    }
    // Explicitly drop the closure so its mutable borrows are released before
    // the venturi / serpentine sections access those variables directly.
    drop(add_rect);

    // ── Planar rectangular venturi section ───────────────────────────────────
    // Real millifluidic venturis narrow the channel *width* while keeping the
    // channel height constant (photolithographic fabrication).  This gives a
    // velocity ratio of w_ch/w_throat ≈ 3–16, far lower than the circular
    // (D_inlet/D_throat)² = 100 ratio that the VenturiModel assumes.
    //
    // We bypass VenturiModel entirely and use a direct analytical model:
    //   ΔP_total = ΔP_throat_viscous + ΔP_contraction_loss + ΔP_diffuser_loss
    //
    // Where:
    //   ΔP_throat_viscous = Darcy-Weisbach in rectangular throat
    //   ΔP_contraction    = K_c · ½ρV_throat²   (K_c ≈ 0.04 for a smooth entry)
    //   ΔP_diffuser       = K_exp · ½ρ(V_throat − V_inlet)²  (K_exp ≈ 0.14 at 5°)
    if candidate.topology.has_venturi() {
        // Read throat geometry from the blueprint cross-section.
        // Falls back to candidate fields if the blueprint has no throat channel
        // (shouldn't happen for any topology with has_venturi() == true).
        let (w_throat, h_throat) = match throat_cs {
            Some(CrossSectionSpec::Rectangular { width_m, height_m }) => (width_m, height_m),
            _ => (candidate.throat_diameter_m, candidate.channel_height_m),
        };
        let l_throat = candidate.throat_length_m.max(w_throat * 2.0);

        // Cross-sectional areas from blueprint CrossSectionSpec
        let a_throat_rect = w_throat * h_throat;
        let a_inlet_rect  = w * h;

        // Velocities
        let v_thr = q_per_venturi / a_throat_rect.max(1e-18);
        let v_in  = q_per_venturi / a_inlet_rect.max(1e-18);

        // Throat hydraulic diameter D_h = 2wh/(w+h)
        let d_h_thr = 2.0 * w_throat * h_throat / (w_throat + h_throat);

        // Wall shear rate at throat (thin rectangular duct: γ ≈ 6V/h_throat)
        let gamma_thr = 6.0 * v_thr / h_throat;
        let mu_thr    = blood.apparent_viscosity(gamma_thr);

        // Throat Reynolds number + Darcy friction factor
        let re_thr = BLOOD_DENSITY_KG_M3 * v_thr * d_h_thr / mu_thr.max(1e-9);
        let f_thr  = if re_thr < 2300.0 {
            64.0 / re_thr.max(1e-6)
        } else {
            0.3164 / re_thr.powf(0.25)
        };

        // Throat viscous ΔP (Darcy-Weisbach)
        let dp_thr_visc = f_thr * (l_throat / d_h_thr.max(1e-9))
            * 0.5 * BLOOD_DENSITY_KG_M3 * v_thr * v_thr;

        // Contraction minor loss (smooth converging entry, K_c ≈ 0.04)
        let dp_contraction = 0.04 * 0.5 * BLOOD_DENSITY_KG_M3 * v_thr * v_thr;

        // Diffuser expansion loss (5° half-angle, K_exp ≈ 0.14, Idelchik 2007)
        let dv = (v_thr - v_in).max(0.0);
        let dp_diffuser = 0.14 * 0.5 * BLOOD_DENSITY_KG_M3 * dv * dv;

        // Net venturi pressure drop (pump must supply this)
        let dp_venturi = (dp_thr_visc + dp_contraction + dp_diffuser).max(0.0);

        // Cavitation number: σ = (P_inlet_abs − P_vapor) / (½ρV_throat²)
        let p_abs_inlet = (candidate.inlet_gauge_pa - total_dp).max(0.0) + P_ATM_PA;
        let dyn_p_thr   = 0.5 * BLOOD_DENSITY_KG_M3 * v_thr * v_thr;
        let sigma = if v_thr > 1e-3 {
            (p_abs_inlet - BLOOD_VAPOR_PRESSURE_PA) / dyn_p_thr
        } else {
            f64::INFINITY
        };

        let shear_pa = mu_thr * gamma_thr;

        venturi_sigma      = venturi_sigma.min(sigma);
        venturi_shear_rate = venturi_shear_rate.max(gamma_thr);
        venturi_shear_pa   = venturi_shear_pa.max(shear_pa);

        // Accumulate serial ΔP, HI, path length and residence time
        total_dp         += dp_venturi;
        let t_throat      = if v_thr > 1e-9 { l_throat / v_thr } else { 0.0 };
        total_hi         += giersiepen_hi(shear_pa, t_throat);
        // Approximate total venturi path: converging cone + throat + diffuser ≈ 10×d_inlet
        total_path_len_m += candidate.inlet_diameter_m * 10.0 + l_throat;
        residence_time_s += t_throat;

        // ── SerialDoubleVenturi: second venturi stage on the same flow path ──
        // The second stage sees the same Q but a lower available pressure
        // (inlet already reduced by the first venturi drop).  Pre-conditioned
        // fluid lowers cavitation onset: σ₂ < σ₁ is the design intent.
        if candidate.topology == DesignTopology::SerialDoubleVenturi {
            let p2_abs  = (candidate.inlet_gauge_pa - total_dp).max(0.0) + P_ATM_PA;
            let sigma_2 = if v_thr > 1e-3 {
                (p2_abs - BLOOD_VAPOR_PRESSURE_PA) / dyn_p_thr
            } else {
                f64::INFINITY
            };
            venturi_sigma    = venturi_sigma.min(sigma_2);
            total_dp         += dp_venturi;
            total_hi         += giersiepen_hi(shear_pa, t_throat);
            total_path_len_m += candidate.inlet_diameter_m * 10.0 + l_throat;
            residence_time_s += t_throat;
        }
    }

    // ── Serpentine section ────────────────────────────────────────────────────
    // For BifurcationSerpentine / TrifurcationSerpentine each arm carries
    // q_per_outlet (= Q/2 or Q/3).  For VenturiSerpentine / SerpentineGrid
    // there is only one arm so q_per_outlet == Q.
    if candidate.topology.has_serpentine() {
        let q_ser = q_per_outlet;
        let v_ch  = q_ser / (w * h);
        let model = SerpentineModel::<f64>::millifluidic_rectangular(
            w, h,
            candidate.segment_length_m,
            candidate.serpentine_segments,
            candidate.bend_radius_m,
        );
        let mut cond = FlowConditions::<f64>::new(v_ch);
        cond.flow_rate = Some(q_ser);

        let analysis = model.analyze(&blood, &cond).map_err(|e| OptimError::PhysicsError {
            id: candidate.id.clone(),
            reason: format!("SerpentineModel failed: {e}"),
        })?;

        let shear_rate = analysis.wall_shear_rate;
        let shear_pa   = blood.apparent_viscosity(shear_rate) * shear_rate;
        max_main_shear_pa = max_main_shear_pa.max(shear_pa);
        total_dp        += analysis.dp_total;

        let n      = candidate.serpentine_segments as f64;
        let bends  = (n - 1.0).max(0.0);
        let path_m = n * candidate.segment_length_m
            + bends * std::f64::consts::PI * candidate.bend_radius_m;
        let res    = path_m * w * h / q_ser.max(1e-12);
        total_hi        += giersiepen_hi(shear_pa, res);
        total_path_len_m += path_m;
        residence_time_s += res;
    }

    // ── Outlet collection tree (mirror of inlet — symmetric split-merge) ──────
    // Every tree topology is symmetric: the outlet collection tree is the
    // mirror image of the inlet distribution tree.  Blood traverses both halves,
    // so we add the same segment costs again in reverse order (leaf → trunk).
    //
    // IMPORTANT: this block is placed *after* the venturi section so that the σ
    // calculation (which uses `total_dp` at the venturi) only sees inlet-side
    // pressure — adding the outlet tree beforehand would incorrectly reduce the
    // apparent available pressure and understate cavitation potential.
    //
    // Redefine `add_rect` here (inlet version was drop()ed above to release borrows
    // before the venturi/serpentine sections).
    let mut add_rect = |q_ch: f64, w_ch: f64, h_ch: f64, len: f64| {
        let (dp, shear) = rect_metrics(q_ch, w_ch, h_ch, len);
        total_dp          += dp;
        max_main_shear_pa  = max_main_shear_pa.max(shear);
        let t = len * w_ch * h_ch / q_ch.max(1e-12);
        total_hi          += giersiepen_hi(shear, t);
        total_path_len_m  += len;
        residence_time_s  += t;
    };
    match candidate.topology {
        // No symmetric outlet merge tree for these topologies:
        DesignTopology::SingleVenturi
        | DesignTopology::SerialDoubleVenturi
        | DesignTopology::VenturiSerpentine
        | DesignTopology::SerpentineGrid
        | DesignTopology::CellSeparationVenturi
        | DesignTopology::WbcCancerSeparationVenturi => {
            // intentional: either a straight path or a deliberate multi-outlet split
        }

        DesignTopology::BifurcationVenturi => {
            // outlet mirror: 2 branches (Q/2) → merge → trunk (Q)
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            let trunk_len  = TREATMENT_HEIGHT_MM * 0.25e-3;
            add_rect(q / 2.0, w, h, branch_len); // outlet branch
            add_rect(q,       w, h, trunk_len);   // outlet trunk
        }

        DesignTopology::TrifurcationVenturi => {
            // outlet mirror: 3 branches (Q/3) → merge → trunk (Q)
            let branch_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            let trunk_len  = TREATMENT_HEIGHT_MM * 0.50e-3;
            add_rect(q / 3.0, w, h, branch_len); // outlet branch
            add_rect(q,       w, h, trunk_len);   // outlet trunk
        }

        DesignTopology::DoubleBifurcationVenturi => {
            // outlet mirror: leaf(Q/4) → branch1(Q/2) → trunk(Q)
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let trunk_len   = TREATMENT_HEIGHT_MM * 0.20e-3;
            add_rect(q / 4.0, w, h, branch2_len); // outlet leaf
            add_rect(q / 2.0, w, h, branch1_len); // outlet branch
            add_rect(q,       w, h, trunk_len);    // outlet trunk
        }

        DesignTopology::TripleBifurcationVenturi => {
            // outlet mirror: leaf(Q/8) → branch2(Q/4) → branch1(Q/2) → trunk(Q)
            let branch3_len = TREATMENT_HEIGHT_MM * 0.08e-3;
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let trunk_len   = TREATMENT_HEIGHT_MM * 0.15e-3;
            add_rect(q / 8.0, w, h, branch3_len); // outlet leaf
            add_rect(q / 4.0, w, h, branch2_len); // outlet branch 2
            add_rect(q / 2.0, w, h, branch1_len); // outlet branch 1
            add_rect(q,       w, h, trunk_len);    // outlet trunk
        }

        DesignTopology::DoubleTrifurcationVenturi => {
            // outlet mirror: leaf(Q/9) → branch1(Q/3) → trunk(Q)
            let branch2_len = TREATMENT_HEIGHT_MM * 0.10e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.15e-3;
            let trunk_len   = TREATMENT_HEIGHT_MM * 0.20e-3;
            add_rect(q / 9.0, w, h, branch2_len); // outlet leaf
            add_rect(q / 3.0, w, h, branch1_len); // outlet branch
            add_rect(q,       w, h, trunk_len);    // outlet trunk
        }

        DesignTopology::BifurcationTrifurcationVenturi => {
            // outlet mirror: leaf(Q/6) → branch1(Q/2) → trunk(Q)
            let branch2_len = TREATMENT_HEIGHT_MM * 0.12e-3;
            let branch1_len = TREATMENT_HEIGHT_MM * 0.20e-3;
            let trunk_len   = TREATMENT_HEIGHT_MM * 0.25e-3;
            add_rect(q / 6.0, w, h, branch2_len); // outlet leaf
            add_rect(q / 2.0, w, h, branch1_len); // outlet branch
            add_rect(q,       w, h, trunk_len);    // outlet trunk
        }

        DesignTopology::BifurcationSerpentine => {
            // outlet mirror: trunk (Q) collects from 2 serpentine arms
            let trunk_len = TREATMENT_HEIGHT_MM * 0.25e-3;
            add_rect(q, w, h, trunk_len);
        }

        DesignTopology::TrifurcationSerpentine => {
            // outlet mirror: trunk (Q) collects from 3 serpentine arms
            let trunk_len = TREATMENT_HEIGHT_MM * 0.33e-3;
            add_rect(q, w, h, trunk_len);
        }

        DesignTopology::AdaptiveTree { levels, split_types } => {
            // Outlet mirror: iterate levels in reverse (leaf → trunk).
            // Compute total fan product = divisor at deepest level first,
            // then peel back one fan stage per iteration.
            if levels > 0 {
                let level_len = TREATMENT_HEIGHT_MM * 1e-3 * 0.5 / (levels as f64 + 1.0);
                let mut total_fan = 1usize;
                for i in 0..levels as usize {
                    let fan = if (split_types >> i) & 1 == 0 { 2 } else { 3 };
                    total_fan *= fan;
                }
                let mut divisor = total_fan;
                for i in (0..levels as usize).rev() {
                    add_rect(q / divisor as f64, w, h, level_len); // outlet leaf / branch
                    let fan = if (split_types >> i) & 1 == 0 { 2 } else { 3 };
                    divisor /= fan;
                }
                add_rect(q, w, h, level_len); // outlet trunk
            }
        }
    }

    // ── Derived / summary metrics ─────────────────────────────────────────────
    // Flow uniformity: by construction all parallel branches are symmetric → 1.0
    let flow_uniformity = 1.0_f64;

    let cav_potential = if venturi_sigma < SIGMA_CRIT && venturi_sigma >= 0.0 {
        (1.0 - venturi_sigma / SIGMA_CRIT).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let well_coverage_fraction = candidate.topology.nominal_well_coverage();
    let pressure_feasible      = total_dp <= candidate.inlet_gauge_pa;

    // ── 2-population separation (cancer vs RBC, CellSeparationVenturi) ───────
    let sep_metrics = if candidate.topology == DesignTopology::CellSeparationVenturi {
        let cancer = cfd_1d::cell_separation::CellProperties::mcf7_breast_cancer();
        let rbc    = cfd_1d::cell_separation::CellProperties::red_blood_cell();
        let model  = cfd_1d::cell_separation::CellSeparationModel::new(
            w, h, Some(candidate.bend_radius_m),
        );
        let mean_v   = q_per_outlet / (w * h);
        let shear_est = 6.0 * mean_v / h;
        match model.analyze(&cancer, &rbc, blood.density, blood.apparent_viscosity(shear_est), mean_v) {
            Some(a) => (a.separation_efficiency, a.target_center_fraction, a.background_peripheral_fraction),
            None    => (0.0, 0.0, 0.0),
        }
    } else {
        (0.0, 0.0, 0.0)
    };

    // ── 3-population separation (WBC+cancer→center, RBC→periphery) ──────────
    let three_pop_metrics = if candidate.topology == DesignTopology::WbcCancerSeparationVenturi {
        three_population_separation(candidate, &blood)
    } else {
        ThreePopMetrics::default()
    };

    Ok(SdtMetrics {
        cavitation_number: venturi_sigma,
        cavitation_potential: cav_potential,
        throat_shear_rate_inv_s: venturi_shear_rate,
        throat_shear_pa: venturi_shear_pa,
        throat_exceeds_fda: venturi_shear_pa > FDA_MAX_WALL_SHEAR_PA,
        max_main_channel_shear_pa: max_main_shear_pa,
        fda_main_compliant: max_main_shear_pa <= FDA_MAX_WALL_SHEAR_PA,
        hemolysis_index_per_pass: total_hi,
        flow_uniformity,
        well_coverage_fraction,
        mean_residence_time_s: residence_time_s,
        total_pressure_drop_pa: total_dp,
        total_path_length_mm: total_path_len_m * 1000.0,
        pressure_feasible,
        cell_separation_efficiency: sep_metrics.0,
        cancer_center_fraction: sep_metrics.1,
        rbc_peripheral_fraction: sep_metrics.2,
        three_pop_sep_efficiency: three_pop_metrics.sep_efficiency,
        wbc_center_fraction: three_pop_metrics.wbc_center_fraction,
        rbc_peripheral_fraction_three_pop: three_pop_metrics.rbc_peripheral_fraction,
        wbc_equilibrium_pos: three_pop_metrics.wbc_eq_pos,
        cancer_equilibrium_pos: three_pop_metrics.cancer_eq_pos,
        rbc_equilibrium_pos: three_pop_metrics.rbc_eq_pos,
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

// ── Three-population separation ──────────────────────────────────────────────

/// Intermediate results from the three-population inertial separation model.
///
/// Computed by [`three_population_separation`] for the
/// [`DesignTopology::WbcCancerSeparationVenturi`] topology only.
#[derive(Debug, Clone)]
struct ThreePopMetrics {
    /// `x̃_rbc − max(x̃_cancer, x̃_wbc)` clamped to `[−1, 1]`.
    /// Positive = RBCs are more peripheral than both WBCs and cancer cells.
    sep_efficiency: f64,
    /// Fraction of WBCs focused into the center channel (`x̃ < 0.3`).
    wbc_center_fraction: f64,
    /// Fraction of RBCs in the peripheral channels (`x̃ > 0.3`).
    rbc_peripheral_fraction: f64,
    /// WBC normalised equilibrium position `x̃ ∈ [0, 1]`.
    wbc_eq_pos: f64,
    /// Cancer cell normalised equilibrium position `x̃ ∈ [0, 1]`.
    cancer_eq_pos: f64,
    /// RBC normalised equilibrium position `x̃ ∈ [0, 1]`.
    rbc_eq_pos: f64,
}

impl Default for ThreePopMetrics {
    fn default() -> Self {
        Self {
            sep_efficiency: 0.0,
            wbc_center_fraction: 0.0,
            rbc_peripheral_fraction: 0.0,
            wbc_eq_pos: 0.5,
            cancer_eq_pos: 0.5,
            rbc_eq_pos: 0.5,
        }
    }
}

/// Compute three-population inertial focusing equilibrium positions for
/// WBC + cancer → center and RBC → periphery.
///
/// Calls the same [`lateral_equilibrium`] solver used by the 2-population model
/// (Di Carlo 2009 lift + Gossett 2009 Dean drag) independently for each of the
/// three cell types.
///
/// # Channel-sizing rationale
/// WBC focusing requires κ_WBC = D_WBC / D_h > 0.07, i.e. D_h < 143 µm.
/// The effective width is capped at **200 µm** so D_h ≤ 133 µm
/// (κ_WBC ≈ 0.075) regardless of the candidate's nominal channel width.
/// This models the narrower serpentine section used in `WbcCancerSeparationVenturi`.
///
/// Cells with κ < 0.07 are treated as dispersed at x̃ = 0.5.
fn three_population_separation(
    candidate: &DesignCandidate,
    blood: &CassonBlood<f64>,
) -> ThreePopMetrics {
    use cfd_1d::cell_separation::{lateral_equilibrium, CellProperties};

    let cancer = CellProperties::mcf7_breast_cancer();
    let wbc    = CellProperties::white_blood_cell();
    let rbc    = CellProperties::red_blood_cell();

    let w = candidate.channel_width_m;
    let h = candidate.channel_height_m;

    // Cap width so that κ_WBC = 10µm / D_h > 0.07.
    // At w_eff = 200µm, h = 100µm: D_h = 2·200·100/300 ≈ 133µm → κ_WBC ≈ 0.075 ✓
    let w_eff = w.min(200e-6_f64);

    let mean_v = candidate.flow_rate_m3_s / (w_eff * h);
    let shear_est = 6.0 * mean_v / h;
    let mu = blood.apparent_viscosity(shear_est);
    let rho = BLOOD_DENSITY_KG_M3;

    // Avoid extreme curvature: bend radius ≥ 5× D_h.
    let dh = 2.0 * w_eff * h / (w_eff + h);
    let r_bend = candidate.bend_radius_m.max(dh * 5.0);

    const DISPERSED: f64 = 0.5;
    const SPLIT: f64 = 0.3; // center-channel boundary

    let eq = |cell: &CellProperties| -> f64 {
        lateral_equilibrium(cell, rho, mu, mean_v, w_eff, h, Some(r_bend))
            .map(|r| if r.will_focus { r.x_tilde_eq } else { DISPERSED })
            .unwrap_or(DISPERSED)
    };

    let cancer_x = eq(&cancer);
    let wbc_x    = eq(&wbc);
    let rbc_x    = eq(&rbc);

    // Efficiency: how far RBCs are from the more-central of the two treatment populations.
    let max_central_x = cancer_x.max(wbc_x);
    let sep_efficiency = (rbc_x - max_central_x).clamp(-1.0, 1.0);

    // WBC center fraction: fraction focused inside center-channel boundary.
    let wbc_center_fraction = if wbc_x < SPLIT {
        (1.0 - wbc_x / SPLIT).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // RBC peripheral fraction: fraction outside center-channel boundary.
    let rbc_peripheral_fraction = if rbc_x > SPLIT {
        ((rbc_x - SPLIT) / (1.0 - SPLIT)).clamp(0.0, 1.0)
    } else {
        0.0
    };

    ThreePopMetrics {
        sep_efficiency,
        wbc_center_fraction,
        rbc_peripheral_fraction,
        wbc_eq_pos: wbc_x,
        cancer_eq_pos: cancer_x,
        rbc_eq_pos: rbc_x,
    }
}


#[cfg(test)]
mod diag {
    use super::*;
    use crate::design::build_candidate_space;
    use crate::scoring::{score_candidate, OptimMode, SdtWeights};
    #[test]
    fn show_first_metrics() {
        let candidates = build_candidate_space();
        let weights = SdtWeights::default();
        for c in candidates.iter().take(9) {
            if let Ok(m) = compute_metrics(c) {
                let score = score_candidate(&m, OptimMode::SdtCavitation, &weights);
                let id_trunc = if c.id.len() > 30 { &c.id[..30] } else { &c.id };
                println!(
                    "{} | dp={:.0}Pa gauge={:.0}Pa feasible={} fda={} sigma={:.3} score={:.4}",
                    id_trunc, m.total_pressure_drop_pa, c.inlet_gauge_pa,
                    m.pressure_feasible, m.fda_main_compliant, m.cavitation_number, score
                );
            }
        }
    }
}
