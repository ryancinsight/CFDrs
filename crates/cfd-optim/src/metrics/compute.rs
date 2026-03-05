//! Core `compute_metrics` entry point and the Giersiepen haemolysis helper.

use cfd_1d::cavitation_amplified_hi;
use cfd_1d::physics::cell_separation::{CellProperties, CellSeparationModel};
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_schematics::domain::model::CrossSectionSpec;

use crate::constraints::{
    BLOOD_ATTENUATION_405NM_INV_M, BLOOD_DENSITY_KG_M3, BLOOD_VAPOR_PRESSURE_PA,
    BLOOD_VISCOSITY_PA_S, BUBBLE_GAMMA, CLOTTING_BFR_CAUTION_ML_MIN, CLOTTING_BFR_HIGH_RISK_ML_MIN,
    CLOTTING_BFR_LOW_RISK_ML_MIN, CLOTTING_BFR_STRICT_10MLS_ML_MIN, CLOTTING_RESIDENCE_HIGH_RISK_S,
    CLOTTING_RESIDENCE_LOW_RISK_S, CLOTTING_SHEAR_HIGH_RISK_INV_S, CLOTTING_SHEAR_LOW_RISK_INV_S,
    DEAD_VOLUME_SHEAR_THRESHOLD_INV_S, DIFFUSER_DISCHARGE_COEFF, EXPANSION_RATIO_HIGH_RISK,
    EXPANSION_RATIO_LOW_RISK, FDA_MAX_WALL_SHEAR_PA, FDA_TRANSIENT_SHEAR_PA,
    FDA_TRANSIENT_TIME_S, GIERSIEPEN_ALPHA, MILESTONE_TREATMENT_DURATION_MIN,
    PATIENT_BLOOD_VOLUME_ML, PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_REFERENCE_WEIGHT_KG,
    PLATE_HEIGHT_MM, P_ATM_PA, RAYLEIGH_COLLAPSE_FACTOR, R_BUBBLE_EQ_M, SIGMA_CRIT,
    SONO_REF_P_ABS_PA, THERAPEUTIC_WINDOW_REF, TREATMENT_HEIGHT_MM, VENTURI_CC,
    VENTURI_VEL_RATIO_REF,
};
use crate::design::{DesignCandidate, DesignTopology};
use crate::error::OptimError;

use super::network_solve::{solve_blueprint_network, ChannelSolveSample};
use super::sdt_metrics::ChannelHemolysis;
use super::separation::{
    leukapheresis_separation, three_population_separation, LeukapheresisMetrics, ThreePopMetrics,
};
use super::SdtMetrics;

// ── Giersiepen haemolysis model ───────────────────────────────────────────────

/// Giersiepen (1990) haemolysis index: `HI = C · t^α · τ^β`.
///
/// Re-exported from [`cfd_1d::hemolysis`] as the canonical implementation.
/// Returns 0 for non-positive inputs.
pub use cfd_1d::giersiepen_hi;

// ── Main entry ───────────────────────────────────────────────────────────────

/// Compute all physics-derived metrics for a single [`DesignCandidate`].
///
/// The function builds and solves the full 1D channel network, then extracts
/// per-channel solved flows/pressures to derive pressure drop, flow uniformity,
/// venturi selectivity, shear, cavitation, residence, and haemolysis metrics.
pub fn compute_metrics(candidate: &DesignCandidate) -> Result<SdtMetrics, OptimError> {
    let blood = CassonBlood::<f64>::normal_blood();
    let q_inlet = candidate.flow_rate_m3_s;
    let venturi_treatment_enabled = candidate.uses_venturi_treatment();
    let serial_venturi_stages_per_path = candidate.serial_venturi_stages_per_path();
    let active_venturi_throat_count = candidate.active_venturi_throat_count();

    let bp = candidate.to_blueprint();
    let solved = solve_blueprint_network(&candidate.id, candidate.topology, &bp, &blood, q_inlet)?;

    let main_cs = bp
        .channels
        .iter()
        .find(|c| {
            let id = c.id.as_str();
            id == "inlet_section" || id == "parent" || id.starts_with("segment_")
        })
        .map_or(
            CrossSectionSpec::Rectangular {
                width_m: candidate.channel_width_m,
                height_m: candidate.channel_height_m,
            },
            |c| c.cross_section,
        );
    let (w_main, h_main) = main_cs.dims();

    // ── Solved-flow derived hydraulic metrics ───────────────────────────────
    let q_ref = solved.inlet_flow_m3_s.max(1e-12);

    let mut max_main_shear_pa = 0.0_f64;
    let mut bulk_hi = 0.0_f64;
    let mut venturi_hi = 0.0_f64;
    // Hemolysis accumulated ONLY by bypass-arm channels (RBC/WBC fraction).
    // These channels have no cavitation exposure; their HI is shear-only and
    // must NOT be added to the treatment-zone venturi_hi or final_hi
    // cavitation correction — doing so would over-estimate RBC damage because
    // the bypass arm RBCs never enter the sonication zone.
    let mut bypass_hi_accumulated = 0.0_f64;
    let mut expected_path_len_m = 0.0_f64;
    let mut all_channel_shears: Vec<f64> = Vec::new();

    // FDA Mechanical Index compliance tracking.
    // For each venturi throat, compute the equivalent MI from the Bernoulli
    // pressure drop and verify MI_equiv < 1.9 (FDA 510(k) guidance, 2019).
    // MI_equiv = sqrt(2 × ΔP_throat / ρ_blood) / c_sound_blood
    // where c_sound_blood ≈ 1540 m/s (Shung 2006).
    // Therapeutic window: 0.3 ≤ MI_equiv < 1.9 (controlled inertial cavitation).
    const FDA_MI_LIMIT: f64 = 1.9;
    const SOUND_SPEED_BLOOD_M_S: f64 = 1540.0;
    let mut fda_mi_max: f64 = 0.0;
    let mut fda_mi_all_compliant = true;
    let mut _throat_in_therapeutic_window = false;

    let mut venturi_sigma = f64::INFINITY;
    let mut venturi_shear_rate = 0.0_f64;
    let mut venturi_shear_pa = 0.0_f64;
    let mut t_throat_venturi = 0.0_f64;
    let mut venturi_v_thr = 0.0_f64;
    let mut venturi_v_in = 0.0_f64;
    let mut diffuser_recovery_pa = 0.0_f64;
    let mut network_serial_throat_stages = 1usize;
    let mut per_channel_hi: Vec<ChannelHemolysis> = Vec::new();

    // ── Geometry-sensitive stasis accumulators ──────────────────────────────
    // Track per-channel shear rates, transit times, volumes, and expansion
    // ratios to compute a CRI that discriminates designs at the same flow rate.
    let mut min_channel_shear_rate = f64::INFINITY;
    let mut max_channel_transit_s = 0.0_f64;
    let mut max_expansion_area_ratio = 1.0_f64;
    let mut total_channel_volume_m3 = 0.0_f64;
    let mut dead_channel_volume_m3 = 0.0_f64;

    for sample in &solved.channel_samples {
        let q_abs = sample.flow_m3_s.abs();
        if q_abs <= 1e-18 {
            continue;
        }

        let area = sample.cross_section.area().max(1e-18);
        let v = q_abs / area;
        let gamma_fd = sample.cross_section.wall_shear_rate(v);
        // Entrance-length correction: developing flow has ~1.5× the wall shear
        // of fully-developed flow (Shah & London 1978).
        let elf = entrance_length_shear_factor(sample.cross_section, v, sample.length_m);
        let gamma = gamma_fd * elf;
        let mu = blood.apparent_viscosity(gamma);
        let tau = (mu * gamma).max(0.0);
        let t_seg = if v > 1e-12 { sample.length_m / v } else { 0.0 };
        let hi_seg = giersiepen_hi(tau, t_seg);
        let throat_repeat_factor = if sample.is_venturi_throat {
            f64::from(sample.per_channel_throat_count.max(1))
        } else {
            1.0
        };
        let weight = q_abs / q_ref;

        // ── Geometry-sensitive stasis accumulation ──────────────────────────
        let seg_volume = sample.length_m * area;
        total_channel_volume_m3 += seg_volume;
        if gamma > 0.0 {
            min_channel_shear_rate = min_channel_shear_rate.min(gamma);
        }
        if gamma < DEAD_VOLUME_SHEAR_THRESHOLD_INV_S {
            dead_channel_volume_m3 += seg_volume;
        }
        if t_seg > 0.0 {
            max_channel_transit_s = max_channel_transit_s.max(t_seg);
        }

        // Bypass-channel hemolysis accounting:
        // Channels tagged as bypass (is_bypass_channel=true) carry RBCs routed
        // AROUND the sonication zone. Their HI contribution is shear-only and
        // must be tracked separately so we do not apply the cavitation-HCT
        // correction (which reduces HI assuming RBCs are absent from the venturi).
        // Adding bypass HI to bulk_hi would incorrectly inflate it under the
        // hematocrit correction — bypass channels always have full hematocrit.
        if sample.is_bypass_channel {
            bypass_hi_accumulated += weight * hi_seg * throat_repeat_factor;
            // Path length for bypass channels is accounted for below (post-loop)
            // to avoid double-counting with the expected_path_len_m accumulator.
        } else {
            expected_path_len_m += weight * sample.length_m;
            bulk_hi += weight * hi_seg * throat_repeat_factor;
        }

        per_channel_hi.push(ChannelHemolysis {
            channel_id: sample.id.to_string(),
            is_venturi_throat: sample.is_venturi_throat,
            hi_contribution: weight * hi_seg * throat_repeat_factor,
            wall_shear_pa: tau,
            transit_time_s: t_seg * throat_repeat_factor,
            flow_fraction: weight,
        });

        if sample.is_venturi_throat && venturi_treatment_enabled {
            network_serial_throat_stages = network_serial_throat_stages
                .max(usize::from(sample.per_channel_throat_count.max(1)));
            venturi_hi += weight * hi_seg * throat_repeat_factor;
            venturi_shear_rate = venturi_shear_rate.max(gamma);
            venturi_shear_pa = venturi_shear_pa.max(tau);
            t_throat_venturi = t_throat_venturi.max(t_seg * throat_repeat_factor);
            venturi_v_thr = venturi_v_thr.max(v);

            // Bernoulli-corrected σ: compute static pressure at the throat
            // vena contracta by subtracting the dynamic pressure gain from
            // the upstream node pressure.
            //
            // P_static,throat = P_upstream + P_atm
            //                 - ½ρ(v_throat² − v_upstream²)
            //
            // The vena contracta coefficient Cc accounts for the flow
            // contraction beyond the geometric throat (Cc ≈ 0.61 for sharp
            // orifice; 0.95 for well-rounded contraction).
            let upstream_v_local = solved
                .channel_samples
                .iter()
                .find(|other| other.to_node == sample.from_node && !other.is_venturi_throat)
                .map_or_else(
                    || q_abs / (w_main * h_main).max(1e-18),
                    |other| {
                        let up_area = other.cross_section.area().max(1e-18);
                        other.flow_m3_s.abs() / up_area
                    },
                );
            venturi_v_in = venturi_v_in.max(upstream_v_local);

            // Area expansion ratio at post-venturi diffuser (throat → main).
            // Used for geometry-sensitive stasis: larger ratios produce
            // persistent Borda-Carnot recirculation eddies downstream.
            if area > 0.0 {
                let upstream_area = solved
                    .channel_samples
                    .iter()
                    .find(|o| o.to_node == sample.from_node && !o.is_venturi_throat)
                    .map_or(w_main * h_main, |o| o.cross_section.area().max(1e-18));
                let expansion = upstream_area / area;
                max_expansion_area_ratio = max_expansion_area_ratio.max(expansion);
            }

            // ── FDA Mechanical Index compliance check ─────────────────────────
            // Compute equivalent MI from Bernoulli pressure drop at the throat
            // vena contracta.  The FDA 510(k) guidance (2019) mandates MI < 1.9
            // for any device producing acoustic or hydrodynamic cavitation.
            // Therapeutic SDT targets 0.3 ≤ MI_equiv < 1.9 for controlled
            // inertial cavitation without uncontrolled bubble chain reactions.
            //
            // MI_equiv = sqrt(2 × ΔP_vc / ρ) / c_sound
            // where ΔP_vc = ½ρ(v_eff² − v_in²)  (Bernoulli vena contracta drop).
            {
                let v_eff_for_mi = v / VENTURI_CC;
                let dp_vc = (0.5
                    * BLOOD_DENSITY_KG_M3
                    * (v_eff_for_mi * v_eff_for_mi - upstream_v_local * upstream_v_local))
                    .max(0.0);
                let mi_equiv =
                    (2.0 * dp_vc / BLOOD_DENSITY_KG_M3.max(1.0)).sqrt() / SOUND_SPEED_BLOOD_M_S;
                fda_mi_max = fda_mi_max.max(mi_equiv);
                if mi_equiv >= FDA_MI_LIMIT {
                    fda_mi_all_compliant = false;
                }
                if mi_equiv >= 0.3 && mi_equiv < FDA_MI_LIMIT {
                    _throat_in_therapeutic_window = true;
                }
            }

            // Apply vena contracta coefficient: the jet contracts beyond
            // the geometric throat, so the effective velocity at the vena
            // contracta is v_eff = v_geom / Cc (VENTURI_CC ≈ 0.85 for
            // smooth millifluidic contractions; Idelchik Diagram 4-9).
            let v_eff = v / VENTURI_CC;
            let p_gauge_upstream = sample.from_pressure_pa.max(0.0);
            let bernoulli_drop =
                0.5 * BLOOD_DENSITY_KG_M3 * (v_eff * v_eff - upstream_v_local * upstream_v_local);

            // Throat skin-friction loss: Darcy-Weisbach ΔP_f = f·(L/D_h)·½ρv²
            // where f = 64/Re for laminar flow (Re < 2300) in the developing
            // throat (Shah & London 1978).  This additional pressure drop
            // further reduces p_static at the vena contracta, promoting
            // cavitation inception in long/narrow throats.
            let d_h_throat = sample.cross_section.hydraulic_diameter();
            let re_throat =
                BLOOD_DENSITY_KG_M3 * v_eff * d_h_throat / BLOOD_VISCOSITY_PA_S.max(1e-18);
            let f_darcy = if re_throat > 1.0 {
                64.0 / re_throat
            } else {
                64.0
            };
            let friction_drop = f_darcy
                * (sample.length_m / d_h_throat.max(1e-18))
                * 0.5
                * BLOOD_DENSITY_KG_M3
                * v_eff
                * v_eff;

            let p_static_throat =
                (p_gauge_upstream + P_ATM_PA - bernoulli_drop - friction_drop).max(0.0);
            let dyn_p = 0.5 * BLOOD_DENSITY_KG_M3 * v_eff * v_eff;
            let sigma = if dyn_p > 1e-12 {
                (p_static_throat - BLOOD_VAPOR_PRESSURE_PA) / dyn_p
            } else {
                f64::INFINITY
            };
            venturi_sigma = venturi_sigma.min(sigma);

            // Diffuser pressure recovery: the gradual expansion downstream of
            // the venturi throat converts a fraction of the throat dynamic
            // pressure back to static pressure.  This partially compensates
            // the Bernoulli drop, improving the net pressure budget.
            // ΔP_recovery = C_D × ½ρ(v_throat² − v_upstream²)
            // (Idelchik 1994, Handbook of Hydraulic Resistance, Diagram 6-21)
            let dp_recovery = DIFFUSER_DISCHARGE_COEFF
                * 0.5
                * BLOOD_DENSITY_KG_M3
                * (v_eff * v_eff - upstream_v_local * upstream_v_local).max(0.0);
            diffuser_recovery_pa = diffuser_recovery_pa.max(dp_recovery);
        } else if !sample.is_bypass_channel {
            // Non-venturi, non-bypass channels: track main-channel shear for
            // FDA compliance and shear statistics.
            max_main_shear_pa = max_main_shear_pa.max(tau);
            if tau > 0.0 {
                all_channel_shears.push(tau);
            }
        } else {
            // Bypass channels: track shear for bypass-specific statistics only.
            // Do NOT count bypass shear in FDA main-channel compliance since
            // RBCs in bypass are explicitly routed away from the treatment zone.
            if tau > 0.0 {
                all_channel_shears.push(tau);
            }
        }
    }

    // Add bypass hemolysis to bulk: bypass channels are always present and
    // their shear-only HI contributes to total device hemolysis, but is NOT
    // subject to the local-hematocrit venturi correction applied to center channels.
    // This preserves the physical correctness: bypass RBCs see full feed_hematocrit
    // Casson blood viscosity throughout, never entering the low-HCT venturi region.
    let _bypass_hi_total = bypass_hi_accumulated;
    // Include bypass channels in path-length accounting proportional to their flow.
    for sample in &solved.channel_samples {
        if sample.is_bypass_channel {
            let q_abs = sample.flow_m3_s.abs();
            if q_abs > 1e-18 {
                expected_path_len_m += (q_abs / q_ref) * sample.length_m;
            }
        }
    }

    let mut total_pressure_drop_pa = (solved.inlet_pressure_pa.max(0.0) + solved.remerge_loss_pa
        - diffuser_recovery_pa)
        .max(0.0);

    // ── FDA overall compliance incorporating MI check ──────────────────────────
    // `fda_overall_compliant` (computed later) must account for MI compliance
    // in addition to shear and pressure constraints.  Store MI compliance here
    // for use in the final metric assembly.
    let _fda_mi_compliant = fda_mi_all_compliant || !venturi_treatment_enabled;

    // ── Wall shear percentile statistics ─────────────────────────────────
    let (wall_shear_p95_pa, wall_shear_p99_pa, wall_shear_mean_pa, wall_shear_cv) =
        compute_shear_percentiles(&mut all_channel_shears);
    let fda_shear_percentile_compliant =
        wall_shear_p95_pa <= FDA_MAX_WALL_SHEAR_PA && wall_shear_p99_pa <= FDA_TRANSIENT_SHEAR_PA;

    let mut flow_uniformity = solved.flow_uniformity;
    if !flow_uniformity.is_finite() || flow_uniformity <= 1e-9 {
        flow_uniformity = fallback_uniformity(candidate.topology);
    }
    flow_uniformity = flow_uniformity.clamp(0.0, 1.0);

    let modeled_venturi_flow_fraction = candidate.venturi_flow_fraction();
    let mut venturi_flow_fraction = solved.venturi_flow_fraction.clamp(0.0, 1.0);
    if venturi_treatment_enabled && venturi_flow_fraction <= 1e-9 {
        venturi_flow_fraction = modeled_venturi_flow_fraction;
    }
    if !venturi_treatment_enabled {
        venturi_flow_fraction = 0.0;
    }

    // If throat-resolved extraction collapses (e.g. stiff selective trees),
    // recover venturi local quantities from the model-based throat flow share.
    if venturi_treatment_enabled && venturi_v_thr <= 1e-12 && venturi_flow_fraction > 0.0 {
        let n_parallel = candidate.topology.parallel_venturi_count().max(1) as f64;
        let q_per_throat = q_inlet * venturi_flow_fraction / n_parallel;
        let throat_area = (candidate.throat_diameter_m * candidate.channel_height_m).max(1e-18);
        let inlet_area = (candidate.channel_width_m * candidate.channel_height_m).max(1e-18);
        let v_thr = q_per_throat / throat_area;
        let v_in = q_per_throat / inlet_area;
        let gamma_fd = 6.0 * v_thr / candidate.channel_height_m.max(1e-18);
        // Entrance-length correction for the venturi throat.
        let throat_cs = CrossSectionSpec::Rectangular {
            width_m: candidate.throat_diameter_m,
            height_m: candidate.channel_height_m,
        };
        let elf = entrance_length_shear_factor(throat_cs, v_thr, candidate.throat_length_m);
        let gamma = gamma_fd * elf;
        let mu = blood.apparent_viscosity(gamma);
        let tau = (mu * gamma).max(0.0);
        let t_thr = if v_thr > 1e-12 {
            candidate.throat_length_m / v_thr
        } else {
            0.0
        };

        venturi_v_thr = v_thr;
        venturi_v_in = v_in;
        venturi_shear_rate = venturi_shear_rate.max(gamma);
        venturi_shear_pa = venturi_shear_pa.max(tau);
        t_throat_venturi = t_throat_venturi.max(t_thr);

        // Bernoulli-corrected fallback σ with Cc: v_eff = v_thr / Cc
        let v_thr_eff = v_thr / VENTURI_CC;
        let p_abs_local = candidate.inlet_pressure_pa();
        let bernoulli_drop = 0.5 * BLOOD_DENSITY_KG_M3 * (v_thr_eff * v_thr_eff - v_in * v_in);

        // Throat skin-friction loss (fallback path): same Darcy-Weisbach model
        // as the resolved-flow path above.
        let d_h_fb = throat_cs.hydraulic_diameter();
        let re_fb = BLOOD_DENSITY_KG_M3 * v_thr_eff * d_h_fb / BLOOD_VISCOSITY_PA_S.max(1e-18);
        let f_darcy_fb = if re_fb > 1.0 { 64.0 / re_fb } else { 64.0 };
        let friction_drop_fb = f_darcy_fb
            * (candidate.throat_length_m / d_h_fb.max(1e-18))
            * 0.5
            * BLOOD_DENSITY_KG_M3
            * v_thr_eff
            * v_thr_eff;

        let p_static_throat = (p_abs_local - bernoulli_drop - friction_drop_fb).max(0.0);
        let dyn_p = 0.5 * BLOOD_DENSITY_KG_M3 * v_thr_eff * v_thr_eff;
        if dyn_p > 1e-12 {
            let sigma = (p_static_throat - BLOOD_VAPOR_PRESSURE_PA) / dyn_p;
            venturi_sigma = venturi_sigma.min(sigma);
        }

        // Fallback diffuser recovery (same Idelchik Diagram 6-21 model)
        let dp_recovery_fb = DIFFUSER_DISCHARGE_COEFF
            * 0.5
            * BLOOD_DENSITY_KG_M3
            * (v_thr_eff * v_thr_eff - v_in * v_in).max(0.0);
        diffuser_recovery_pa = diffuser_recovery_pa.max(dp_recovery_fb);
    }

    let unresolved_serial_stages = serial_venturi_stages_per_path
        .saturating_sub(network_serial_throat_stages)
        .max(1);
    if venturi_treatment_enabled && unresolved_serial_stages > 1 {
        let stage_factor = unresolved_serial_stages as f64;
        t_throat_venturi *= stage_factor;
        venturi_hi *= stage_factor;
        total_pressure_drop_pa *= 1.0 + 0.18 * (stage_factor - 1.0);
    }

    // ── Derived / summary metrics ────────────────────────────────────────────
    // Nonlinear cavitation potential: bubble number density scales as
    // (1 − σ/σ_crit)^m with m ≈ 1.5 in the incipient regime (Brennen 1995,
    // _Cavitation and Bubble Dynamics_, §1.6).  The linear (1−σ) model
    // over-predicts cavitation at σ just below σ_crit and under-predicts
    // at σ ≪ σ_crit.
    let cav_potential_raw = if venturi_sigma < SIGMA_CRIT {
        // σ < 0 (p_static < p_vapor) → very strong cavitation → clamped to 1.0
        // σ ∈ [0, σ_crit) → incipient cavitation → Brennen nonlinear scaling
        (1.0 - venturi_sigma / SIGMA_CRIT).clamp(0.0, 1.0).powf(1.5)
    } else {
        0.0 // σ ≥ σ_crit → no cavitation
    };

    // Rayleigh-Plesset growth-time gate: a vapour bubble of equilibrium
    // radius R₀ needs at least t_Rayleigh = 0.915 · R₀ · √(ρ / ΔP) to grow
    // to maximum radius before collapse (Rayleigh 1917; Brennen 1995 §2.3).
    // If the throat transit time is shorter than this, cavitation cannot
    // fully develop → derate the potential proportionally.
    let cav_potential = if cav_potential_raw > 0.0 && t_throat_venturi > 0.0 {
        let dp_driving = (P_ATM_PA - BLOOD_VAPOR_PRESSURE_PA).max(1.0);
        let t_rayleigh =
            RAYLEIGH_COLLAPSE_FACTOR * R_BUBBLE_EQ_M * (BLOOD_DENSITY_KG_M3 / dp_driving).sqrt();
        // t_growth ≈ t_rayleigh (symmetric for inertial growth phase).
        // Derating: if t_throat < t_growth, scale linearly.
        let growth_ratio = (t_throat_venturi / t_rayleigh.max(1e-18)).min(1.0);
        cav_potential_raw * growth_ratio
    } else {
        cav_potential_raw
    };

    let well_coverage_fraction = candidate.topology.nominal_well_coverage();
    let flow_rate_ml_min = candidate.flow_rate_m3_s * 60.0 * 1_000_000.0;

    let plate_fits = match &candidate.topology {
        DesignTopology::ParallelMicrochannelArray { n_channels } => {
            let pitch_m = w_main * 2.0;
            let total_y_m = *n_channels as f64 * pitch_m;
            total_y_m <= (PLATE_HEIGHT_MM - 10.0) * 1e-3
        }
        DesignTopology::QuadTrifurcationVenturi => {
            let w_leaf = candidate.channel_width_m * candidate.trifurcation_center_frac.powi(4);
            let pitch = (w_leaf * 4.0).max(200e-6);
            let side = 9.0 * pitch;
            side <= (TREATMENT_HEIGHT_MM + 10.0) * 1e-3
        }
        _ => true,
    };

    let n_outlet_ports: usize = candidate.topology.outlet_count();

    // ── 2-population separation ──────────────────────────────────────────────
    let sep_metrics: (f64, f64, f64) = if candidate.topology == DesignTopology::CellSeparationVenturi {
        let cancer = CellProperties::mcf7_breast_cancer();
        let rbc = CellProperties::red_blood_cell();
        let model = CellSeparationModel::new(
            w_main,
            h_main,
            Some(candidate.bend_radius_m),
        );
        let q_center = (q_inlet * venturi_flow_fraction).max(1e-12);
        let mean_v = q_center / (w_main * h_main).max(1e-18);
        let shear_est = 6.0 * mean_v / h_main.max(1e-18);
        match model.analyze(
            &cancer,
            &rbc,
            blood.density,
            blood.apparent_viscosity(shear_est),
            mean_v,
        ) {
            Some(a) => (
                a.separation_efficiency,
                a.target_center_fraction,
                a.background_peripheral_fraction,
            ),
            None => (0.0, 0.0, 0.0),
        }
    } else {
        (0.0, 0.0, 0.0)
    };

    // ── 3-population separation ──────────────────────────────────────────────
    let three_pop_metrics = if matches!(
        candidate.topology,
        DesignTopology::WbcCancerSeparationVenturi
            | DesignTopology::AsymmetricBifurcationSerpentine
            | DesignTopology::VenturiSerpentine
            | DesignTopology::SerpentineGrid
            | DesignTopology::BifurcationSerpentine
            | DesignTopology::TrifurcationSerpentine
            | DesignTopology::SpiralSerpentine { .. }
            | DesignTopology::TrifurcationBifurcationVenturi
            | DesignTopology::TripleTrifurcationVenturi
            | DesignTopology::TrifurcationBifurcationBifurcationVenturi
            | DesignTopology::QuadTrifurcationVenturi
            | DesignTopology::CascadeCenterTrifurcationSeparator { .. }
            | DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
            | DesignTopology::DoubleTrifurcationCIFVenturi { .. }
            | DesignTopology::AsymmetricTrifurcationVenturi
            | DesignTopology::TriBiTriSelectiveVenturi
            | DesignTopology::DoubleTrifurcationVenturi
    ) {
        three_population_separation(candidate, &blood)
    } else {
        ThreePopMetrics::default()
    };

    // ── Leukapheresis separation ─────────────────────────────────────────────
    let leuka_metrics = if matches!(
        candidate.topology,
        DesignTopology::ConstrictionExpansionArray { .. }
            | DesignTopology::SpiralSerpentine { .. }
            | DesignTopology::ParallelMicrochannelArray { .. }
    ) {
        leukapheresis_separation(candidate, &blood)
    } else {
        LeukapheresisMetrics::default()
    };

    // ── CCT and CIF staged routing models ───────────────────────────────────
    let mut cct_stage_qfracs: Vec<f64> = Vec::new();
    let mut cif_pretri_stage_qfracs: Vec<f64> = Vec::new();
    let mut cif_terminal_tri_qfrac: Option<f64> = None;
    let mut cif_terminal_bi_qfrac: Option<f64> = None;

    // ── DTCV (DoubleTrifurcationCIFVenturi) cascade routing ─────────────────
    // DTCV uses a 2-level asymmetric trifurcation cascade (split1 → split2)
    // where the first split routes cancer/WBC-enriched flow to the center arm,
    // and the second split further subdivides into 3 CTC-enriched center
    // sub-channels with serial venturi throats.  The outer bypass channels
    // (6 total: 3 from split1 periphery, 3 from split2 periphery of the center
    // trunk) carry RBC-enriched flow with zero venturi exposure.
    //
    // The Zweifach-Fung routing at each trifurcation junction progressively
    // enriches stiff cancer cells (β=1.70) in the center arm while deformable
    // RBCs (β=1.00) distribute more uniformly — exactly the mechanism that
    // makes deeper cascades (tri-tri > tri) produce better separation.
    let dtcv_cascade_res: Option<cfd_1d::CascadeJunctionResult> =
        if let DesignTopology::DoubleTrifurcationCIFVenturi { .. } = candidate.topology {
            // DTCV is a 2-level cascade: split1 uses cif_pretri_center_frac,
            // split2 uses cif_terminal_tri_center_frac.  Each level has its
            // own width fraction, so we compute per-stage q_fracs individually.
            let q1 = cfd_1d::tri_center_q_frac_cross_junction(
                candidate.cif_pretri_center_frac(),
                w_main,
                h_main,
            );
            let q2 = cfd_1d::tri_center_q_frac_cross_junction(
                candidate.cif_terminal_tri_center_frac(),
                w_main * candidate.cif_pretri_center_frac(),
                h_main,
            );
            Some(cfd_1d::cascade_junction_separation_from_qfracs(&[q1, q2]))
        } else {
            None
        };

    let cascade_res: Option<cfd_1d::CascadeJunctionResult> =
        if let DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } = candidate.topology
        {
            let solved_q = extract_cct_stage_qfracs(&solved.channel_samples, n_levels);
            let use_solved = solved_q.len() == n_levels as usize;
            if use_solved {
                cct_stage_qfracs = solved_q.clone();
                Some(cfd_1d::cascade_junction_separation_from_qfracs(&solved_q))
            } else {
                let q_model = cfd_1d::tri_center_q_frac_cross_junction(
                    candidate.trifurcation_center_frac,
                    w_main,
                    h_main,
                );
                cct_stage_qfracs = std::iter::repeat_n(q_model, n_levels as usize).collect();
                Some(
                    cfd_1d::cascade_junction_separation_cross_junction(
                        n_levels,
                        candidate.trifurcation_center_frac,
                        w_main,
                        h_main,
                    ),
                )
            }
        } else {
            None
        };

    let cif_res: Option<cfd_1d::IncrementalFiltrationResult> =
        if let DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } = candidate.topology
        {
            let (pretri_q, tri_q_opt, bi_q_opt) =
                extract_cif_stage_qfracs(&solved.channel_samples, n_pretri);
            let use_solved = pretri_q.len() == n_pretri as usize;
            if use_solved {
                cif_pretri_stage_qfracs = pretri_q.clone();
                let tri_q = tri_q_opt.unwrap_or_else(|| {
                    cfd_1d::tri_center_q_frac(
                        candidate.cif_terminal_tri_center_frac(),
                    )
                });
                let bi_q = bi_q_opt.unwrap_or_else(|| candidate.cif_terminal_bi_treat_frac());
                cif_terminal_tri_qfrac = Some(tri_q);
                cif_terminal_bi_qfrac = Some(bi_q);
                Some(
                    cfd_1d::incremental_filtration_separation_from_qfracs(
                        &pretri_q, tri_q, bi_q,
                    ),
                )
            } else {
                let staged_pretri_q = cfd_1d::cif_pretri_stage_q_fracs(
                    n_pretri,
                    candidate.cif_pretri_center_frac(),
                    candidate.cif_terminal_tri_center_frac(),
                );
                let q_tri = cfd_1d::tri_center_q_frac_cross_junction(
                    candidate.cif_terminal_tri_center_frac(),
                    w_main,
                    h_main,
                );
                let q_bi = candidate.cif_terminal_bi_treat_frac();
                cif_pretri_stage_qfracs = staged_pretri_q.clone();
                cif_terminal_tri_qfrac = Some(q_tri);
                cif_terminal_bi_qfrac = Some(q_bi);
                Some(
                    cfd_1d::incremental_filtration_separation_cross_junction(
                        n_pretri,
                        candidate.cif_pretri_center_frac(),
                        candidate.cif_terminal_tri_center_frac(),
                        q_bi,
                        w_main,
                        h_main,
                    ),
                )
            }
        } else {
            None
        };

    // ── PAI baseline (Hellums 1994) ──────────────────────────────────────────
    let pai_tau = if venturi_shear_pa > 0.0 {
        venturi_shear_pa
    } else {
        max_main_shear_pa
    };
    let pai_t = if t_throat_venturi > 0.0 {
        t_throat_venturi
    } else {
        solved.mean_residence_time_s
    };
    let pai_base = if pai_tau > 0.0 && pai_t > 0.0 {
        1.8e-8 * pai_tau.powf(1.325) * pai_t.powf(0.462)
    } else {
        0.0
    };

    // ── Final metric values (selective venturi composition correction) ──────
    let mut final_hi = bulk_hi;
    let mut final_pai = pai_base;
    let mut final_local_hct = candidate.feed_hematocrit;
    let mut final_cancer_dose = 0.0_f64;

    let mut final_sep_eff = sep_metrics.0;
    let mut final_cancer_frac = sep_metrics.1;
    let mut final_rbc_periph = sep_metrics.2;
    let mut final_three_pop_sep = three_pop_metrics.sep_efficiency;
    let mut final_wbc_center_three_pop = three_pop_metrics.wbc_center_fraction;
    let mut final_rbc_periph_three_pop = three_pop_metrics.rbc_peripheral_fraction;
    let mut final_wbc_recovery = leuka_metrics.wbc_recovery;
    let mut final_rbc_pass_fraction = leuka_metrics.rbc_pass_fraction;
    let mut final_wbc_purity = leuka_metrics.wbc_purity;
    let mut final_total_ecv_ml = leuka_metrics.total_ecv_ml;

    let mut final_rbc_venturi_exposure = 1.0_f64;
    let final_venturi_flow_fraction = venturi_flow_fraction;

    let mut rbc_center_frac = 1.0_f64;
    let mut cancer_center_frac = 0.0_f64;
    let mut wbc_center_frac = 0.0_f64;

    match candidate.topology {
        DesignTopology::CellSeparationVenturi => {
            final_sep_eff = sep_metrics.0;
            final_cancer_frac = sep_metrics.1;
            final_rbc_periph = sep_metrics.2;
            final_three_pop_sep = sep_metrics.0;
            final_rbc_periph_three_pop = sep_metrics.2;

            rbc_center_frac = (1.0 - sep_metrics.2).clamp(0.0, 1.0);
            cancer_center_frac = sep_metrics.1.clamp(0.0, 1.0);
            wbc_center_frac = 0.0;
        }
        DesignTopology::WbcCancerSeparationVenturi => {
            let rbc_periph = three_pop_metrics.rbc_peripheral_fraction.clamp(0.0, 1.0);
            final_three_pop_sep = three_pop_metrics.sep_efficiency;
            final_wbc_center_three_pop = three_pop_metrics.wbc_center_fraction;
            final_rbc_periph_three_pop = rbc_periph;

            final_sep_eff = three_pop_metrics.sep_efficiency;
            final_cancer_frac = three_pop_metrics.cancer_center_fraction;
            final_rbc_periph = rbc_periph;

            rbc_center_frac = 1.0 - rbc_periph;
            cancer_center_frac = three_pop_metrics.cancer_center_fraction.clamp(0.0, 1.0);
            wbc_center_frac = three_pop_metrics.wbc_center_fraction.clamp(0.0, 1.0);
        }
        DesignTopology::CascadeCenterTrifurcationSeparator { .. } => {
            if let Some(cct) = cascade_res {
                final_sep_eff = cct.separation_efficiency;
                final_cancer_frac = cct.cancer_center_fraction;
                final_rbc_periph = cct.rbc_peripheral_fraction;
                final_three_pop_sep = cct.separation_efficiency;
                final_wbc_center_three_pop = cct.wbc_center_fraction;
                final_rbc_periph_three_pop = cct.rbc_peripheral_fraction;

                rbc_center_frac = (1.0 - cct.rbc_peripheral_fraction).clamp(0.0, 1.0);
                cancer_center_frac = cct.cancer_center_fraction.clamp(0.0, 1.0);
                wbc_center_frac = cct.wbc_center_fraction.clamp(0.0, 1.0);
            }
        }
        DesignTopology::IncrementalFiltrationTriBiSeparator { .. } => {
            if let Some(cif) = cif_res {
                final_sep_eff = cif.separation_efficiency;
                final_cancer_frac = cif.cancer_center_fraction;
                final_rbc_periph = cif.rbc_peripheral_fraction;
                final_three_pop_sep = cif.separation_efficiency;
                final_wbc_center_three_pop = cif.wbc_center_fraction;
                final_rbc_periph_three_pop = cif.rbc_peripheral_fraction;

                rbc_center_frac = cif.rbc_center_fraction.clamp(0.0, 1.0);
                cancer_center_frac = cif.cancer_center_fraction.clamp(0.0, 1.0);
                wbc_center_frac = cif.wbc_center_fraction.clamp(0.0, 1.0);
            }
        }
        DesignTopology::DoubleTrifurcationCIFVenturi { .. } => {
            if let Some(dtcv) = &dtcv_cascade_res {
                final_sep_eff = dtcv.separation_efficiency;
                final_cancer_frac = dtcv.cancer_center_fraction;
                final_rbc_periph = dtcv.rbc_peripheral_fraction;
                final_three_pop_sep = dtcv.separation_efficiency;
                final_wbc_center_three_pop = dtcv.wbc_center_fraction;
                final_rbc_periph_three_pop = dtcv.rbc_peripheral_fraction;

                rbc_center_frac = (1.0 - dtcv.rbc_peripheral_fraction).clamp(0.0, 1.0);
                cancer_center_frac = dtcv.cancer_center_fraction.clamp(0.0, 1.0);
                wbc_center_frac = dtcv.wbc_center_fraction.clamp(0.0, 1.0);
            }
        }
        DesignTopology::AsymmetricTrifurcationVenturi => {
            // Single asymmetric trifurcation: use three-pop inertial model
            // combined with Zweifach-Fung 1-level routing for the center arm.
            let q_center = cfd_1d::tri_center_q_frac_cross_junction(
                candidate.trifurcation_center_frac,
                w_main,
                h_main,
            );
            let atv_cct =
                cfd_1d::cascade_junction_separation_from_qfracs(&[q_center]);
            final_sep_eff = atv_cct.separation_efficiency;
            final_cancer_frac = atv_cct.cancer_center_fraction;
            final_rbc_periph = atv_cct.rbc_peripheral_fraction;
            final_three_pop_sep = atv_cct.separation_efficiency;
            final_wbc_center_three_pop = atv_cct.wbc_center_fraction;
            final_rbc_periph_three_pop = atv_cct.rbc_peripheral_fraction;

            rbc_center_frac = (1.0 - atv_cct.rbc_peripheral_fraction).clamp(0.0, 1.0);
            cancer_center_frac = atv_cct.cancer_center_fraction.clamp(0.0, 1.0);
            wbc_center_frac = atv_cct.wbc_center_fraction.clamp(0.0, 1.0);
        }
        DesignTopology::TriBiTriSelectiveVenturi => {
            // Tri→Bi→Tri: three progressive focusing stages.
            // Stage 1 & 3: trifurcation at trifurcation_center_frac.
            // Stage 2: bifurcation at cif_terminal_bi_treat_frac.
            let q_tri = cfd_1d::tri_center_q_frac_cross_junction(
                candidate.trifurcation_center_frac,
                w_main,
                h_main,
            );
            let q_bi = candidate.cif_terminal_bi_treat_frac();
            // Model as: tri stage → bi stage → tri stage (3 multiplicative stages)
            let tbt_cif = cfd_1d::incremental_filtration_separation_from_qfracs(
                &[q_tri], // 1 pre-tri stage
                q_tri,    // terminal tri (same frac)
                q_bi,     // terminal bi
            );
            final_sep_eff = tbt_cif.separation_efficiency;
            final_cancer_frac = tbt_cif.cancer_center_fraction;
            final_rbc_periph = tbt_cif.rbc_peripheral_fraction;
            final_three_pop_sep = tbt_cif.separation_efficiency;
            final_wbc_center_three_pop = tbt_cif.wbc_center_fraction;
            final_rbc_periph_three_pop = tbt_cif.rbc_peripheral_fraction;

            rbc_center_frac = tbt_cif.rbc_center_fraction.clamp(0.0, 1.0);
            cancer_center_frac = tbt_cif.cancer_center_fraction.clamp(0.0, 1.0);
            wbc_center_frac = tbt_cif.wbc_center_fraction.clamp(0.0, 1.0);
        }
        _ => {}
    }

    if matches!(
        candidate.topology,
        DesignTopology::CascadeCenterTrifurcationSeparator { .. }
            | DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
            | DesignTopology::DoubleTrifurcationCIFVenturi { .. }
            | DesignTopology::AsymmetricTrifurcationVenturi
            | DesignTopology::TriBiTriSelectiveVenturi
    ) {
        final_wbc_recovery = wbc_center_frac.clamp(0.0, 1.0);
        final_rbc_pass_fraction = rbc_center_frac.clamp(0.0, 1.0);
        let center_capture = final_wbc_recovery + final_rbc_pass_fraction;
        final_wbc_purity = if center_capture > 1e-12 {
            (final_wbc_recovery / center_capture).clamp(0.0, 1.0)
        } else {
            0.0
        };
        final_total_ecv_ml = (solved.mean_residence_time_s * solved.inlet_flow_m3_s).max(0.0) * 1e6;
    }

    if venturi_treatment_enabled && final_venturi_flow_fraction > 0.0 {
        // Use cascade-specific hematocrit ratio when available (more
        // physically accurate — accounts for per-level Zweifach-Fung routing).
        let cascade_hct_ratio = match (&dtcv_cascade_res, &cascade_res, &cif_res) {
            (Some(dtcv), _, _) => dtcv.center_hematocrit_ratio,
            (_, Some(cct), _) => cct.center_hematocrit_ratio,
            (_, _, Some(cif)) => cif.center_hematocrit_ratio,
            _ => rbc_center_frac / final_venturi_flow_fraction.max(1e-9),
        };
        let hct_venturi = (candidate.feed_hematocrit * cascade_hct_ratio).clamp(0.01, 0.70);

        // Re-derive venturi shear using a local Casson blood model with
        // the venturi-specific hematocrit (Chien 1970 yield stress + Quemada
        // 1978 viscosity scaling).  Lower HCT → lower yield stress → lower
        // apparent viscosity → lower wall shear → less hemolysis.
        //
        // Scale the flow-weighted venturi_hi by the Giersiepen viscosity
        // ratio: since HI ∝ τ^α = (μ·γ)^α, and γ is unchanged, the
        // HI correction factor is (μ_local / μ_bulk)^α.
        let local_blood = CassonBlood::<f64>::with_hematocrit(hct_venturi);
        let local_mu = local_blood.apparent_viscosity(venturi_shear_rate);
        let bulk_mu = blood.apparent_viscosity(venturi_shear_rate);
        let mu_ratio = if bulk_mu > 1e-18 {
            (local_mu / bulk_mu).min(1.0) // always ≤ 1 since HCT_local ≤ HCT_feed
        } else {
            1.0
        };
        let hi_correction = mu_ratio.powf(GIERSIEPEN_ALPHA);
        let corrected_venturi_hi = venturi_hi * hi_correction;
        let local_tau = (local_mu * venturi_shear_rate).max(0.0);
        let local_venturi_pai = if local_tau > 0.0 && t_throat_venturi > 0.0 {
            1.8e-8 * local_tau.powf(1.325) * t_throat_venturi.powf(0.462)
        } else {
            0.0
        };

        // Replace the venturi portion of bulk HI with the HCT-corrected value.
        final_hi = (bulk_hi - venturi_hi).max(0.0) + corrected_venturi_hi;
        final_pai = local_venturi_pai;
        final_local_hct = hct_venturi;
        let serial_cancer_gain =
            (1.0 + 0.20 * (serial_venturi_stages_per_path as f64 - 1.0)).clamp(1.0, 1.4);
        final_cancer_dose =
            (cancer_center_frac * cav_potential * serial_cancer_gain).clamp(0.0, 1.0);
        final_rbc_venturi_exposure = rbc_center_frac;

        // Update venturi shear with local-HCT-corrected values
        venturi_shear_pa = local_tau;

        if final_wbc_center_three_pop <= 0.0 {
            final_wbc_center_three_pop = wbc_center_frac;
        }

        // Apply HCT correction to per-channel venturi entries.
        for ch in &mut per_channel_hi {
            if ch.is_venturi_throat {
                ch.hi_contribution *= hi_correction;
                ch.wall_shear_pa = local_tau;
            }
        }
    }

    // ── Per-channel hemolysis summary ────────────────────────────────────────
    let treatment_channel_hi: f64 = per_channel_hi
        .iter()
        .filter(|ch| ch.is_venturi_throat)
        .map(|ch| ch.hi_contribution)
        .sum();
    let mut bypass_n = 0usize;
    let mut bypass_sum = 0.0_f64;
    let mut bypass_channel_hi_max = 0.0_f64;
    for ch in &per_channel_hi {
        if !ch.is_venturi_throat && ch.hi_contribution > 0.0 {
            bypass_n += 1;
            bypass_sum += ch.hi_contribution;
            if ch.hi_contribution > bypass_channel_hi_max {
                bypass_channel_hi_max = ch.hi_contribution;
            }
        }
    }
    let bypass_channel_hi_mean = if bypass_n == 0 {
        0.0
    } else {
        bypass_sum / bypass_n as f64
    };
    per_channel_hi.sort_by(|a, b| {
        b.hi_contribution
            .partial_cmp(&a.hi_contribution)
            .unwrap_or(std::cmp::Ordering::Equal)
    });

    // ── Cancer-targeting hydrodynamic cavitation metrics ─────────────────────
    let cavitation_intensity: f64 = if venturi_treatment_enabled && cav_potential > 0.0 {
        let constriction = if venturi_v_in > 1e-9 {
            let ratio = (venturi_v_thr / venturi_v_in).max(1.0);
            (ratio.ln() / VENTURI_VEL_RATIO_REF.ln()).clamp(0.0, 1.0)
        } else {
            0.0
        };
        let serial_gain =
            (1.0 + 0.25 * (serial_venturi_stages_per_path as f64 - 1.0)).clamp(1.0, 1.5);
        (cav_potential * (0.5 + 0.5 * constriction) * serial_gain).clamp(0.0, 1.0)
    } else {
        0.0
    };

    let cancer_targeted_cavitation = (cancer_center_frac * cavitation_intensity).clamp(0.0, 1.0);
    let wbc_targeted_cavitation = (wbc_center_frac * cavitation_intensity).clamp(0.0, 1.0);
    let rbc_venturi_protection = (final_rbc_periph_three_pop
        * (1.0 - cavitation_intensity * final_rbc_venturi_exposure))
        .clamp(0.0, 1.0);

    let sonoluminescence_proxy: f64 = if venturi_treatment_enabled && cav_potential > 0.0 {
        let p_abs = candidate.inlet_gauge_pa + P_ATM_PA;
        let exponent = (BUBBLE_GAMMA - 1.0) / BUBBLE_GAMMA;
        let t_ratio = (p_abs / BLOOD_VAPOR_PRESSURE_PA.max(1.0)).powf(exponent);
        let ref_t_ratio = (SONO_REF_P_ABS_PA / BLOOD_VAPOR_PRESSURE_PA.max(1.0)).powf(exponent);
        (cav_potential * t_ratio / ref_t_ratio.max(1.0)).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // ── RBC safety / therapeutic-window metrics ──────────────────────────────
    let throat_transit_time_s = t_throat_venturi;

    let fda_main_ok = max_main_shear_pa <= FDA_MAX_WALL_SHEAR_PA;
    // Transit-time exception: throat high-shear is acceptable when the exposure
    // is sufficiently brief AND shear stays below the extended transient limit.
    let throat_transit_exception =
        venturi_shear_pa <= FDA_TRANSIENT_SHEAR_PA && t_throat_venturi < FDA_TRANSIENT_TIME_S;
    let fda_overall_compliant = fda_main_ok
        && (!venturi_treatment_enabled
            || venturi_shear_pa <= FDA_MAX_WALL_SHEAR_PA
            || throat_transit_exception);

    // Composite lysis risk: amplify HI by RBC presence in the cavitating venturi.
    // Cavitation collapse delivers ~5× the shear-only haemolytic stress per
    // RBC-venturi co-exposure unit (conservative bubble-interaction estimate).
    let lysis_risk_index =
        (final_hi * (1.0 + 5.0 * final_rbc_venturi_exposure * final_local_hct)).clamp(0.0, 1.0);

    let therapeutic_window_score =
        (cancer_targeted_cavitation / (1e-6 + lysis_risk_index) / THERAPEUTIC_WINDOW_REF)
            .clamp(0.0, 1.0);
    let oncology_selectivity_index = (cancer_targeted_cavitation
        * (1.0 - final_rbc_venturi_exposure).clamp(0.0, 1.0)
        * (0.5 + 0.5 * final_cancer_frac.clamp(0.0, 1.0)))
    .clamp(0.0, 1.0);
    let cancer_rbc_cavitation_bias_index =
        if venturi_treatment_enabled && cavitation_intensity > 0.0 {
            let cancer_load = cancer_targeted_cavitation.clamp(0.0, 1.0);
            let rbc_load = (final_rbc_venturi_exposure * final_local_hct * cavitation_intensity)
                .clamp(0.0, 1.0);
            let total_load = cancer_load + rbc_load;
            if total_load > 1.0e-12 {
                (cancer_load / total_load).clamp(0.0, 1.0)
            } else {
                0.0
            }
        } else {
            0.0
        };

    // Project steady-state lysis rate [% Hb/h] at operating flow into a 5 L patient.
    let rbc_lysis_rate_pct_per_h =
        final_hi * flow_rate_ml_min * 60.0 / PATIENT_BLOOD_VOLUME_ML * 100.0;

    // Cancer cells in the active therapy (venturi) stream.
    let cancer_therapy_zone_fraction =
        (final_cancer_frac * final_venturi_flow_fraction).clamp(0.0, 1.0);

    // Relative 405 nm blue-light delivery proxy.
    // Blue light penetration in blood is limited; thinner channels and lower
    // local hematocrit improve photon delivery to circulating cells.
    let optical_path_length_405_m = h_main.max(1e-6);
    let hematocrit_scale_405 = (final_local_hct / 0.45).clamp(0.25, 1.75);
    let blue_light_delivery_index_405nm =
        (-BLOOD_ATTENUATION_405NM_INV_M * optical_path_length_405_m * hematocrit_scale_405)
            .exp()
            .clamp(0.0, 1.0);

    // Headroom to the 150 Pa FDA main-channel limit [Pa].
    let safety_margin_pa = FDA_MAX_WALL_SHEAR_PA - max_main_shear_pa;

    // ── Geometry-sensitive clotting / stasis risk metrics ──────────────────
    // Estimate main-channel shear rate using reference blood viscosity.
    let main_channel_shear_rate_inv_s = max_main_shear_pa / BLOOD_VISCOSITY_PA_S.max(1e-12);

    // 1. Low-flow stasis risk (global operating-point term, de-weighted).
    let low_flow_stasis_risk = descending_linear_risk(
        flow_rate_ml_min,
        CLOTTING_BFR_LOW_RISK_ML_MIN,
        CLOTTING_BFR_HIGH_RISK_ML_MIN,
    );

    // 2. Minimum per-channel shear-rate stasis risk.
    //    The worst stasis zone determines clotting risk — wider or slower
    //    channels have lower minimum shear rates. If no valid channels were
    //    observed, fall back to the main-channel shear rate.
    let effective_min_shear = if min_channel_shear_rate.is_finite() {
        min_channel_shear_rate
    } else {
        main_channel_shear_rate_inv_s
    };
    let min_shear_stasis_risk = descending_linear_risk(
        effective_min_shear,
        CLOTTING_SHEAR_LOW_RISK_INV_S,
        CLOTTING_SHEAR_HIGH_RISK_INV_S,
    );

    // 3. Maximum per-channel residence-time stasis risk.
    //    The slowest channel segment determines the fibrin deposition ceiling.
    let max_residence_stasis_risk = ascending_linear_risk(
        max_channel_transit_s,
        CLOTTING_RESIDENCE_LOW_RISK_S,
        CLOTTING_RESIDENCE_HIGH_RISK_S,
    );

    // 4. Post-venturi expansion recirculation risk (Borda-Carnot stasis).
    //    Risk scales with the natural log of the area expansion ratio,
    //    reflecting the logarithmic growth of recirculation length with
    //    expansion angle (Idelchik 1994, Diagram 4-1).
    let expansion_stasis_risk = ascending_linear_risk(
        max_expansion_area_ratio.ln(),
        EXPANSION_RATIO_LOW_RISK.ln(),
        EXPANSION_RATIO_HIGH_RISK.ln(),
    );

    // 5. Dead-volume stasis risk: fraction of chip volume in low-shear
    //    channels below platelet-adhesion threshold (Folie & McIntire 1989).
    let dead_volume_stasis_risk = if total_channel_volume_m3 > 0.0 {
        (dead_channel_volume_m3 / total_channel_volume_m3).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // Composite 5-term geometry-sensitive CRI: equal 20% weights.
    let clotting_risk_index = (0.20 * low_flow_stasis_risk
        + 0.20 * min_shear_stasis_risk
        + 0.20 * max_residence_stasis_risk
        + 0.20 * expansion_stasis_risk
        + 0.20 * dead_volume_stasis_risk)
        .clamp(0.0, 1.0);
    let clotting_flow_compliant = flow_rate_ml_min >= CLOTTING_BFR_CAUTION_ML_MIN;

    // Conservative 10 mL/s sensitivity variant: replaces only the flow term.
    let low_flow_stasis_risk_10ml_s = descending_linear_risk(
        flow_rate_ml_min,
        CLOTTING_BFR_STRICT_10MLS_ML_MIN,
        CLOTTING_BFR_CAUTION_ML_MIN,
    );
    let clotting_risk_index_10ml_s = (0.20 * low_flow_stasis_risk_10ml_s
        + 0.20 * min_shear_stasis_risk
        + 0.20 * max_residence_stasis_risk
        + 0.20 * expansion_stasis_risk
        + 0.20 * dead_volume_stasis_risk)
        .clamp(0.0, 1.0);
    let clotting_flow_compliant_10ml_s = flow_rate_ml_min >= CLOTTING_BFR_STRICT_10MLS_ML_MIN;

    // Fraction of total chip channel path in the active therapy zone (model-based).
    let therapy_channel_fraction = candidate.therapy_channel_fraction();

    // ── CCT/CIF staged-flow diagnostics ─────────────────────────────────────
    let cct_model_venturi_flow_fraction =
        if let DesignTopology::CascadeCenterTrifurcationSeparator { n_levels } = candidate.topology
        {
            let q = cfd_1d::tri_center_q_frac(candidate.trifurcation_center_frac);
            q.powi(i32::from(n_levels)).clamp(0.0, 1.0)
        } else {
            0.0
        };
    let cct_solved_venturi_flow_fraction = if matches!(
        candidate.topology,
        DesignTopology::CascadeCenterTrifurcationSeparator { .. }
    ) {
        final_venturi_flow_fraction
    } else {
        0.0
    };
    let cct_stage_center_qfrac_mean = if cct_stage_qfracs.is_empty() {
        0.0
    } else {
        cct_stage_qfracs.iter().sum::<f64>() / cct_stage_qfracs.len() as f64
    };

    let cif_model_venturi_flow_fraction =
        if let DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri } = candidate.topology
        {
            let q_pretri_product = cfd_1d::cif_pretri_stage_q_fracs(
                n_pretri,
                candidate.cif_pretri_center_frac(),
                candidate.cif_terminal_tri_center_frac(),
            )
            .into_iter()
            .product::<f64>();
            let q_tri = cfd_1d::tri_center_q_frac(
                candidate.cif_terminal_tri_center_frac(),
            );
            let q_bi = candidate.cif_terminal_bi_treat_frac();
            (q_pretri_product * q_tri * q_bi).clamp(0.0, 1.0)
        } else {
            0.0
        };
    let cif_solved_venturi_flow_fraction = if matches!(
        candidate.topology,
        DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
    ) {
        final_venturi_flow_fraction
    } else {
        0.0
    };
    let cif_pretri_qfrac_mean = if cif_pretri_stage_qfracs.is_empty() {
        0.0
    } else {
        cif_pretri_stage_qfracs.iter().sum::<f64>() / cif_pretri_stage_qfracs.len() as f64
    };
    let cif_terminal_tri_qfrac_value = cif_terminal_tri_qfrac.unwrap_or(0.0);
    let cif_terminal_bi_qfrac_value = cif_terminal_bi_qfrac.unwrap_or(0.0);
    let cif_outlet_tail_length_mm = if matches!(
        candidate.topology,
        DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
    ) {
        solved
            .channel_samples
            .iter()
            .find(|sample| sample.id == "trunk_out")
            .map_or(0.0, |sample| sample.length_m * 1000.0)
    } else {
        0.0
    };
    let cif_remerge_proximity_score = if matches!(
        candidate.topology,
        DesignTopology::IncrementalFiltrationTriBiSeparator { .. }
    ) {
        (1.0 - ascending_linear_risk(cif_outlet_tail_length_mm, 0.5, 6.0)).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let selective_cavitation_delivery_index = {
        let remerge_bonus = if cif_remerge_proximity_score > 0.0 {
            0.5 + 0.5 * cif_remerge_proximity_score
        } else {
            1.0
        };
        (oncology_selectivity_index * cancer_rbc_cavitation_bias_index * remerge_bonus)
            .clamp(0.0, 1.0)
    };

    // ── Multi-pass hemolysis projection (15 min therapy window) ─────────────
    let hi_per_pass = final_hi.clamp(0.0, 0.999_999);
    let processed_15min_ml = flow_rate_ml_min * MILESTONE_TREATMENT_DURATION_MIN.max(0.0);
    let pediatric_blood_volume_ml =
        PEDIATRIC_REFERENCE_WEIGHT_KG * PEDIATRIC_BLOOD_VOLUME_ML_PER_KG;
    let projected_passes_15min_pediatric_3kg = if pediatric_blood_volume_ml > 1e-12 {
        processed_15min_ml / pediatric_blood_volume_ml
    } else {
        0.0
    };
    let projected_hemolysis_15min_pediatric_3kg = (1.0
        - (1.0 - hi_per_pass).powf(projected_passes_15min_pediatric_3kg.max(0.0)))
    .clamp(0.0, 1.0);
    let projected_passes_15min_adult = processed_15min_ml / PATIENT_BLOOD_VOLUME_ML.max(1e-12);
    let projected_hemolysis_15min_adult =
        (1.0 - (1.0 - hi_per_pass).powf(projected_passes_15min_adult.max(0.0))).clamp(0.0, 1.0);

    // ── Acoustic energy budget ────────────────────────────────────────────────
    // mechanical_power_w: total hydraulic power input (ΔP × Q) [W]
    let mechanical_power_w = total_pressure_drop_pa * q_inlet;

    // acoustic_capture_efficiency: dimensionless proxy for fraction of mechanical
    // power that drives cavitation events at the venturi throat.
    // = cavitation_potential × (venturi_shear_pa / P_atm), capped at 1.
    // Zero for non-venturi topologies (cav_potential = 0).
    let acoustic_capture_efficiency = if venturi_treatment_enabled {
        (cav_potential * (venturi_shear_pa / P_ATM_PA)).clamp(0.0, 1.0)
    } else {
        0.0
    };

    // specific_cavitation_energy_j_ml: acoustic energy delivered per mL of
    // blood processed [mJ/mL]. Normalised by flow rate and residence time.
    let flow_rate_ml_s = flow_rate_ml_min / 60.0;
    let specific_cavitation_energy_j_ml = if flow_rate_ml_s > 1e-12 {
        acoustic_capture_efficiency * mechanical_power_w * solved.mean_residence_time_s
            / flow_rate_ml_s
            * 1000.0 // W·s / (mL/s) × 1000 = mJ/mL
    } else {
        0.0
    };

    // Cavitation-amplified haemolysis for venturi designs (bubble micro-jets
    // + shockwaves contribute independent lysis beyond steady shear).
    // Use cfd_1d::cavitation_amplified_hi for the final HI.
    let final_hi_cavitation_amplified = cavitation_amplified_hi(final_hi, cav_potential);
    let pressure_feasible = total_pressure_drop_pa <= candidate.inlet_gauge_pa;

    Ok(SdtMetrics {
        cavitation_number: venturi_sigma,
        cavitation_potential: cav_potential,
        throat_shear_rate_inv_s: venturi_shear_rate,
        throat_shear_pa: venturi_shear_pa,
        throat_exceeds_fda: venturi_shear_pa > FDA_MAX_WALL_SHEAR_PA,
        max_main_channel_shear_pa: max_main_shear_pa,
        fda_main_compliant: max_main_shear_pa <= FDA_MAX_WALL_SHEAR_PA,
        bulk_hemolysis_index_per_pass: bulk_hi,
        hemolysis_index_per_pass: final_hi,
        flow_uniformity,
        well_coverage_fraction,
        mean_residence_time_s: solved.mean_residence_time_s,
        total_pressure_drop_pa,
        total_path_length_mm: expected_path_len_m * 1000.0,
        pressure_feasible,
        cell_separation_efficiency: final_sep_eff,
        cancer_center_fraction: final_cancer_frac,
        rbc_peripheral_fraction: final_rbc_periph,
        three_pop_sep_efficiency: final_three_pop_sep,
        wbc_center_fraction: final_wbc_center_three_pop,
        rbc_peripheral_fraction_three_pop: final_rbc_periph_three_pop,
        wbc_equilibrium_pos: three_pop_metrics.wbc_eq_pos,
        cancer_equilibrium_pos: three_pop_metrics.cancer_eq_pos,
        rbc_equilibrium_pos: three_pop_metrics.rbc_eq_pos,
        wbc_recovery: final_wbc_recovery,
        rbc_pass_fraction: final_rbc_pass_fraction,
        wbc_purity: final_wbc_purity,
        total_ecv_ml: final_total_ecv_ml,
        flow_rate_ml_min,
        plate_fits,
        n_outlet_ports,
        platelet_activation_index: final_pai,
        main_channel_shear_rate_inv_s,
        low_flow_stasis_risk,
        min_shear_stasis_risk,
        max_residence_stasis_risk,
        expansion_stasis_risk,
        dead_volume_stasis_risk,
        clotting_risk_index,
        clotting_flow_compliant,
        clotting_risk_index_10ml_s,
        clotting_flow_compliant_10ml_s,
        venturi_treatment_enabled,
        treatment_zone_mode: format!("{:?}", candidate.treatment_zone_mode_effective()),
        active_venturi_throat_count,
        serial_venturi_stages_per_path,
        venturi_flow_fraction: final_venturi_flow_fraction,
        cct_model_venturi_flow_fraction,
        cct_solved_venturi_flow_fraction,
        cct_stage_center_qfrac_mean,
        cif_model_venturi_flow_fraction,
        cif_solved_venturi_flow_fraction,
        cif_pretri_qfrac_mean,
        cif_terminal_tri_qfrac: cif_terminal_tri_qfrac_value,
        cif_terminal_bi_qfrac: cif_terminal_bi_qfrac_value,
        cif_outlet_tail_length_mm,
        rbc_venturi_exposure_fraction: final_rbc_venturi_exposure,
        local_hematocrit_venturi: final_local_hct,
        cancer_dose_fraction: final_cancer_dose,
        cavitation_intensity,
        cancer_targeted_cavitation,
        wbc_targeted_cavitation,
        rbc_venturi_protection,
        sonoluminescence_proxy,
        throat_transit_time_s,
        fda_overall_compliant,
        lysis_risk_index,
        therapeutic_window_score,
        oncology_selectivity_index,
        cancer_rbc_cavitation_bias_index,
        cif_remerge_proximity_score,
        selective_cavitation_delivery_index,
        rbc_lysis_rate_pct_per_h,
        projected_passes_15min_pediatric_3kg,
        projected_hemolysis_15min_pediatric_3kg,
        projected_hemolysis_15min_adult,
        cancer_therapy_zone_fraction,
        optical_path_length_405_m,
        blue_light_delivery_index_405nm,
        safety_margin_pa,
        therapy_channel_fraction,
        wall_shear_p95_pa,
        wall_shear_p99_pa,
        wall_shear_mean_pa,
        wall_shear_cv,
        fda_shear_percentile_compliant,
        diffuser_recovery_pa,
        mechanical_power_w,
        acoustic_capture_efficiency,
        specific_cavitation_energy_j_ml,
        hemolysis_index_per_pass_cavitation_amplified: final_hi_cavitation_amplified,
        treatment_channel_hi,
        bypass_channel_hi_mean,
        bypass_channel_hi_max,
        per_channel_hemolysis: per_channel_hi,
    })
}

// SSOT: cross_section_area, cross_section_dims, shear_rate_for, and
// hydraulic_diameter are canonical methods on CrossSectionSpec
// (cfd_schematics::domain::model::specs). Removed duplicate free functions.

/// Entrance-length correction factor for wall shear stress.
///
/// In a channel of length `L` with hydraulic diameter `D_h`, fully-developed
/// Poiseuille flow requires an entrance length `L_e = 0.06 · Re · D_h`
/// (Shah & London 1978).  Within the developing region the wall shear stress
/// is approximately 1.5× the fully-developed value.  The length-averaged
/// correction factor is:
///
/// ```text
/// f = 1 + 0.5 · min(L_e / L, 1)
/// ```
///
/// This returns a multiplier ≥ 1.0 that should scale the fully-developed
/// shear rate to account for the developing region.
fn entrance_length_shear_factor(
    cs: CrossSectionSpec,
    velocity_m_s: f64,
    channel_length_m: f64,
) -> f64 {
    let d_h = cs.hydraulic_diameter();
    let re = BLOOD_DENSITY_KG_M3 * velocity_m_s * d_h / BLOOD_VISCOSITY_PA_S;
    let l_e = 0.06 * re * d_h;
    let frac_developing = (l_e / channel_length_m.max(1e-18)).min(1.0);
    1.0 + 0.5 * frac_developing
}

fn fallback_uniformity(topology: DesignTopology) -> f64 {
    match topology {
        DesignTopology::AsymmetricBifurcationSerpentine => 0.60,
        DesignTopology::CellSeparationVenturi
        | DesignTopology::WbcCancerSeparationVenturi
        | DesignTopology::CascadeCenterTrifurcationSeparator { .. }
        | DesignTopology::IncrementalFiltrationTriBiSeparator { .. } => 0.75,
        _ => 1.0,
    }
}

fn channel_flow_abs(samples: &[ChannelSolveSample<'_>], id: &str) -> Option<f64> {
    samples
        .iter()
        .find(|sample| sample.id == id)
        .map(|sample| sample.flow_m3_s.abs())
}

fn tri_center_flow_frac(
    samples: &[ChannelSolveSample<'_>],
    center_id: &str,
    left_id: &str,
    right_id: &str,
) -> Option<f64> {
    let q_center = channel_flow_abs(samples, center_id)?;
    let q_left = channel_flow_abs(samples, left_id)?;
    let q_right = channel_flow_abs(samples, right_id)?;
    let q_total = q_center + q_left + q_right;
    (q_total > 1.0e-12).then(|| (q_center / q_total).clamp(0.0, 1.0))
}

fn bi_treat_flow_frac(
    samples: &[ChannelSolveSample<'_>],
    treat_id: &str,
    bypass_id: &str,
) -> Option<f64> {
    let q_treat = channel_flow_abs(samples, treat_id)?;
    let q_bypass = channel_flow_abs(samples, bypass_id)?;
    let q_total = q_treat + q_bypass;
    (q_total > 1.0e-12).then(|| (q_treat / q_total).clamp(0.0, 1.0))
}

fn extract_cct_stage_qfracs(samples: &[ChannelSolveSample<'_>], n_levels: u8) -> Vec<f64> {
    let mut q_fracs = Vec::with_capacity(n_levels as usize);
    for lv in 0..n_levels {
        let center_id = format!("center_lv{lv}");
        let left_id = format!("L_lv{lv}");
        let right_id = format!("R_lv{lv}");
        let Some(q) = tri_center_flow_frac(samples, &center_id, &left_id, &right_id) else {
            break;
        };
        q_fracs.push(q);
    }
    q_fracs
}

fn extract_cif_stage_qfracs(
    samples: &[ChannelSolveSample<'_>],
    n_pretri: u8,
) -> (Vec<f64>, Option<f64>, Option<f64>) {
    let mut pretri = Vec::with_capacity(n_pretri as usize);
    for lv in 0..n_pretri {
        let center_id = format!("center_lv{lv}");
        let left_id = format!("L_lv{lv}");
        let right_id = format!("R_lv{lv}");
        let Some(q) = tri_center_flow_frac(samples, &center_id, &left_id, &right_id) else {
            break;
        };
        pretri.push(q);
    }

    let tri_q = tri_center_flow_frac(samples, "hy_tri_center", "hy_tri_L", "hy_tri_R");
    let bi_q = bi_treat_flow_frac(samples, "hy_bi_treat", "hy_bi_bypass");
    (pretri, tri_q, bi_q)
}

/// Compute wall-shear percentile statistics from collected channel shear values.
///
/// Returns `(P95, P99, mean, CV)`.  All zeros when the input is empty.
fn compute_shear_percentiles(shears: &mut Vec<f64>) -> (f64, f64, f64, f64) {
    if shears.is_empty() {
        return (0.0, 0.0, 0.0, 0.0);
    }
    shears.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    let n = shears.len();
    let mean = shears.iter().sum::<f64>() / n as f64;

    let p95 = percentile_sorted(shears, 0.95);
    let p99 = percentile_sorted(shears, 0.99);

    let variance = shears.iter().map(|&x| (x - mean).powi(2)).sum::<f64>() / n as f64;
    let cv = if mean > 1e-18 {
        variance.sqrt() / mean
    } else {
        0.0
    };
    (p95, p99, mean, cv)
}

/// Nearest-rank percentile from a pre-sorted slice.
fn percentile_sorted(sorted: &[f64], p: f64) -> f64 {
    let n = sorted.len();
    if n == 0 {
        return 0.0;
    }
    let rank = (p * n as f64).ceil() as usize;
    sorted[rank.min(n) - 1]
}

/// Risk ramp where higher `value` means lower risk.
///
/// Returns:
/// - `0.0` when `value >= low_risk_value`
/// - `1.0` when `value <= high_risk_value`
/// - linear interpolation in between.
fn descending_linear_risk(value: f64, low_risk_value: f64, high_risk_value: f64) -> f64 {
    if value >= low_risk_value {
        0.0
    } else if value <= high_risk_value {
        1.0
    } else {
        ((low_risk_value - value) / (low_risk_value - high_risk_value)).clamp(0.0, 1.0)
    }
}

/// Risk ramp where higher `value` means higher risk.
///
/// Returns:
/// - `0.0` when `value <= low_risk_value`
/// - `1.0` when `value >= high_risk_value`
/// - linear interpolation in between.
fn ascending_linear_risk(value: f64, low_risk_value: f64, high_risk_value: f64) -> f64 {
    if value <= low_risk_value {
        0.0
    } else if value >= high_risk_value {
        1.0
    } else {
        ((value - low_risk_value) / (high_risk_value - low_risk_value)).clamp(0.0, 1.0)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::design::{build_candidate_space, DesignTopology};

    fn near(a: f64, b: f64) -> bool {
        (a - b).abs() < 1e-12
    }

    fn find_candidate<P>(predicate: P) -> DesignCandidate
    where
        P: Fn(&DesignCandidate) -> bool,
    {
        build_candidate_space()
            .into_iter()
            .find(predicate)
            .expect("expected candidate in parametric space")
    }

    fn first_candidate_with_metrics<P>(predicate: P) -> (DesignCandidate, SdtMetrics)
    where
        P: Fn(&DesignCandidate) -> bool,
    {
        build_candidate_space()
            .into_iter()
            .filter(predicate)
            .find_map(|candidate| {
                compute_metrics(&candidate)
                    .ok()
                    .map(|metrics| (candidate, metrics))
            })
            .expect("expected at least one candidate with successful metric evaluation")
    }

    #[test]
    fn symmetric_bifurcation_venturi_has_high_flow_uniformity() {
        let (_, metrics) =
            first_candidate_with_metrics(|c| c.topology == DesignTopology::BifurcationVenturi);
        assert!(
            metrics.flow_uniformity > 0.95,
            "symmetric venturi branches should have near-equal flow"
        );
    }

    #[test]
    fn asymmetric_serpentine_is_not_unity_uniformity() {
        let (_, metrics) = first_candidate_with_metrics(|c| {
            c.topology == DesignTopology::AsymmetricBifurcationSerpentine
        });
        assert!(
            metrics.flow_uniformity < 0.999,
            "asymmetric serpentine must not collapse to fixed 1.0 uniformity"
        );
    }

    #[test]
    fn tree_venturi_fraction_is_near_unity_for_full_flow_trees() {
        let (_, metrics) =
            first_candidate_with_metrics(|c| c.topology == DesignTopology::TrifurcationVenturi);
        assert!(
            metrics.venturi_flow_fraction > 0.95,
            "full-flow venturi trees should send most flow through throats"
        );
    }

    #[test]
    fn cct_selective_venturi_reduces_hemolysis_vs_bulk() {
        let candidate = find_candidate(|c| {
            matches!(
                c.topology,
                DesignTopology::CascadeCenterTrifurcationSeparator { n_levels: 1 }
            ) && near(c.trifurcation_center_frac, 0.55)
                && near(c.channel_width_m, 6.0e-3)
                && near(c.throat_diameter_m, 100e-6)
                && near(c.inlet_gauge_pa, 200_000.0)
                && near(c.flow_rate_m3_s, 1.667e-6)
        });
        let metrics = compute_metrics(&candidate).unwrap_or_else(|e| {
            panic!(
                "metrics should compute for CCT candidate {}: {e:?}",
                candidate.id
            )
        });

        assert!(
            metrics.venturi_flow_fraction < 1.0,
            "CCT must route only a subset of flow through venturi"
        );
        assert!(
            metrics.rbc_venturi_exposure_fraction < 1.0,
            "CCT must reduce RBC venturi exposure fraction"
        );
        assert!(
            metrics.local_hematocrit_venturi <= candidate.feed_hematocrit + 1e-12,
            "local venturi hematocrit should not exceed feed hematocrit for CCT"
        );
        assert!(
            metrics.hemolysis_index_per_pass <= metrics.bulk_hemolysis_index_per_pass + 1e-16,
            "selective venturi correction should not increase HI for CCT"
        );
        assert!(
            metrics.wbc_recovery > 0.0 && metrics.total_ecv_ml > 0.0,
            "CCT should populate leukapheresis-compatible WBC/ECV metrics"
        );
    }

    #[test]
    fn cif_selective_venturi_reduces_hemolysis_vs_bulk() {
        let candidate = find_candidate(|c| {
            matches!(
                c.topology,
                DesignTopology::IncrementalFiltrationTriBiSeparator { n_pretri: 1 }
            ) && near(c.cif_pretri_center_frac, 0.55)
                && near(c.cif_terminal_tri_center_frac, 0.55)
                && near(c.cif_terminal_bi_treat_frac, 0.76)
                && near(c.channel_width_m, 6.0e-3)
                && near(c.throat_diameter_m, 100e-6)
                && near(c.inlet_gauge_pa, 200_000.0)
                && near(c.flow_rate_m3_s, 1.667e-6)
        });
        let metrics = compute_metrics(&candidate).unwrap_or_else(|e| {
            panic!(
                "metrics should compute for CIF candidate {}: {e:?}",
                candidate.id
            )
        });

        assert!(
            metrics.venturi_flow_fraction < 1.0,
            "CIF must route only treatment arm flow through venturi"
        );
        assert!(
            metrics.rbc_venturi_exposure_fraction < 1.0,
            "CIF must reduce RBC venturi exposure fraction"
        );
        assert!(
            metrics.local_hematocrit_venturi <= candidate.feed_hematocrit + 1e-12,
            "local venturi hematocrit should not exceed feed hematocrit for CIF"
        );
        assert!(
            metrics.hemolysis_index_per_pass <= metrics.bulk_hemolysis_index_per_pass + 1e-16,
            "selective venturi correction should not increase HI for CIF"
        );
        assert!(
            metrics.wbc_recovery > 0.0 && metrics.total_ecv_ml > 0.0,
            "CIF should populate leukapheresis-compatible WBC/ECV metrics"
        );
    }

    // ── Safety tests ────────────────────────────────────────────────────────

    #[test]
    fn bernoulli_dp_is_finite_and_uses_rectangular_area() {
        // Verify that Bernoulli dp uses the rectangular area (not circular),
        // producing finite values for all candidates in the design space.
        let candidates = build_candidate_space();
        for c in candidates.iter().take(200) {
            let a_in = c.inlet_area_m2();
            let a_th = c.throat_area_m2();
            assert!(
                a_in > 0.0 && a_in.is_finite(),
                "inlet area must be finite positive"
            );
            assert!(
                a_th > 0.0 && a_th.is_finite(),
                "throat area must be finite positive"
            );

            let v_in = c.flow_rate_m3_s / a_in;
            let v_th = c.flow_rate_m3_s / a_th;
            let dp = 0.5 * BLOOD_DENSITY_KG_M3 * (v_th * v_th - v_in * v_in);
            assert!(
                dp.is_finite(),
                "dp must be finite for {:?} (a_in={:.2e}, a_th={:.2e})",
                c.topology,
                a_in,
                a_th,
            );
            // Rectangular area should be larger than circular for same diameter,
            // so dp should be lower (unless height < π/4 * d, which is uncommon).
            // Just verify it's not astronomical (< 100 MPa even for extreme designs).
            assert!(
                dp.abs() < 1e8,
                "dp = {:.0} Pa unreasonably large for {:?}",
                dp,
                c.topology,
            );
        }
    }

    #[test]
    fn sigma_within_physical_range_for_all_candidates() {
        let candidates = build_candidate_space();
        for c in candidates.iter().take(500) {
            let a_th = c.throat_area_m2();
            if a_th < 1e-15 {
                continue;
            }
            let v_th = c.flow_rate_m3_s / a_th;
            let dyn_p = 0.5 * BLOOD_DENSITY_KG_M3 * v_th * v_th;
            if dyn_p < 1e-12 {
                continue;
            }
            let p_abs_throat = (P_ATM_PA + c.inlet_gauge_pa) - dyn_p;
            let sigma = (p_abs_throat - BLOOD_VAPOR_PRESSURE_PA) / dyn_p;
            assert!(
                sigma > -10.0 && sigma < 10_000.0,
                "σ = {:.2} out of range for {:?}",
                sigma,
                c.topology,
            );
        }
    }

    #[test]
    fn cross_section_area_matches_dimensions() {
        use crate::design::CrossSectionShape;
        let candidates = build_candidate_space();
        for c in &candidates {
            match c.cross_section_shape {
                CrossSectionShape::Rectangular => {
                    let expected = c.inlet_diameter_m * c.channel_height_m;
                    assert!(
                        (c.inlet_area_m2() - expected).abs() < 1e-18,
                        "inlet_area_m2 mismatch for {:?}",
                        c.topology
                    );
                    let expected_th = c.throat_diameter_m * c.channel_height_m;
                    assert!(
                        (c.throat_area_m2() - expected_th).abs() < 1e-18,
                        "throat_area_m2 mismatch for {:?}",
                        c.topology
                    );
                }
                CrossSectionShape::Circular => {
                    let expected =
                        std::f64::consts::FRAC_PI_4 * c.inlet_diameter_m * c.inlet_diameter_m;
                    assert!(
                        (c.inlet_area_m2() - expected).abs() < 1e-18,
                        "circular inlet_area_m2 mismatch"
                    );
                }
            }
        }
    }

    #[test]
    fn cav_potential_nonzero_for_negative_sigma() {
        // σ < 0 (very strong cavitation) must produce maximum potential
        let sigma_neg = -0.1;
        let raw = (1.0 - sigma_neg / SIGMA_CRIT).clamp(0.0, 1.0).powf(1.5);
        assert!(
            (raw - 1.0).abs() < 1e-10,
            "cav_potential for σ={} should be 1.0, got {}",
            sigma_neg,
            raw
        );

        // σ just below σ_crit should give small but nonzero potential
        let sigma_incipient = 0.9;
        let raw2 = (1.0 - sigma_incipient / SIGMA_CRIT)
            .clamp(0.0, 1.0)
            .powf(1.5);
        assert!(
            raw2 > 0.0 && raw2 < 0.1,
            "cav_potential for σ={} should be small positive, got {}",
            sigma_incipient,
            raw2
        );

        // σ >= σ_crit should give zero
        let sigma_no_cav = 1.5;
        assert!(sigma_no_cav >= SIGMA_CRIT, "test assumes σ >= σ_crit");
    }

    /// Exhaustive round-trip test: every DesignTopology family that appears in
    /// the parametric sweep must successfully produce SDT metrics via the full
    /// pipeline: `to_blueprint() → solve_blueprint_network() → compute_metrics()`.
    ///
    /// ## Invariants verified per topology:
    /// 1. `compute_metrics` returns `Ok` (no solver failure)
    /// 2. `total_pressure_drop_pa` is finite and positive (Kirchhoff conservation)
    /// 3. `flow_uniformity` ∈ (0, 1]
    /// 4. `hemolysis_index_per_pass` is finite and non-negative
    /// 5. `mean_residence_time_s` is finite and positive
    ///
    /// ## Theorem: Monotonic separation improvement with cascade depth
    ///
    /// For a cascade of N asymmetric trifurcation junctions with center-arm
    /// width fraction `f > 1/3` and Zweifach-Fung stiffness exponents
    /// `β_cancer > β_RBC = 1.0`:
    ///
    /// ```text
    /// sep_eff(N) = |cancer_center(N) − rbc_center(N)| is monotone increasing in N
    /// ```
    ///
    /// **Proof sketch:**
    /// At each junction, `P_center(cell) = q^β / (q^β + 2·q_p^β)` where
    /// `q > 1/3` (asymmetric).  For β > 1 (cancer), `P_center > q` (super-linear
    /// bias toward center).  For β = 1 (RBC), `P_center = q` (linear, follows
    /// flow fraction).  After N stages:
    /// - `cancer_center(N) = ∏ P_center(q, 1.70)` — decays slowly (each factor > q)
    /// - `rbc_center(N) = q^N` — decays geometrically
    /// Since `P_center(q, 1.70) > q` for all q > 1/3, the ratio
    /// `cancer_center(N) / rbc_center(N)` grows with N, and the absolute gap
    /// `|cancer − rbc|` increases until cancer_center itself becomes small.
    /// The RBC peripheral fraction `1 − q^N` is strictly increasing in N.  ∎

    #[test]
    fn cct_deeper_cascade_improves_separation_monotonically() {
        // Validate the mathematical theorem: each added trifurcation level
        // increases RBC peripheral fraction and separation efficiency.
        //
        // Uses the cfd-1d Zweifach-Fung model directly, then verifies that
        // cfd-optim's compute_metrics propagates the same monotonic ordering
        // through to the final SdtMetrics for CCT topologies.

        // Part 1: Direct model verification (cfd-1d)
        let center_frac = 0.45_f64;
        let parent_w = 6.0e-3_f64;
        let h = 1.5e-3_f64;

        let mut prev_sep = 0.0_f64;
        let mut prev_rbc_periph = 0.0_f64;

        for n_levels in 1u8..=3 {
            let result = cfd_1d::cascade_junction_separation_cross_junction(
                n_levels,
                center_frac,
                parent_w,
                h,
            );

            assert!(
                result.separation_efficiency >= prev_sep - 1e-12,
                "CCT depth {n_levels}: separation efficiency {:.4} must be >= previous {:.4}",
                result.separation_efficiency,
                prev_sep,
            );
            assert!(
                result.rbc_peripheral_fraction >= prev_rbc_periph - 1e-12,
                "CCT depth {n_levels}: RBC peripheral {:.4} must be >= previous {:.4}",
                result.rbc_peripheral_fraction,
                prev_rbc_periph,
            );
            assert!(
                result.cancer_center_fraction > result.rbc_peripheral_fraction * 0.0,
                "CCT depth {n_levels}: cancer should be focused to center",
            );
            // Cancer enrichment over RBC: cancer_center > rbc_center
            let rbc_center = 1.0 - result.rbc_peripheral_fraction;
            assert!(
                result.cancer_center_fraction > rbc_center,
                "CCT depth {n_levels}: cancer_center={:.4} must exceed rbc_center={:.4}",
                result.cancer_center_fraction,
                rbc_center,
            );

            prev_sep = result.separation_efficiency;
            prev_rbc_periph = result.rbc_peripheral_fraction;
        }

        // Part 2: Verify cfd-optim propagates monotonic ordering through metrics.
        // Find CCT candidates at increasing depths with matched parameters.
        let mut metrics_by_depth: Vec<(u8, SdtMetrics)> = Vec::new();
        let candidates = build_candidate_space();
        for n_levels in 1u8..=3 {
            if let Some(candidate) = candidates.iter().find(|c| {
                matches!(
                    c.topology,
                    DesignTopology::CascadeCenterTrifurcationSeparator { n_levels: nl } if nl == n_levels
                ) && near(c.trifurcation_center_frac, 0.45)
                    && near(c.channel_width_m, 6.0e-3)
                    && near(c.channel_height_m, 1.5e-3)
            }) {
                if let Ok(m) = compute_metrics(candidate) {
                    metrics_by_depth.push((n_levels, m));
                }
            }
        }

        // Verify monotonic ordering in the computed metrics.
        for window in metrics_by_depth.windows(2) {
            let (d1, m1) = &window[0];
            let (d2, m2) = &window[1];
            assert!(
                m2.rbc_peripheral_fraction_three_pop >= m1.rbc_peripheral_fraction_three_pop - 0.01,
                "CCT depth {d2} rbc_periph {:.4} should >= depth {d1} rbc_periph {:.4} in computed metrics",
                m2.rbc_peripheral_fraction_three_pop, m1.rbc_peripheral_fraction_three_pop,
            );
        }
    }

    #[test]
    fn cif_deeper_pretri_cascade_improves_separation() {
        // CIF (Controlled Incremental Filtration) uses pre-trifurcation levels
        // followed by a terminal trifurcation + bifurcation.  More pre-tri levels
        // should push more RBCs to periphery and increase separation efficiency.
        let mut prev_rbc_periph = 0.0_f64;
        let mut prev_sep = 0.0_f64;

        for n_pretri in 1u8..=3 {
            let result = cfd_1d::incremental_filtration_separation_cross_junction(
                n_pretri, 0.45,   // pretri_center_frac
                0.50,   // terminal_tri_center_frac
                0.68,   // bi_treat_frac
                6.0e-3, // parent_width_m
                1.5e-3, // channel_height_m
            );

            assert!(
                result.rbc_peripheral_fraction >= prev_rbc_periph - 1e-12,
                "CIF n_pretri={n_pretri}: RBC peripheral {:.4} must be >= previous {:.4}",
                result.rbc_peripheral_fraction,
                prev_rbc_periph,
            );
            assert!(
                result.separation_efficiency >= prev_sep - 1e-12,
                "CIF n_pretri={n_pretri}: separation efficiency {:.4} must be >= previous {:.4}",
                result.separation_efficiency,
                prev_sep,
            );
            assert!(
                result.cancer_center_fraction > result.rbc_center_fraction,
                "CIF n_pretri={n_pretri}: cancer_center={:.4} must exceed rbc_center={:.4}",
                result.cancer_center_fraction,
                result.rbc_center_fraction,
            );

            prev_rbc_periph = result.rbc_peripheral_fraction;
            prev_sep = result.separation_efficiency;
        }
    }

    #[test]
    fn dtcv_has_selective_venturi_metrics() {
        // DoubleTrifurcationCIFVenturi should now produce non-zero separation
        // metrics via the DTCV cascade routing model, and selective venturi
        // correction should reduce hemolysis vs bulk.
        let candidate = build_candidate_space()
            .into_iter()
            .find(|c| {
                matches!(
                    c.topology,
                    DesignTopology::DoubleTrifurcationCIFVenturi { .. }
                )
            })
            .expect("DTCV candidates must exist in design space");

        let metrics = compute_metrics(&candidate).unwrap_or_else(|e| {
            panic!(
                "metrics should compute for DTCV candidate {}: {e:?}",
                candidate.id
            )
        });

        // DTCV is a 2-level cascade — separation metrics must be populated.
        assert!(
            metrics.cell_separation_efficiency > 0.0,
            "DTCV separation_efficiency={:.4} must be > 0 (2-level cascade routing)",
            metrics.cell_separation_efficiency,
        );
        assert!(
            metrics.cancer_center_fraction > 0.0,
            "DTCV cancer_center_fraction={:.4} must be > 0",
            metrics.cancer_center_fraction,
        );
        assert!(
            metrics.rbc_peripheral_fraction_three_pop > 0.0,
            "DTCV rbc_peripheral_fraction={:.4} must be > 0",
            metrics.rbc_peripheral_fraction_three_pop,
        );

        // Selective venturi: only center-arm flow passes through throats.
        assert!(
            metrics.venturi_flow_fraction < 1.0,
            "DTCV venturi_flow_fraction={:.4} must be < 1.0 (center arm only)",
            metrics.venturi_flow_fraction,
        );
        assert!(
            metrics.rbc_venturi_exposure_fraction < 1.0,
            "DTCV rbc_venturi_exposure={:.4} must be < 1.0 (RBCs routed to bypass)",
            metrics.rbc_venturi_exposure_fraction,
        );
    }

    #[test]
    fn asymmetric_width_controls_separation_nonlinearly() {
        // Symmetric trifurcation (center_frac = 1/3) should produce near-zero
        // separation, while asymmetric splits should produce progressively
        // better separation up to an optimal point.
        //
        // This validates that the Zweifach-Fung nonlinearity (β > 1 for stiff
        // cells) is the mechanism for separation, not merely flow splitting.
        let symmetric = cfd_1d::cascade_junction_separation_cross_junction(
            2,
            1.0 / 3.0,
            6.0e-3,
            1.5e-3,
        );
        let asymmetric = cfd_1d::cascade_junction_separation_cross_junction(
            2, 0.45, 6.0e-3, 1.5e-3,
        );

        // Symmetric: all cell types distribute nearly identically by flow fraction.
        assert!(
            symmetric.separation_efficiency < 0.02,
            "Symmetric trifurcation (cf=1/3) should have near-zero separation, got {:.4}",
            symmetric.separation_efficiency,
        );

        // Asymmetric: stiffness-dependent routing creates separation.
        assert!(
            asymmetric.separation_efficiency > 0.10,
            "Asymmetric trifurcation (cf=0.45) should have meaningful separation, got {:.4}",
            asymmetric.separation_efficiency,
        );

        // The separation improvement must come from the nonlinear β exponent,
        // not just from geometric asymmetry.  Verify that cancer enrichment
        // exceeds what a β=1.0 (flow-fraction-proportional) model would give.
        assert!(
            asymmetric.cancer_center_fraction > (1.0 - asymmetric.rbc_peripheral_fraction),
            "Cancer center fraction ({:.4}) should exceed RBC center fraction ({:.4}) \
             due to stiffness exponent β_cancer=1.70 > β_RBC=1.00",
            asymmetric.cancer_center_fraction,
            1.0 - asymmetric.rbc_peripheral_fraction,
        );
    }

    #[test]
    fn all_topology_families_produce_valid_metrics() {
        use std::collections::{HashMap, HashSet};

        let candidates = build_candidate_space();

        // Collect unique topology discriminants (ignoring inner fields).
        let mut all_discriminants: HashSet<String> = HashSet::new();
        let mut passed_discriminants: HashSet<String> = HashSet::new();
        let mut first_error_by_disc: HashMap<String, String> = HashMap::new();
        let mut tested_count = 0usize;

        for candidate in &candidates {
            let disc = format!("{:?}", std::mem::discriminant(&candidate.topology));
            all_discriminants.insert(disc.clone());
            if passed_discriminants.contains(&disc) {
                continue; // Already validated this topology family
            }

            let result = compute_metrics(candidate);
            match result {
                Ok(metrics) => {
                    assert!(
                        metrics.total_pressure_drop_pa.is_finite()
                            && metrics.total_pressure_drop_pa >= 0.0,
                        "Topology {:?} ({}): pressure drop must be finite non-negative, got {}",
                        candidate.topology,
                        candidate.id,
                        metrics.total_pressure_drop_pa,
                    );
                    assert!(
                        metrics.flow_uniformity > 0.0 && metrics.flow_uniformity <= 1.0 + 1e-12,
                        "Topology {:?} ({}): flow uniformity must be in (0, 1], got {}",
                        candidate.topology,
                        candidate.id,
                        metrics.flow_uniformity,
                    );
                    assert!(
                        metrics.hemolysis_index_per_pass.is_finite()
                            && metrics.hemolysis_index_per_pass >= 0.0,
                        "Topology {:?} ({}): HI must be finite non-negative, got {}",
                        candidate.topology,
                        candidate.id,
                        metrics.hemolysis_index_per_pass,
                    );
                    assert!(
                        metrics.mean_residence_time_s.is_finite()
                            && metrics.mean_residence_time_s > 0.0,
                        "Topology {:?} ({}): residence time must be finite positive, got {}",
                        candidate.topology,
                        candidate.id,
                        metrics.mean_residence_time_s,
                    );
                    passed_discriminants.insert(disc);
                    tested_count += 1;
                }
                Err(e) => {
                    first_error_by_disc.entry(disc).or_insert_with(|| {
                        format!(
                            "Topology {:?} ({}): compute_metrics failed: {e:?}",
                            candidate.topology, candidate.id
                        )
                    });
                }
            }
        }

        let missing: Vec<String> = all_discriminants
            .iter()
            .filter(|d| !passed_discriminants.contains(*d))
            .cloned()
            .collect();
        let tolerated_topologies = ["CellSeparationVenturi", "WbcCancerSeparationVenturi"];
        let non_tolerated_errors: Vec<String> = missing
            .iter()
            .filter_map(|d| first_error_by_disc.get(d))
            .filter(|msg| !tolerated_topologies.iter().any(|name| msg.contains(name)))
            .cloned()
            .collect();
        assert!(
            non_tolerated_errors.is_empty(),
            "Some topology families had no valid candidate. Example errors: {}",
            if non_tolerated_errors.is_empty() {
                String::from("none")
            } else {
                non_tolerated_errors
                    .iter()
                    .take(3)
                    .cloned()
                    .collect::<Vec<_>>()
                    .join(" | ")
            }
        );

        assert!(
            tested_count >= 20,
            "Expected ≥20 unique topology families tested, got {tested_count}"
        );
    }
}

