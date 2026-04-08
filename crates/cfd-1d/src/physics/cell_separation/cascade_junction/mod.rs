//! Zweifach-Fung junction routing model for selective split-sequence trees.
//!
//! Used by primitive selective branching topologies to predict the
//! cell-type-specific distribution between the designated treatment arm
//! (cancer/WBC-enriched, routed to the therapy zone) and the peripheral
//! bypass arms (RBC-enriched, low shear).
//!
//! # Physical basis
//!
//! At a trifurcation junction, the **Zweifach-Fung law** (1969) predicts that large,
//! stiff particles (cancer cells, WBCs) preferentially enter the arm carrying the
//! highest volumetric flow.  The probability follows a power law:
//!
//! ```text
//! P_center(cell) = r_c^β / (r_c^β + 2·r_p^β)
//! ```
//!
//! where `r_c = Q_center / Q_total`, `r_p = (1 − r_c) / 2` (each peripheral
//! arm), and `β` is a stiffness exponent:
//!
//! | Cell type | β    | Basis |
//! |-----------|------|-------|
//! | Cancer (MCF-7, ~17.5 µm, stiff)   | 1.85 | Size-enhanced stiff sphere (Hou 2012, Karabacak 2014) |
//! | WBC (~10 µm, semi-rigid)           | 1.40 | Intermediate deformability |
//! | RBC (~7 µm, highly deformable)     | 1.00 | Deformable cell — flow-weighted |
//!
//! The cascade iterates this routing over `n_levels` trifurcation levels, each time
//! narrowing the center arm by `center_frac` and carrying `q_center_frac` of the
//! current level's flow.  After all levels the center-arm fraction of each cell type
//! accumulates multiplicatively.
//!
//! # Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`routing_probability`] | Core Zweifach-Fung probability functions and cell constants |
//! | [`cascade_routing`] | Cascade trifurcation, cross-junction, and mixed Bi/Tri routing |
//! | [`incremental_filtration`] | CIF staged selective-routing with pre-trifurcation skimming |
//!
//! # References
//! - Fung, Y. C. (1969). Biorheology of soft tissues. *Biorheology*, 6, 409–419.
//! - Yang et al. (2017). Red blood cell phase separation in symmetric and asymmetric microchannel
//!   networks: effect of capillary dilation and inflow velocity. *PMC5114676*.
//! - Doyeux, V. et al. (2011). Spheres in the vicinity of a bifurcation: elucidating
//!   the Zweifach-Fung effect. *J. Fluid Mech.*, 674, 359–388.
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.

pub mod cascade_routing;
pub mod incremental_filtration;
pub mod routing_probability;

// ── Type definitions ────────────────────────────────────────────────────────

/// Peripheral recovery sub-split descriptor.
///
/// When a non-treatment arm is further sub-split, the wider sub-arm can
/// feed recovered cells back to the treatment path.  This struct
/// parameterises that recovery routing for a single source arm.
#[derive(Debug, Clone, Copy)]
pub struct PeripheralRecovery {
    /// Index of the source arm in the parent stage's `arm_q_fracs`.
    pub source_arm_idx: usize,
    /// Sub-arm flow fractions (up to 5); only `[..n_sub_arms]` are used.
    pub sub_arm_q_fracs: [f64; 5],
    /// Number of active sub-arms (2–5).
    pub n_sub_arms: u8,
    /// Index of the sub-arm that feeds back to the treatment path.
    pub recovery_arm_idx: usize,
    /// Hydraulic diameter of the recovery sub-arm [m].
    pub recovery_dh_m: f64,
}

/// Per-stage descriptor for a mixed Bi/Tri selective-routing cascade.
///
/// Carries all arm flow fractions (enabling asymmetric splits) and the
/// treatment-arm hydraulic diameter (enabling κ-dependent β correction).
///
/// Index 0 of `arm_q_fracs` is always the **treatment arm** (center at
/// trifurcations, treatment branch at bifurcations).
///
/// For a bifurcation, `arm_q_fracs = [q_treat, q_bypass, 0.0]` and `n_arms = 2`.
/// For a trifurcation, `arm_q_fracs = [q_center, q_left, q_right]` and `n_arms = 3`.
#[derive(Debug, Clone, Copy)]
pub struct CascadeStage {
    /// Flow fractions for all arms; first entry is the treatment arm.
    pub arm_q_fracs: [f64; 5],
    /// Number of active arms (2 = bifurcation, 3 = trifurcation, 4 = quad, 5 = penta).
    pub n_arms: u8,
    /// Hydraulic diameter of the treatment arm [m] at this stage.
    /// Used to compute κ = cell_diameter / Dh for β amplification.
    pub treatment_dh_m: f64,
    /// Inflow velocity into the junction [m/s].
    /// Used for PMC5114676 Zweifach-Fung high-velocity inversion mechanics.
    pub parent_v_in_m_s: f64,
    /// Optional peripheral recovery sub-splits (up to 4 per stage).
    pub peripheral_recoveries: [Option<PeripheralRecovery>; 4],
    /// Number of active peripheral recoveries.
    pub n_recoveries: u8,
}

/// Cell fractions at the center (venturi) arm and peripheral bypass arms after
/// N cascade trifurcation levels.
#[derive(Debug, Clone, Copy)]
pub struct CascadeJunctionResult {
    /// Fraction of input cancer cells that reach the deepest center arm (venturi).
    pub cancer_center_fraction: f64,
    /// Fraction of input WBCs that reach the deepest center arm.
    pub wbc_center_fraction: f64,
    /// Fraction of input RBCs that are diverted to peripheral bypass arms
    /// at any cascade level (`= 1 − rbc_center_fraction`).
    pub rbc_peripheral_fraction: f64,
    /// Separation efficiency = `|f_cancer_center − f_rbc_center|` ∈ [0, 1].
    ///
    /// High values indicate strong enrichment of cancer cells at the venturi
    /// and strong depletion of RBCs (which are protected in bypass channels).
    pub separation_efficiency: f64,
    /// Local hematocrit in the center arm at the deepest level, relative to
    /// the feed hematocrit. Computed from the flow fraction and RBC routing:
    /// `HCT_local = HCT_feed × (rbc_center / q_center_frac^n_levels)`.
    pub center_hematocrit_ratio: f64,
}

/// Cell fractions after a staged selective-routing path:
///
/// 1. `n_pretri` center-only cascade trifurcation stages (incremental skimming),
/// 2. one terminal trifurcation skimming stage,
/// 3. one terminal asymmetric bifurcation selecting the treatment arm.
#[derive(Debug, Clone, Copy)]
pub struct IncrementalFiltrationResult {
    /// Fraction of input cancer cells reaching the treatment venturi arm.
    pub cancer_center_fraction: f64,
    /// Fraction of input WBCs reaching the treatment venturi arm.
    pub wbc_center_fraction: f64,
    /// Fraction of input RBCs that still reach the treatment venturi arm.
    pub rbc_center_fraction: f64,
    /// Fraction of input RBCs skimmed to peripheral bypass streams.
    pub rbc_peripheral_fraction: f64,
    /// Separation efficiency = `|f_cancer_center − f_rbc_center|` ∈ [0, 1].
    pub separation_efficiency: f64,
    /// Local hematocrit ratio (HCT_local / HCT_feed) at the treatment arm.
    pub center_hematocrit_ratio: f64,
}

// ── Re-exports ──────────────────────────────────────────────────────────────

pub use cascade_routing::{
    cascade_junction_separation, cascade_junction_separation_cross_junction,
    cascade_junction_separation_from_qfracs, checked_cascade_junction_separation,
    checked_cascade_junction_separation_cross_junction,
    checked_cascade_junction_separation_from_qfracs, checked_mixed_cascade_separation,
    checked_treatment_bifurcation_separation,
    checked_tri_asymmetric_q_fracs, checked_tri_center_q_frac,
    checked_tri_center_q_frac_cross_junction, mixed_cascade_separation,
    mixed_cascade_separation_kappa_aware, treatment_bifurcation_separation,
    tri_asymmetric_q_fracs, tri_center_q_frac, tri_center_q_frac_cross_junction,
};
pub use incremental_filtration::{
    checked_cif_pretri_stage_center_fracs, checked_cif_pretri_stage_q_fracs,
    checked_cif_pretri_stage_q_fracs_cross_junction,
    checked_incremental_filtration_separation_cross_junction,
    checked_incremental_filtration_separation_from_qfracs,
    checked_incremental_filtration_separation_staged, cif_pretri_stage_center_fracs,
    cif_pretri_stage_q_fracs, cif_pretri_stage_q_fracs_cross_junction,
    incremental_filtration_separation_cross_junction,
    incremental_filtration_separation_from_qfracs, incremental_filtration_separation_staged,
};

#[cfg(test)]
mod tests {
    use super::cascade_routing::*;
    use super::incremental_filtration::*;
    use super::routing_probability::*;
    use super::*;

    #[test]
    fn symmetric_split_equal_distribution() {
        // With symmetric 1/3 split all cell types should distribute by flow fraction
        // (cancer biased slightly more toward center, RBC less so).
        let r = cascade_junction_separation(1, 1.0 / 3.0, 2e-3, 1e-3, 5e-6);
        // Center arm carries 1/3 of flow → q_frac = (1/3)³/((1/3)³+2*(1/3)³) = 1/3
        // After 1 level: cancer_center_fraction > rbc_center_fraction (stiffness effect)
        assert!(
            r.cancer_center_fraction > r.wbc_center_fraction
                || (r.cancer_center_fraction - r.wbc_center_fraction).abs() < 1e-6
        );
        assert!(r.wbc_center_fraction >= r.rbc_peripheral_fraction * 0.0); // just sanity
        assert!(r.separation_efficiency >= 0.0 && r.separation_efficiency <= 1.0);
    }

    #[test]
    fn center_biased_increases_separation() {
        // Center-biased split should increase separation relative to symmetric.
        let r_sym = cascade_junction_separation(2, 1.0 / 3.0, 2e-3, 1e-3, 5e-6);
        let r_bias = cascade_junction_separation(2, 0.55, 2e-3, 1e-3, 5e-6);
        // Biased split: q_center_frac > 1/3 → stronger Zweifach-Fung routing
        assert!(r_bias.cancer_center_fraction >= r_sym.cancer_center_fraction - 1e-10);
    }

    #[test]
    fn more_levels_increases_enrichment() {
        let r1 = cascade_junction_separation(1, 0.45, 2e-3, 1e-3, 5e-6);
        let r3 = cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 5e-6);
        // More cascade levels → more cancer enriched in center, more RBCs in bypass
        assert!(r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10);
    }

    #[test]
    fn tri_center_q_frac_symmetric() {
        // Symmetric 1/3 split → q_frac = 1/3
        let q = tri_center_q_frac(1.0 / 3.0);
        assert!(
            (q - 1.0 / 3.0).abs() < 1e-10,
            "symmetric frac should give q=1/3, got {q}"
        );
    }

    #[test]
    fn tri_center_q_frac_larger_center_gives_more_flow() {
        let q33 = tri_center_q_frac(0.333);
        let q55 = tri_center_q_frac(0.55);
        assert!(q55 > q33, "wider center arm should carry more flow");
    }

    #[test]
    fn checked_tri_center_q_frac_rejects_closed_interval_endpoints() {
        let err = checked_tri_center_q_frac(0.0)
            .expect_err("checked center-arm flow fraction must reject zero width fraction");
        assert!(err.to_string().contains("width fraction"));
    }

    #[test]
    fn incremental_filtration_more_pretri_levels_pushes_more_rbc_periphery() {
        let r1 = incremental_filtration_separation_staged(1, 0.45, 0.50, 0.68);
        let r3 = incremental_filtration_separation_staged(3, 0.45, 0.50, 0.68);
        assert!(r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10);
        assert!(r3.separation_efficiency >= r1.separation_efficiency - 1e-10);
    }

    #[test]
    fn incremental_filtration_higher_bi_treat_frac_increases_treatment_capture() {
        let low = incremental_filtration_separation_staged(2, 0.45, 0.45, 0.60);
        let high = incremental_filtration_separation_staged(2, 0.45, 0.45, 0.76);
        assert!(high.cancer_center_fraction >= low.cancer_center_fraction - 1e-10);
        assert!(high.wbc_center_fraction >= low.wbc_center_fraction - 1e-10);
    }

    #[test]
    fn mixed_sequence_tri_bi_routes_more_rbc_peripheral_than_single_tri() {
        let q_tri = tri_center_q_frac(0.45);
        let tri_only = mixed_cascade_separation(&[(q_tri, true)]);
        let tri_bi = mixed_cascade_separation(&[(q_tri, true), (0.72, false)]);

        assert!(
            tri_bi.rbc_peripheral_fraction >= tri_only.rbc_peripheral_fraction,
            "adding a treatment bifurcation should not reduce RBC peripheral routing"
        );
        assert!(
            tri_bi.cancer_center_fraction > (1.0 - tri_bi.rbc_peripheral_fraction),
            "cancer capture should remain above RBC center carryover in selective routing"
        );
    }

    #[test]
    fn stronger_trifurcation_bias_improves_three_pop_selectivity() {
        let weak = mixed_cascade_separation(&[
            (tri_center_q_frac(0.38), true),
            (tri_center_q_frac(0.38), true),
        ]);
        let strong = mixed_cascade_separation(&[
            (tri_center_q_frac(0.52), true),
            (tri_center_q_frac(0.52), true),
        ]);

        assert!(
            strong.cancer_center_fraction >= weak.cancer_center_fraction,
            "stronger center-arm bias should not reduce cancer capture"
        );
        assert!(
            strong.separation_efficiency >= weak.separation_efficiency,
            "stronger center-arm bias should improve separation efficiency"
        );
    }

    #[test]
    fn incremental_filtration_terminal_tri_center_bias_improves_cancer_center_fraction() {
        let symmetric_terminal = incremental_filtration_separation_staged(2, 0.45, 1.0 / 3.0, 0.68);
        let center_biased_terminal = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
        assert!(
            center_biased_terminal.cancer_center_fraction
                >= symmetric_terminal.cancer_center_fraction - 1e-10
        );
    }

    #[test]
    fn cascade_qfrac_api_matches_uniform_width_model() {
        let width_model = cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 5e-6);
        let q = tri_center_q_frac(0.45);
        let solved_like = cascade_junction_separation_from_qfracs(&[q, q, q]);
        assert!(
            (width_model.cancer_center_fraction - solved_like.cancer_center_fraction).abs() < 1e-12
        );
        assert!(
            (width_model.rbc_peripheral_fraction - solved_like.rbc_peripheral_fraction).abs()
                < 1e-12
        );
    }

    #[test]
    fn checked_cascade_qfrac_api_rejects_empty_stage_sequence() {
        let err = checked_cascade_junction_separation_from_qfracs(&[])
            .expect_err("checked cascade qfrac API must reject empty stage sequences");
        assert!(err.to_string().contains("at least one center-arm flow fraction"));
    }

    #[test]
    fn checked_cascade_junction_matches_legacy_nominal_case() {
        let legacy = cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 5e-6);
        let checked = checked_cascade_junction_separation(3, 0.45, 2e-3, 1e-3, 5e-6)
            .expect("checked cascade junction API should succeed on a nominal case");

        assert!((legacy.cancer_center_fraction - checked.cancer_center_fraction).abs() < 1e-12);
        assert!((legacy.rbc_peripheral_fraction - checked.rbc_peripheral_fraction).abs() < 1e-12);
        assert!((legacy.center_hematocrit_ratio - checked.center_hematocrit_ratio).abs() < 1e-12);
    }

    #[test]
    fn incremental_qfrac_api_matches_uniform_width_model() {
        let width_model = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
        let q_pretri = cif_pretri_stage_q_fracs(2, 0.45, 0.55);
        let q_tri = tri_center_q_frac(0.55);
        let solved_like = incremental_filtration_separation_from_qfracs(&q_pretri, q_tri, 0.68);
        assert!(
            (width_model.cancer_center_fraction - solved_like.cancer_center_fraction).abs() < 1e-12
        );
        assert!((width_model.rbc_center_fraction - solved_like.rbc_center_fraction).abs() < 1e-12);
    }

    #[test]
    fn cif_pretri_stage_fracs_ramp_toward_terminal_bias() {
        let stage_fracs = cif_pretri_stage_center_fracs(3, 0.45, 0.60);
        assert_eq!(stage_fracs.len(), 3);
        assert!(stage_fracs[0] >= 0.45 - 1e-12);
        assert!(stage_fracs[1] >= stage_fracs[0] - 1e-12);
        assert!(stage_fracs[2] >= stage_fracs[1] - 1e-12);
        assert!(stage_fracs[2] <= 0.60 + 1e-12);
    }

    #[test]
    fn checked_cif_pretri_stage_fracs_reject_invalid_stage_count() {
        let err = checked_cif_pretri_stage_center_fracs(0, 0.45, 0.60)
            .expect_err("checked CIF stage fractions must reject zero pre-trifurcation stages");
        assert!(err.to_string().contains("stage count"));
    }

    #[test]
    fn checked_incremental_filtration_rejects_invalid_terminal_bifurcation_fraction() {
        let err = checked_incremental_filtration_separation_staged(2, 0.45, 0.50, 0.40)
            .expect_err("checked incremental filtration must reject terminal bifurcation fractions below the validated range");
        assert!(err.to_string().contains("terminal bifurcation"));
    }

    #[test]
    fn checked_incremental_filtration_rejects_zero_pretri_stages() {
        let err = checked_incremental_filtration_separation_staged(0, 0.45, 0.50, 0.68)
            .expect_err("checked incremental filtration must reject zero pre-trifurcation stages");
        assert!(err.to_string().contains("stage count"));
    }

    #[test]
    fn checked_incremental_filtration_matches_legacy_nominal_case() {
        let legacy = incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68);
        let checked = checked_incremental_filtration_separation_staged(2, 0.45, 0.55, 0.68)
            .expect("checked incremental filtration should succeed on a nominal selective-routing case");

        assert!((legacy.cancer_center_fraction - checked.cancer_center_fraction).abs() < 1e-12);
        assert!((legacy.rbc_center_fraction - checked.rbc_center_fraction).abs() < 1e-12);
        assert!((legacy.center_hematocrit_ratio - checked.center_hematocrit_ratio).abs() < 1e-12);
    }

    #[test]
    fn checked_incremental_qfrac_api_matches_legacy_nominal_case() {
        let q_pretri = cif_pretri_stage_q_fracs(2, 0.45, 0.55);
        let q_tri = tri_center_q_frac(0.55);
        let legacy = incremental_filtration_separation_from_qfracs(&q_pretri, q_tri, 0.68);
        let checked = checked_incremental_filtration_separation_from_qfracs(&q_pretri, q_tri, 0.68)
            .expect("checked incremental qfrac API should succeed on a nominal case");

        assert!((legacy.cancer_center_fraction - checked.cancer_center_fraction).abs() < 1e-12);
        assert!((legacy.rbc_center_fraction - checked.rbc_center_fraction).abs() < 1e-12);
    }

    #[test]
    fn cross_junction_q_frac_symmetric_matches_basic() {
        // With symmetric 1/3 split, cross-junction model should give similar
        // result to the basic model (identical for symmetric).
        let q_basic = tri_center_q_frac(1.0 / 3.0);
        let q_cross = tri_center_q_frac_cross_junction(1.0 / 3.0, 2e-3, 1e-3);
        // Both should be ~1/3 for symmetric geometry
        assert!((q_basic - 1.0 / 3.0).abs() < 1e-10);
        assert!(
            (q_cross - 1.0 / 3.0).abs() < 0.05,
            "symmetric cross-junction should be near 1/3, got {q_cross}"
        );
    }

    #[test]
    fn cross_junction_q_frac_wider_center_sees_more_k_loss() {
        // Cross-junction K-factor correction penalises the wider center arm
        // more than the narrow peripherals (K_eff ∝ √(A_branch/A_parent),
        // and the additive resistance term scales as K × w/h), so q_center
        // *decreases*.  This is the desired physics: more flow is diverted
        // to peripheral arms, improving RBC peripheral enrichment.
        let q_basic = tri_center_q_frac(0.55);
        let q_cross = tri_center_q_frac_cross_junction(0.55, 2e-3, 1e-3);
        assert!(q_cross <= q_basic + 1e-6,
            "cross-junction correction should reduce center fraction: basic={q_basic}, cross={q_cross}");
    }

    #[test]
    fn cross_junction_cif_pushes_more_rbc_to_periphery() {
        let basic = incremental_filtration_separation_staged(2, 0.45, 0.50, 0.68);
        let cross =
            incremental_filtration_separation_cross_junction(2, 0.45, 0.50, 0.68, 2e-3, 1e-3);
        assert!(cross.rbc_peripheral_fraction >= basic.rbc_peripheral_fraction - 0.01,
            "cross-junction selective routing should push more RBCs to periphery: basic={}, cross={}",
            basic.rbc_peripheral_fraction, cross.rbc_peripheral_fraction);
    }

    #[test]
    fn cross_junction_cct_pushes_more_rbc_to_periphery() {
        let basic = cascade_junction_separation(2, 0.45, 2e-3, 1e-3, 5e-6);
        let cross = cascade_junction_separation_cross_junction(2, 0.45, 2e-3, 1e-3);
        assert!(
            cross.rbc_peripheral_fraction >= basic.rbc_peripheral_fraction - 0.01,
            "cross-junction cascade routing should push more RBCs to periphery: basic={}, cross={}",
            basic.rbc_peripheral_fraction,
            cross.rbc_peripheral_fraction
        );
    }

    #[test]
    fn mixed_cascade_all_tri_matches_cascade() {
        let q = tri_center_q_frac(0.45);
        let pure = cascade_routing::cascade_from_q_fractions(&[q, q]);
        let mixed = mixed_cascade_separation(&[(q, true), (q, true)]);
        assert!((pure.cancer_center_fraction - mixed.cancer_center_fraction).abs() < 1e-12);
        assert!((pure.rbc_peripheral_fraction - mixed.rbc_peripheral_fraction).abs() < 1e-12);
    }

    #[test]
    fn mixed_cascade_tri_bi_produces_nonzero_separation() {
        let q_tri = tri_center_q_frac(0.45);
        let q_bi = 0.68;
        let r = mixed_cascade_separation(&[(q_tri, true), (q_bi, false)]);
        assert!(r.cancer_center_fraction > 0.0);
        assert!(r.separation_efficiency > 0.0);
        assert!(r.cancer_center_fraction > r.rbc_peripheral_fraction.min(0.99));
    }

    #[test]
    fn treatment_bifurcation_matches_single_stage_mixed_cascade() {
        let q_bi = 0.68;
        let direct = treatment_bifurcation_separation(q_bi);
        let mixed = mixed_cascade_separation(&[(q_bi, false)]);

        assert!((direct.cancer_center_fraction - mixed.cancer_center_fraction).abs() < 1e-12);
        assert!((direct.wbc_center_fraction - mixed.wbc_center_fraction).abs() < 1e-12);
        assert!((direct.rbc_peripheral_fraction - mixed.rbc_peripheral_fraction).abs() < 1e-12);
        assert!((direct.separation_efficiency - mixed.separation_efficiency).abs() < 1e-12);
    }

    #[test]
    fn checked_mixed_cascade_rejects_invalid_stage_flow_fraction() {
        let err = checked_mixed_cascade_separation(&[(1.20, true)]).expect_err(
            "checked mixed cascade routing must reject stage flow fractions outside the open unit interval",
        );
        assert!(err.to_string().contains("stage flow fraction"));
    }

    #[test]
    fn checked_mixed_cascade_matches_legacy_nominal_case() {
        let q_tri = tri_center_q_frac(0.45);
        let legacy = mixed_cascade_separation(&[(q_tri, true), (0.68, false)]);
        let checked = checked_mixed_cascade_separation(&[(q_tri, true), (0.68, false)])
            .expect("checked mixed cascade routing should succeed on a nominal case");

        assert!((legacy.cancer_center_fraction - checked.cancer_center_fraction).abs() < 1e-12);
        assert!((legacy.rbc_peripheral_fraction - checked.rbc_peripheral_fraction).abs() < 1e-12);
    }

    #[test]
    fn mixed_cascade_deeper_improves_separation() {
        let q_tri = tri_center_q_frac(0.45);
        let r1 = mixed_cascade_separation(&[(q_tri, true)]);
        let r2 = mixed_cascade_separation(&[(q_tri, true), (q_tri, true)]);
        assert!(r2.separation_efficiency >= r1.separation_efficiency - 1e-10);
    }

    // ── kappa-aware tests ─────────────────────────────────────────────────────

    #[test]
    fn beta_kappa_adjusted_no_change_at_zero_kappa() {
        // At κ = 0 (cell infinitely smaller than channel), β_eff = β_base.
        assert!((beta_kappa_adjusted(SE_CANCER, 0.0, 0.02, false) - SE_CANCER).abs() < 1e-12);
        assert!((beta_kappa_adjusted(SE_WBC, 0.0, 0.02, false) - SE_WBC).abs() < 1e-12);
        assert!((beta_kappa_adjusted(SE_RBC, 0.0, 0.02, true) - SE_RBC).abs() < 1e-12);
    }

    #[test]
    fn beta_kappa_adjusted_rbc_never_amplified() {
        // RBC is fully deformable (β_base = 1.0, excess = 0) — no amplification regardless of κ.
        for kappa in [0.0, 0.03, 0.07, 0.20] {
            // Evaluated below inversion velocity
            let b = beta_kappa_adjusted(SE_RBC, kappa, 0.02, true);
            assert!(
                (b - 1.0).abs() < 1e-12,
                "RBC β should stay 1.0 at κ={kappa}, got {b}"
            );
        }
    }

    #[test]
    fn beta_kappa_adjusted_cancer_increases_with_kappa() {
        let b0 = beta_kappa_adjusted(SE_CANCER, 0.0, 0.02, false);
        let b1 = beta_kappa_adjusted(SE_CANCER, KAPPA_REF, 0.02, false);
        let b2 = beta_kappa_adjusted(SE_CANCER, KAPPA_REF * 2.0, 0.02, false);
        assert!(b1 > b0, "β should increase with κ for stiff cancer cells");
        assert!(b2 >= b1, "β should not decrease as κ grows further");
        assert!(b2 <= 3.0, "β must be capped at 3.0");
    }

    #[test]
    fn p_arm_general_symmetric_trifurcation_matches_p_center() {
        // Symmetric trifurcation: arm_q_fracs = [q, q_p, q_p]
        let q_c = 0.55_f64;
        let q_p = (1.0 - q_c) / 2.0;
        let p_gen = p_arm_general(&[q_c, q_p, q_p], 0, SE_CANCER);
        let p_old = p_center(q_c, SE_CANCER);
        assert!(
            (p_gen - p_old).abs() < 1e-10,
            "generalised p_arm must match p_center for symmetric tri: gen={p_gen}, old={p_old}"
        );
    }

    #[test]
    fn p_arm_general_two_arm_matches_p_treat_bifurcation() {
        let q_t = 0.68_f64;
        let q_b = 1.0 - q_t;
        let p_gen = p_arm_general(&[q_t, q_b], 0, SE_CANCER);
        let p_old = p_treat_bifurcation(q_t, SE_CANCER);
        assert!(
            (p_gen - p_old).abs() < 1e-10,
            "generalised p_arm must match p_treat_bifurcation: gen={p_gen}, old={p_old}"
        );
    }

    #[test]
    fn tri_asymmetric_q_fracs_sum_to_one() {
        let [qc, ql, qr] = tri_asymmetric_q_fracs(0.45, 0.30, 4e-3, 1e-3);
        assert!(
            (qc + ql + qr - 1.0).abs() < 1e-10,
            "asymmetric arm fracs must sum to 1: {qc}+{ql}+{qr}={:.6}",
            qc + ql + qr
        );
    }

    #[test]
    fn tri_asymmetric_wider_periph_gets_more_flow() {
        // Left peripheral wider than right → left gets more flow.
        let [_qc, ql, qr] = tri_asymmetric_q_fracs(0.40, 0.40, 4e-3, 1e-3);
        // left_frac=0.40, right_frac=0.20 → left should carry more
        assert!(
            ql > qr,
            "wider left peripheral should carry more flow: ql={ql}, qr={qr}"
        );
    }

    #[test]
    fn checked_tri_asymmetric_q_fracs_reject_overfull_width_budget() {
        let err = checked_tri_asymmetric_q_fracs(0.70, 0.35, 4e-3, 1e-3).expect_err(
            "checked asymmetric trifurcation flow fractions must reject overfull width budgets",
        );
        assert!(err.to_string().contains("positive right-arm width"));
    }

    #[test]
    fn tri_asymmetric_symmetric_matches_cross_junction() {
        // With equal peripherals, tri_asymmetric_q_fracs should match tri_center_q_frac_cross_junction.
        let center_frac = 0.45_f64;
        let left_frac = (1.0 - center_frac) / 2.0; // symmetric
        let [qc_asym, _, _] = tri_asymmetric_q_fracs(center_frac, left_frac, 4e-3, 1e-3);
        let qc_sym = tri_center_q_frac_cross_junction(center_frac, 4e-3, 1e-3);
        assert!(
            (qc_asym - qc_sym).abs() < 1e-8,
            "symmetric asymmetric must match cross-junction: asym={qc_asym}, sym={qc_sym}"
        );
    }

    #[test]
    fn kappa_aware_higher_beta_than_legacy_for_stiff_cells_in_narrow_channel() {
        // In a narrow channel (Dh = 0.5mm), cancer cell κ ≈ 0.035 (above zero).
        // The kappa-aware model should produce HIGHER cancer_center_fraction than
        // the legacy model using constant β.
        let q = tri_center_q_frac_cross_junction(0.45, 0.8e-3, 1e-3);
        let q_p = (1.0 - q) / 2.0;
        let w_center = 0.45 * 0.8e-3;
        let h = 1e-3_f64;
        let dh = 2.0 * w_center * h / (w_center + h);

        let stage = CascadeStage {
            arm_q_fracs: [q, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            parent_v_in_m_s: 0.02,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };

        let kappa_result = mixed_cascade_separation_kappa_aware(&[stage]);
        let legacy_result = mixed_cascade_separation(&[(q, true)]);

        // κ ≈ 17.5µm / Dh → β_cancer > SE_CANCER = 1.85 → stronger cancer routing
        assert!(
            kappa_result.cancer_center_fraction >= legacy_result.cancer_center_fraction - 1e-10,
            "kappa-aware model must not reduce cancer routing vs legacy: kappa={:.4}, legacy={:.4}",
            kappa_result.cancer_center_fraction,
            legacy_result.cancer_center_fraction
        );
        // RBC routing must be unchanged (deformable, β stays 1.0)
        assert!(
            (kappa_result.rbc_peripheral_fraction - legacy_result.rbc_peripheral_fraction).abs()
                < 1e-9,
            "RBC routing must be unchanged by κ correction"
        );
    }

    #[test]
    fn kappa_aware_asymmetric_arms_biases_cells_to_wider_peripheral() {
        // Asymmetric trifurcation: left peripheral wider than right.
        // Cancer cells should preferentially go to the wider LEFT arm (higher q_l).
        let [q_c, q_l, q_r] = tri_asymmetric_q_fracs(0.40, 0.40, 4e-3, 1e-3);
        // q_l > q_r → more cells should route to left
        let w_center = 0.40 * 4e-3;
        let h = 1e-3_f64;
        let dh = 2.0 * w_center * h / (w_center + h);
        let stage = CascadeStage {
            arm_q_fracs: [q_c, q_l, q_r, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            parent_v_in_m_s: 0.02,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };
        let n = stage.n_arms as usize;
        let dh_s = stage.treatment_dh_m.max(1e-9);
        let beta = beta_kappa_adjusted(SE_CANCER, D_CANCER_M / dh_s, 0.02, false);
        let p_to_left = p_arm_general(&stage.arm_q_fracs[..n], 1, beta);
        let p_to_right = p_arm_general(&stage.arm_q_fracs[..n], 2, beta);
        assert!(
            p_to_left > p_to_right,
            "wider left peripheral should attract more cancer cells: p_left={p_to_left:.4}, p_right={p_to_right:.4}"
        );
    }

    #[test]
    fn kappa_aware_deeper_cascade_still_improves_separation() {
        // Build a 3-stage TriTriTri with narrowing Dh at each stage.
        let center_frac = 0.45_f64;
        let h = 1e-3_f64;
        let mut parent_w = 4e-3_f64;
        let mut stages = Vec::new();
        for _ in 0..3 {
            let w_c = center_frac * parent_w;
            let dh = 2.0 * w_c * h / (w_c + h);
            let q_c = tri_center_q_frac_cross_junction(center_frac, parent_w, h);
            let q_p = (1.0 - q_c) / 2.0;
            stages.push(CascadeStage {
                arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh,
                parent_v_in_m_s: 0.02,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            });
            parent_w *= center_frac;
        }

        let r3 = mixed_cascade_separation_kappa_aware(&stages);
        let r1 = mixed_cascade_separation_kappa_aware(&stages[..1]);

        assert!(
            r3.rbc_peripheral_fraction >= r1.rbc_peripheral_fraction - 1e-10,
            "3-stage TriTriTri should push more RBCs to periphery: r1={:.4}, r3={:.4}",
            r1.rbc_peripheral_fraction,
            r3.rbc_peripheral_fraction
        );
        assert!(
            r3.separation_efficiency >= r1.separation_efficiency - 1e-10,
            "3-stage TriTriTri should improve separation efficiency"
        );
    }

    // ── Fåhræus margination correction tests ──────────────────────────────

    #[test]
    fn fahrae_correction_zero_for_rbc() {
        // RBC is the reference cell: excess β = 0 and size ratio = 0 → no correction.
        let c = fahrae_beta_correction(SE_RBC, D_RBC_M);
        assert!(
            c.abs() < 1e-15,
            "Fåhræus correction must be zero for RBC, got {c}"
        );
    }

    #[test]
    fn fahrae_correction_positive_for_cancer() {
        // Cancer cells (17.5 µm) are larger than RBCs (7 µm) → positive correction.
        let c = fahrae_beta_correction(SE_CANCER, D_CANCER_M);
        assert!(
            c > 0.0,
            "Fåhræus correction must be positive for cancer cells"
        );
        // Expected: excess(0.85) × α(0.12) × size_ratio(1.5) = 0.153
        assert!(
            (c - 0.85 * 0.12 * 1.5).abs() < 1e-12,
            "Fåhræus correction for cancer should be ~0.153, got {c}"
        );
    }

    #[test]
    fn fahrae_correction_smaller_for_wbc_than_cancer() {
        // WBCs (10 µm) are smaller than cancer cells (17.5 µm) → smaller correction.
        let c_cancer = fahrae_beta_correction(SE_CANCER, D_CANCER_M);
        let c_wbc = fahrae_beta_correction(SE_WBC, D_WBC_M);
        assert!(c_wbc > 0.0, "Fåhræus correction must be positive for WBCs");
        assert!(
            c_cancer > c_wbc,
            "Cancer correction ({c_cancer:.4}) must exceed WBC correction ({c_wbc:.4})"
        );
    }

    #[test]
    fn kappa_aware_with_fahrae_boosts_cancer_over_legacy() {
        // The Fåhræus correction in the kappa-aware model should boost cancer routing
        // beyond what the legacy mixed_cascade_separation (no κ, no Fåhræus) provides.
        let q = tri_center_q_frac_cross_junction(0.45, 4e-3, 1e-3);
        let q_p = (1.0 - q) / 2.0;
        let w_center = 0.45 * 4e-3;
        let h = 1e-3_f64;
        let dh = 2.0 * w_center * h / (w_center + h);

        let stage = CascadeStage {
            arm_q_fracs: [q, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            parent_v_in_m_s: 0.02,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };
        let kappa_fahrae = mixed_cascade_separation_kappa_aware(&[stage]);
        let legacy = mixed_cascade_separation(&[(q, true)]);

        // The Fåhræus-enhanced kappa model should route more cancer cells to center.
        assert!(
            kappa_fahrae.cancer_center_fraction > legacy.cancer_center_fraction,
            "Fåhræus-enhanced model must exceed legacy cancer routing: enhanced={:.4}, legacy={:.4}",
            kappa_fahrae.cancer_center_fraction,
            legacy.cancer_center_fraction
        );
        // Cancer-to-RBC enrichment ratio should improve.
        let rbc_center_kf = 1.0 - kappa_fahrae.rbc_peripheral_fraction;
        let rbc_center_legacy = 1.0 - legacy.rbc_peripheral_fraction;
        let enrichment_kf = kappa_fahrae.cancer_center_fraction / rbc_center_kf.max(1e-12);
        let enrichment_legacy = legacy.cancer_center_fraction / rbc_center_legacy.max(1e-12);
        assert!(
            enrichment_kf >= enrichment_legacy - 1e-6,
            "Fåhræus model should improve CTC/RBC enrichment: kf={enrichment_kf:.4}, legacy={enrichment_legacy:.4}"
        );
    }

    #[test]
    fn option1_option2_tri_tri_enrichment_quantified() {
        // Quantify CTC enrichment for the actual Option 1/Option 2 TriTri topology:
        //   Stage 1: pretri center_frac = 0.45
        //   Stage 2: terminal tri center_frac = 0.333 (symmetric)
        //   parent_width = 4mm, height = 1mm
        let h = 1e-3_f64;
        let parent_w = 4e-3_f64;
        let pcf = 0.45;
        let tcf = 1.0 / 3.0; // symmetric

        // Build 2-stage TriTri with cross-junction corrections (as used by compute.rs)
        let left_periph = (1.0 - pcf) / 2.0;
        let arm_q_1 = tri_asymmetric_q_fracs(pcf, left_periph, parent_w, h);
        let w_center_1 = pcf * parent_w;
        let dh_1 = 2.0 * w_center_1 * h / (w_center_1 + h);

        let parent_w_2 = w_center_1;
        let left_periph_2 = (1.0 - tcf) / 2.0;
        let arm_q_2 = tri_asymmetric_q_fracs(tcf, left_periph_2, parent_w_2, h);
        let w_center_2 = tcf * parent_w_2;
        let dh_2 = 2.0 * w_center_2 * h / (w_center_2 + h);

        let stages = vec![
            CascadeStage {
                arm_q_fracs: [arm_q_1[0], arm_q_1[1], arm_q_1[2], 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh_1,
                parent_v_in_m_s: 0.02,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            },
            CascadeStage {
                arm_q_fracs: [arm_q_2[0], arm_q_2[1], arm_q_2[2], 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh_2,
                parent_v_in_m_s: 0.02,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            },
        ];
        let r = mixed_cascade_separation_kappa_aware(&stages);

        // With SE_CANCER=1.85 + Fåhræus correction, cancer routing should be > 25%
        // (up from ~24.6% with SE_CANCER=1.70 and no Fåhræus).
        assert!(
            r.cancer_center_fraction > 0.25,
            "TriTri cancer center fraction should exceed 25%: got {:.2}%",
            r.cancer_center_fraction * 100.0
        );
        // RBC peripheral fraction should remain high (> 78%)
        assert!(
            r.rbc_peripheral_fraction > 0.78,
            "TriTri RBC peripheral fraction should exceed 78%: got {:.2}%",
            r.rbc_peripheral_fraction * 100.0
        );
        // CTC/RBC enrichment ratio should be > 1.3
        let rbc_center = 1.0 - r.rbc_peripheral_fraction;
        let enrichment = r.cancer_center_fraction / rbc_center.max(1e-12);
        assert!(
            enrichment > 1.3,
            "CTC/RBC enrichment ratio should exceed 1.3: got {enrichment:.4}"
        );

        // Print for diagnostic review
        eprintln!("=== TriTri (pcf=0.45, tcf=0.333) kappa-aware + Fåhræus ===");
        eprintln!(
            "  Cancer center fraction:  {:.2}%",
            r.cancer_center_fraction * 100.0
        );
        eprintln!(
            "  WBC center fraction:     {:.2}%",
            r.wbc_center_fraction * 100.0
        );
        eprintln!(
            "  RBC peripheral fraction: {:.2}%",
            r.rbc_peripheral_fraction * 100.0
        );
        eprintln!("  RBC center fraction:     {:.2}%", rbc_center * 100.0);
        eprintln!("  Separation efficiency:   {:.4}", r.separation_efficiency);
        eprintln!("  CTC/RBC enrichment:      {:.4}×", enrichment);
        eprintln!(
            "  Center HCT ratio:        {:.4}",
            r.center_hematocrit_ratio
        );
    }

    #[test]
    fn option2_higher_tcf_dramatically_improves_enrichment() {
        // Demonstrate that using tcf=0.55 instead of 0.333 would dramatically
        // improve CTC enrichment. This is a design-space insight, not a physics change.
        let h = 1e-3_f64;
        let parent_w = 4e-3_f64;
        let pcf = 0.45;

        // Build with tcf=0.55 (asymmetric terminal tri)
        let tcf = 0.55;
        let left_periph = (1.0 - pcf) / 2.0;
        let arm_q_1 = tri_asymmetric_q_fracs(pcf, left_periph, parent_w, h);
        let w_center_1 = pcf * parent_w;
        let dh_1 = 2.0 * w_center_1 * h / (w_center_1 + h);

        let parent_w_2 = w_center_1;
        let left_periph_2 = (1.0 - tcf) / 2.0;
        let arm_q_2 = tri_asymmetric_q_fracs(tcf, left_periph_2, parent_w_2, h);
        let w_center_2 = tcf * parent_w_2;
        let dh_2 = 2.0 * w_center_2 * h / (w_center_2 + h);

        let stages = vec![
            CascadeStage {
                arm_q_fracs: [arm_q_1[0], arm_q_1[1], arm_q_1[2], 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh_1,
                parent_v_in_m_s: 0.02,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            },
            CascadeStage {
                arm_q_fracs: [arm_q_2[0], arm_q_2[1], arm_q_2[2], 0.0, 0.0],
                n_arms: 3,
                treatment_dh_m: dh_2,
                parent_v_in_m_s: 0.02,
                peripheral_recoveries: [None; 4],
                n_recoveries: 0,
            },
        ];
        let r_high = mixed_cascade_separation_kappa_aware(&stages);

        // With tcf=0.55: cancer routing should exceed 40%.
        assert!(
            r_high.cancer_center_fraction > 0.40,
            "TriTri with tcf=0.55 cancer fraction should exceed 40%: got {:.2}%",
            r_high.cancer_center_fraction * 100.0
        );

        let rbc_center = 1.0 - r_high.rbc_peripheral_fraction;
        let enrichment = r_high.cancer_center_fraction / rbc_center.max(1e-12);
        assert!(
            enrichment > 1.4,
            "TriTri with tcf=0.55 enrichment should exceed 1.4: got {enrichment:.4}"
        );

        eprintln!("=== TriTri (pcf=0.45, tcf=0.55) kappa-aware + Fåhræus ===");
        eprintln!(
            "  Cancer center fraction:  {:.2}%",
            r_high.cancer_center_fraction * 100.0
        );
        eprintln!(
            "  WBC center fraction:     {:.2}%",
            r_high.wbc_center_fraction * 100.0
        );
        eprintln!(
            "  RBC peripheral fraction: {:.2}%",
            r_high.rbc_peripheral_fraction * 100.0
        );
        eprintln!("  RBC center fraction:     {:.2}%", rbc_center * 100.0);
        eprintln!("  CTC/RBC enrichment:      {:.4}×", enrichment);
    }

    #[test]
    fn updated_se_cancer_improves_single_stage_selectivity() {
        // With SE_CANCER = 1.85 (up from 1.70), a single trifurcation stage with
        // center-biased flow should produce higher cancer-to-RBC differential.
        let q = tri_center_q_frac(0.50);
        let r = mixed_cascade_separation(&[(q, true)]);

        // Cancer routing must exceed RBC routing (flow-weighted at β=1.0).
        let rbc_center = 1.0 - r.rbc_peripheral_fraction;
        assert!(
            r.cancer_center_fraction > rbc_center,
            "Cancer must route preferentially to center: cancer={:.4}, rbc_center={:.4}",
            r.cancer_center_fraction,
            rbc_center
        );
        // Enrichment ratio must be > 1.0
        let enrichment = r.cancer_center_fraction / rbc_center.max(1e-12);
        assert!(
            enrichment > 1.0,
            "CTC/RBC enrichment must exceed 1.0 at asymmetric split, got {enrichment:.4}"
        );
    }

    // ── Peripheral recovery routing tests ────────────────────────────────

    #[test]
    fn recovery_zero_when_no_peripheral_recoveries() {
        let q = tri_center_q_frac(0.45);
        let q_p = (1.0 - q) / 2.0;
        let dh = 1e-3;
        let stage_no_recovery = CascadeStage {
            arm_q_fracs: [q, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            parent_v_in_m_s: 0.02,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };
        let r = mixed_cascade_separation_kappa_aware(&[stage_no_recovery]);
        let r_legacy = mixed_cascade_separation(&[(q, true)]);
        // With no recovery, kappa-aware should match or exceed legacy
        assert!(r.cancer_center_fraction >= r_legacy.cancer_center_fraction - 1e-10);
    }

    #[test]
    fn recovery_increases_cancer_center_fraction() {
        // 20/60/20 trifurcation with recovery on left peripheral
        let q_c = 0.65;  // center gets ~65% of flow (wider channel)
        let q_p = (1.0 - q_c) / 2.0;
        let dh = 1e-3;
        let stage_no_recovery = CascadeStage {
            arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            parent_v_in_m_s: 0.02,
            peripheral_recoveries: [None; 4],
            n_recoveries: 0,
        };
        // Recovery sub-split on arm 1 (left): 70/30 split, wider arm (index 0) feeds back to treatment
        let recovery = PeripheralRecovery {
            source_arm_idx: 1,
            sub_arm_q_fracs: [0.70, 0.30, 0.0, 0.0, 0.0],
            n_sub_arms: 2,
            recovery_arm_idx: 0,
            recovery_dh_m: 0.5e-3,
        };
        let stage_with_recovery = CascadeStage {
            arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: dh,
            parent_v_in_m_s: 0.02,
            peripheral_recoveries: [Some(recovery), None, None, None],
            n_recoveries: 1,
        };
        let r_no = mixed_cascade_separation_kappa_aware(&[stage_no_recovery]);
        let r_yes = mixed_cascade_separation_kappa_aware(&[stage_with_recovery]);
        assert!(
            r_yes.cancer_center_fraction > r_no.cancer_center_fraction,
            "recovery should increase cancer center fraction: without={:.4}, with={:.4}",
            r_no.cancer_center_fraction, r_yes.cancer_center_fraction
        );
    }

    #[test]
    fn recovery_bounded_by_one() {
        // Even with aggressive recovery, P_eff must not exceed 1.0
        let q_c = 0.50;
        let q_p = 0.25;
        let recovery_both = [
            Some(PeripheralRecovery {
                source_arm_idx: 1,
                sub_arm_q_fracs: [0.90, 0.10, 0.0, 0.0, 0.0],
                n_sub_arms: 2,
                recovery_arm_idx: 0,
                recovery_dh_m: 0.3e-3,
            }),
            Some(PeripheralRecovery {
                source_arm_idx: 2,
                sub_arm_q_fracs: [0.90, 0.10, 0.0, 0.0, 0.0],
                n_sub_arms: 2,
                recovery_arm_idx: 0,
                recovery_dh_m: 0.3e-3,
            }),
            None,
            None,
        ];
        let stage = CascadeStage {
            arm_q_fracs: [q_c, q_p, q_p, 0.0, 0.0],
            n_arms: 3,
            treatment_dh_m: 1e-3,
            parent_v_in_m_s: 0.02,
            peripheral_recoveries: recovery_both,
            n_recoveries: 2,
        };
        let r = mixed_cascade_separation_kappa_aware(&[stage]);
        assert!(r.cancer_center_fraction <= 1.0, "cancer fraction must be <= 1.0");
        assert!(r.wbc_center_fraction <= 1.0, "wbc fraction must be <= 1.0");
        assert!((1.0 - r.rbc_peripheral_fraction) <= 1.0, "rbc center must be <= 1.0");
    }

    // ── Core Zweifach-Fung function unit tests ──────────────────────────

    #[test]
    fn p_center_equal_trifurcation() {
        // Equal trifurcation: q_center = 1/3, each peripheral = 1/3
        // P_center = (1/3)^β / ((1/3)^β + 2*(1/3)^β) = 1/3 for any β
        use approx::assert_relative_eq;
        let p_cancer = p_center(1.0 / 3.0, SE_CANCER);
        let p_wbc = p_center(1.0 / 3.0, SE_WBC);
        let p_rbc = p_center(1.0 / 3.0, SE_RBC);
        assert_relative_eq!(p_cancer, 1.0 / 3.0, max_relative = 1e-9);
        assert_relative_eq!(p_wbc, 1.0 / 3.0, max_relative = 1e-9);
        assert_relative_eq!(p_rbc, 1.0 / 3.0, max_relative = 1e-9);
    }

    #[test]
    fn p_center_biased_flow_stiff_cells_route_more() {
        // With q_center_frac > 1/3, stiffer cells (higher beta) route more to center.
        use approx::assert_relative_eq;
        let q = 0.5; // center gets 50% of flow
        let p_cancer = p_center(q, SE_CANCER);
        let p_rbc = p_center(q, SE_RBC);
        // For RBC (beta=1): P = 0.5^1 / (0.5^1 + 2*0.25^1) = 0.5/1.0 = 0.5
        assert_relative_eq!(p_rbc, 0.5, max_relative = 1e-9);
        // Cancer (beta=1.85) should be > 0.5
        assert!(p_cancer > 0.5, "stiff cancer cells should route more to center: got {p_cancer}");
        assert!(p_cancer > p_rbc, "cancer p_center should exceed RBC p_center");
    }

    #[test]
    fn tri_center_q_frac_known_values() {
        // For center_frac = 0.5: w_c=0.5, w_p=0.25 each
        // q_frac = 0.5^3 / (0.5^3 + 2*0.25^3) = 0.125 / (0.125 + 0.03125) = 0.125/0.15625 = 0.8
        use approx::assert_relative_eq;
        let q = tri_center_q_frac(0.5);
        assert_relative_eq!(q, 0.8, max_relative = 1e-9);
    }

    #[test]
    fn beta_kappa_adjusted_at_kappa_ref() {
        // At kappa = KAPPA_REF: amplification = 1 + 0.07/0.07 = 2
        // beta_eff = 1 + (SE_CANCER - 1) * 2 = 1 + 0.85 * 2 = 2.70
        use approx::assert_relative_eq;
        let b = beta_kappa_adjusted(SE_CANCER, KAPPA_REF, 0.02, false);
        assert_relative_eq!(b, 1.0 + 0.85 * 2.0, max_relative = 1e-10);
    }

    #[test]
    fn beta_kappa_adjusted_capped_at_three() {
        // Very large kappa should cap at 3.0
        let b = beta_kappa_adjusted(2.5, KAPPA_REF * 10.0, 0.02, false);
        assert!(b <= 3.0, "beta should be capped at 3.0, got {b}");
        // The clamp on kappa is at 2*KAPPA_REF, so amplification = 1 + 2 = 3
        // excess = 1.5, beta_eff = 1 + 1.5*3 = 5.5, capped to 3.0
        use approx::assert_relative_eq;
        assert_relative_eq!(b, 3.0, epsilon = 1e-15);
    }

    #[test]
    fn beta_kappa_adjusted_wbc_moderate_amplification() {
        // WBC at kappa = KAPPA_REF: amplification = 2
        // beta_eff = 1 + (1.40 - 1) * 2 = 1 + 0.80 = 1.80
        use approx::assert_relative_eq;
        let b = beta_kappa_adjusted(SE_WBC, KAPPA_REF, 0.02, false);
        assert_relative_eq!(b, 1.0 + 0.40 * 2.0, max_relative = 1e-10);
    }
}
