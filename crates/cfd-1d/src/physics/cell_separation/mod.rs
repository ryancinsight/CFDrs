//! Inertial microfluidic cell separation for millifluidic device design.
//!
//! This module provides physics-based models for predicting the lateral
//! equilibrium positions of cell populations in rectangular microchannels,
//! enabling design of devices that separate cancer cells from healthy blood
//! cells using inertial focusing and Dean flow.
//!
//! # Physical basis
//!
//! In a straight rectangular channel at finite Reynolds number, particles
//! experience inertial lift forces that drive them to stable equilibrium
//! positions.  Larger, stiffer particles (cancer cells) focus near the
//! channel center; smaller, more deformable particles (RBCs) focus near
//! the walls.  In curved channels, Dean flow secondary circulation enhances
//! this separation.
//!
//! # Module structure
//!
//! | Module | Contents |
//! |--------|----------|
//! | [`properties`] | Cell physical properties (size, deformability, density) |
//! | [`margination`] | Inertial lift and Dean drag force models |
//! | [`separation_model`] | High-level separation analysis and efficiency metrics |
//! | [`cell_interaction`] | RBC-WBC cell-cell interaction correction (CFL model, HCT ≥ 1%) |
//!
//! # Quick start
//!
//! ```rust,ignore
//! use cfd_1d::cell_separation::{CellProperties, CellSeparationModel};
//!
//! let cancer = CellProperties::mcf7_breast_cancer();
//! let healthy = CellProperties::red_blood_cell();
//!
//! let model = CellSeparationModel::new(500e-6, 200e-6, 0.01, None);
//! let analysis = model.analyze(&cancer, &healthy, 1060.0, 3.5e-3, 0.05)
//!     .expect("cells must focus (κ > 0.07)");
//!
//! println!("Separation efficiency: {:.2}", analysis.separation_efficiency);
//! println!("Cancer center fraction: {:.2}", analysis.target_center_fraction);
//! ```
//!
//! # References
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.
//! - Gossett, D. R. & Di Carlo, D. (2009). *Anal. Chem.*, 81, 8459–8465.
//! - Hur, S. C. et al. (2011). *Lab Chip*, 11, 912–920.
//! - Fedosov, D. A. et al. (2012). Margination of white blood cells in microcapillary flow.
//!   *Phys. Rev. Lett.*, 108, 028104.

pub mod cascade_junction;
pub mod cell_free_layer;
pub mod cell_interaction;
pub mod fahraeus_effect;
pub mod fahraeus_lindqvist;
pub mod margination;
pub mod plasma_skimming;
pub mod properties;
pub mod rouleaux_aggregation;
pub mod separation_model;
pub mod zweifach_fung;

pub use cascade_junction::{
    cascade_junction_separation, cascade_junction_separation_cross_junction,
    cascade_junction_separation_from_qfracs, cif_pretri_stage_center_fracs,
    cif_pretri_stage_q_fracs, cif_pretri_stage_q_fracs_cross_junction,
    incremental_filtration_separation_cross_junction,
    incremental_filtration_separation_from_qfracs, incremental_filtration_separation_staged,
    mixed_cascade_separation, mixed_cascade_separation_kappa_aware, tri_asymmetric_q_fracs,
    treatment_bifurcation_separation, tri_center_q_frac, tri_center_q_frac_cross_junction,
    CascadeJunctionResult, CascadeStage, IncrementalFiltrationResult, PeripheralRecovery,
};
pub use cell_interaction::enhanced_lateral_equilibrium;
pub use margination::{
    amini_confinement_correction, checked_amini_confinement_correction,
    checked_inertial_lift_force_n, checked_lateral_equilibrium, checked_lateral_velocity_m_s,
    dean_drag_force_n, dean_number, inertial_lift_force_n, lateral_equilibrium,
    EquilibriumResult, AMINI_ALPHA_CONFINEMENT, AMINI_KAPPA_REF,
};
pub use fahraeus_lindqvist::{
    fahraeus_lindqvist_viscosity, secomb_network_viscosity, secomb_phase_separation_x0,
};
pub use cell_free_layer::{cfl_width_fedosov, cfl_width_sharan_popel, two_layer_viscosity};
pub use fahraeus_effect::{discharge_hematocrit, tube_hematocrit, tube_hematocrit_ratio};
pub use plasma_skimming::{
    checked_plasma_skimming_hematocrit, checked_pries_phase_separation,
    plasma_skimming_hematocrit, pries_phase_separation, PhaseSeparationResult,
};
pub use properties::CellProperties;
pub use rouleaux_aggregation::{checked_quemada_viscosity, quemada_viscosity};
pub use separation_model::{CellSeparationAnalysis, CellSeparationModel};
pub use cell_interaction::checked_enhanced_lateral_equilibrium;
pub use zweifach_fung::{
    checked_confinement_ratio, checked_critical_fractional_flow,
    checked_zweifach_fung_daughter_hematocrits, checked_zweifach_fung_rbc_fraction,
    confinement_ratio, critical_fractional_flow, zweifach_fung_daughter_hematocrits,
    zweifach_fung_rbc_fraction,
};
use cfd_core::error::{Error, Result};

// ── Three-population simultaneous model ─────────────────────────────────────

/// Simultaneously-computed lateral equilibrium positions for three cell populations.
///
/// All positions are derived from identical flow conditions: same channel geometry,
/// flow rate, fluid properties, and curvature — ensuring self-consistent metrics.
#[derive(Debug, Clone)]
pub struct ThreePopEquilibria {
    /// Normalised equilibrium position of MCF-7 cancer cells [0 = center, 1 = wall].
    pub cancer_eq: f64,
    /// Normalised equilibrium position of WBCs, CFL-corrected at `hematocrit` [0 = center, 1 = wall].
    pub wbc_eq: f64,
    /// Normalised equilibrium position of RBCs [0 = center, 1 = wall].
    pub rbc_eq: f64,
    /// Absolute cancer/WBC separation: `|x̃_cancer − x̃_wbc|`.
    pub cancer_wbc_sep: f64,
    /// RBC vs centre separation: `x̃_rbc − max(x̃_cancer, x̃_wbc)`, clamped to [−1, 1].
    pub rbc_vs_center_sep: f64,
    /// Fraction of cancer cells in the centre channel (`x̃ < 0.3`).
    pub cancer_center_fraction: f64,
    /// Fraction of WBCs in the centre channel (`x̃ < 0.3`).
    pub wbc_center_fraction: f64,
    /// Fraction of RBCs in the peripheral channels (`x̃ > 0.3`).
    pub rbc_peripheral_fraction: f64,
}

/// Compute three-population inertial focusing equilibria simultaneously under
/// identical flow conditions.
///
/// Cancer cells and RBCs use the plain inertial lift model
/// ([`margination::lateral_equilibrium`]).  WBCs additionally receive the
/// empirical CFL correction ([`cell_interaction::enhanced_lateral_equilibrium`])
/// because nucleated cells marginate toward the RBC-depleted layer at
/// physiological haematocrit (Fedosov 2012).
///
/// All three populations are evaluated at the **same** mean velocity, viscosity,
/// and Dean number — avoiding the inconsistency of two separate pairwise calls
/// with potentially different flow conditions.
///
/// # Arguments
/// - `width_m` — channel width (wider dimension) [m]
/// - `height_m` — channel height (shorter dimension) [m]
/// - `flow_rate_m3_s` — volumetric flow rate [m³/s]
/// - `blood_density` — fluid density [kg/m³]
/// - `viscosity` — dynamic viscosity [Pa·s]
/// - `hematocrit` — volumetric RBC fraction (0.0–0.45); used for WBC CFL correction
/// - `bend_radius_m` — radius of curvature [m], or `None` for straight channel
///
/// # Equilibrium convention
/// `x̃ = 0` → channel centre; `x̃ = 1` → channel wall.
/// Centre fraction: cells with `x̃ < 0.3` are in the centre stream.
/// Peripheral fraction: cells with `x̃ > 0.3` are in the peripheral streams.
#[must_use]
pub fn three_population_equilibria(
    width_m: f64,
    height_m: f64,
    flow_rate_m3_s: f64,
    blood_density: f64,
    viscosity: f64,
    hematocrit: f64,
    bend_radius_m: Option<f64>,
) -> ThreePopEquilibria {
    checked_three_population_equilibria(
        width_m,
        height_m,
        flow_rate_m3_s,
        blood_density,
        viscosity,
        hematocrit,
        bend_radius_m,
    )
    .unwrap_or(ThreePopEquilibria {
        cancer_eq: 0.5,
        wbc_eq: 0.5,
        rbc_eq: 0.5,
        cancer_wbc_sep: 0.0,
        rbc_vs_center_sep: 0.0,
        cancer_center_fraction: 0.0,
        wbc_center_fraction: 0.0,
        rbc_peripheral_fraction: 0.0,
    })
}

/// Checked three-population focusing analysis that rejects invalid operating conditions.
pub fn checked_three_population_equilibria(
    width_m: f64,
    height_m: f64,
    flow_rate_m3_s: f64,
    blood_density: f64,
    viscosity: f64,
    hematocrit: f64,
    bend_radius_m: Option<f64>,
) -> Result<ThreePopEquilibria> {
    let cancer = CellProperties::mcf7_breast_cancer();
    let wbc = CellProperties::white_blood_cell();
    let rbc = CellProperties::red_blood_cell();

    for (name, value) in [
        ("width", width_m),
        ("height", height_m),
        ("flow rate", flow_rate_m3_s),
        ("density", blood_density),
        ("viscosity", viscosity),
        ("hematocrit", hematocrit),
    ] {
        if !value.is_finite() {
            return Err(Error::InvalidConfiguration(format!(
                "Three-population cell separation {name} must be finite"
            )));
        }
    }
    if width_m <= 0.0 || height_m <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Three-population cell separation geometry must be positive".to_string(),
        ));
    }
    if flow_rate_m3_s <= 0.0 || blood_density <= 0.0 || viscosity <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Three-population cell separation flow and fluid properties must be positive".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&hematocrit) {
        return Err(Error::InvalidConfiguration(
            "Three-population cell separation hematocrit must lie in [0, 1]".to_string(),
        ));
    }
    if let Some(radius) = bend_radius_m {
        if !radius.is_finite() || radius <= 0.0 {
            return Err(Error::InvalidConfiguration(
                "Three-population cell separation bend radius must be finite and positive".to_string(),
            ));
        }
    }

    let area = width_m * height_m;
    let mean_v = flow_rate_m3_s / area;

    const SPLIT: f64 = 0.3;

    let cancer_x = margination::checked_lateral_equilibrium(
        &cancer,
        blood_density,
        viscosity,
        mean_v,
        width_m,
        height_m,
        bend_radius_m,
    )?
    .x_tilde_eq;

    let wbc_x = cell_interaction::checked_enhanced_lateral_equilibrium(
        &wbc,
        blood_density,
        viscosity,
        mean_v,
        width_m,
        height_m,
        bend_radius_m,
        hematocrit,
    )?
    .x_tilde_eq;

    let rbc_x = margination::checked_lateral_equilibrium(
        &rbc,
        blood_density,
        viscosity,
        mean_v,
        width_m,
        height_m,
        bend_radius_m,
    )?
    .x_tilde_eq;

    let cancer_wbc_sep = (cancer_x - wbc_x).abs();
    let max_central_x = cancer_x.max(wbc_x);
    let rbc_vs_center_sep = (rbc_x - max_central_x).clamp(-1.0, 1.0);

    let cancer_center_fraction = if cancer_x < SPLIT {
        (1.0 - cancer_x / SPLIT).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let wbc_center_fraction = if wbc_x < SPLIT {
        (1.0 - wbc_x / SPLIT).clamp(0.0, 1.0)
    } else {
        0.0
    };
    let rbc_peripheral_fraction = if rbc_x > SPLIT {
        ((rbc_x - SPLIT) / (1.0 - SPLIT)).clamp(0.0, 1.0)
    } else {
        0.0
    };

    Ok(ThreePopEquilibria {
        cancer_eq: cancer_x,
        wbc_eq: wbc_x,
        rbc_eq: rbc_x,
        cancer_wbc_sep,
        rbc_vs_center_sep,
        cancer_center_fraction,
        wbc_center_fraction,
        rbc_peripheral_fraction,
    })
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Typical millifluidic conditions for blood flow.
    /// Width = 500 µm, height = 200 µm, Q = 50 µL/min = 8.33e-10 m³/s.
    const WIDTH: f64 = 500e-6;
    const HEIGHT: f64 = 200e-6;
    const Q: f64 = 8.33e-10;
    const RHO: f64 = 1060.0;
    const MU: f64 = 3.5e-3;
    const HCT: f64 = 0.30;

    /// Ordering invariant: cancer cells focus more centripetally than RBCs.
    /// cancer_eq ≤ rbc_eq (lower x̃ = closer to center).
    #[test]
    fn ordering_cancer_center_rbc_wall() {
        let eq = three_population_equilibria(WIDTH, HEIGHT, Q, RHO, MU, HCT, None);
        assert!(
            eq.cancer_eq <= eq.rbc_eq,
            "Cancer (x̃={:.3}) should focus closer to center than RBC (x̃={:.3})",
            eq.cancer_eq,
            eq.rbc_eq
        );
    }

    /// All fractions must be in [0, 1].
    #[test]
    fn fractions_bounded() {
        let eq = three_population_equilibria(WIDTH, HEIGHT, Q, RHO, MU, HCT, None);
        assert!(
            eq.cancer_center_fraction >= 0.0 && eq.cancer_center_fraction <= 1.0,
            "cancer_center_fraction = {} out of [0,1]",
            eq.cancer_center_fraction
        );
        assert!(
            eq.wbc_center_fraction >= 0.0 && eq.wbc_center_fraction <= 1.0,
            "wbc_center_fraction = {} out of [0,1]",
            eq.wbc_center_fraction
        );
        assert!(
            eq.rbc_peripheral_fraction >= 0.0 && eq.rbc_peripheral_fraction <= 1.0,
            "rbc_peripheral_fraction = {} out of [0,1]",
            eq.rbc_peripheral_fraction
        );
    }

    /// Cancer-WBC separation is always non-negative.
    #[test]
    fn cancer_wbc_separation_nonneg() {
        let eq = three_population_equilibria(WIDTH, HEIGHT, Q, RHO, MU, HCT, None);
        assert!(
            eq.cancer_wbc_sep >= 0.0,
            "cancer_wbc_sep = {} should be >= 0",
            eq.cancer_wbc_sep
        );
    }

    /// At zero hematocrit, WBC CFL correction is minimal.
    /// WBC equilibrium should still be bounded.
    #[test]
    fn zero_hematocrit_wbc_bounded() {
        let eq = three_population_equilibria(WIDTH, HEIGHT, Q, RHO, MU, 0.0, None);
        assert!(
            eq.wbc_eq >= 0.0 && eq.wbc_eq <= 1.0,
            "WBC x̃ = {} out of [0,1] at HCT=0",
            eq.wbc_eq
        );
    }

    #[test]
    fn checked_three_population_rejects_nonphysical_geometry() {
        let err = checked_three_population_equilibria(0.0, HEIGHT, Q, RHO, MU, HCT, None)
            .expect_err("checked three-population analysis must reject zero width");
        assert!(err.to_string().contains("geometry"));
    }

    #[test]
    fn checked_three_population_matches_legacy_nominal_case() {
        let legacy = three_population_equilibria(WIDTH, HEIGHT, Q, RHO, MU, HCT, None);
        let checked = checked_three_population_equilibria(WIDTH, HEIGHT, Q, RHO, MU, HCT, None)
            .expect("checked three-population analysis should succeed");

        assert!((legacy.cancer_eq - checked.cancer_eq).abs() < 1e-12);
        assert!((legacy.wbc_eq - checked.wbc_eq).abs() < 1e-12);
        assert!((legacy.rbc_eq - checked.rbc_eq).abs() < 1e-12);
    }
}
