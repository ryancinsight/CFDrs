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
//! let model = CellSeparationModel::new(500e-6, 200e-6, None);
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
pub mod cell_interaction;
pub mod margination;
pub mod properties;
pub mod separation_model;

pub use cascade_junction::{
    cascade_junction_separation, cascade_junction_separation_cross_junction,
    cascade_junction_separation_from_qfracs, cif_pretri_stage_center_fracs,
    cif_pretri_stage_q_fracs, incremental_filtration_separation_cross_junction,
    incremental_filtration_separation_from_qfracs, incremental_filtration_separation_staged,
    mixed_cascade_separation, mixed_cascade_separation_kappa_aware, tri_asymmetric_q_fracs,
    tri_center_q_frac, tri_center_q_frac_cross_junction, CascadeJunctionResult, CascadeStage,
    IncrementalFiltrationResult,
};
pub use cell_interaction::enhanced_lateral_equilibrium;
pub use margination::{
    dean_drag_force_n, dean_number, inertial_lift_force_n, lateral_equilibrium, EquilibriumResult,
};
pub use properties::CellProperties;
pub use separation_model::{CellSeparationAnalysis, CellSeparationModel};

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
    use cell_interaction::enhanced_lateral_equilibrium;

    let cancer = CellProperties::mcf7_breast_cancer();
    let wbc = CellProperties::white_blood_cell();
    let rbc = CellProperties::red_blood_cell();

    let area = (width_m * height_m).max(1e-18);
    let mean_v = flow_rate_m3_s / area;

    const DISPERSED: f64 = 0.5;
    const SPLIT: f64 = 0.3;

    let cancer_x = margination::lateral_equilibrium(
        &cancer,
        blood_density,
        viscosity,
        mean_v,
        width_m,
        height_m,
        bend_radius_m,
    )
    .map_or(DISPERSED, |r| r.x_tilde_eq);

    let wbc_x = enhanced_lateral_equilibrium(
        &wbc,
        blood_density,
        viscosity,
        mean_v,
        width_m,
        height_m,
        bend_radius_m,
        hematocrit,
    )
    .map_or(DISPERSED, |r| r.x_tilde_eq);

    let rbc_x = margination::lateral_equilibrium(
        &rbc,
        blood_density,
        viscosity,
        mean_v,
        width_m,
        height_m,
        bend_radius_m,
    )
    .map_or(DISPERSED, |r| r.x_tilde_eq);

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

    ThreePopEquilibria {
        cancer_eq: cancer_x,
        wbc_eq: wbc_x,
        rbc_eq: rbc_x,
        cancer_wbc_sep,
        rbc_vs_center_sep,
        cancer_center_fraction,
        wbc_center_fraction,
        rbc_peripheral_fraction,
    }
}
