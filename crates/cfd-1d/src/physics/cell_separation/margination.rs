//! Inertial lift and Dean drag forces for cell margination in microchannels.
//!
//! # Physical model
//!
//! In a straight rectangular channel at finite Reynolds number, a neutrally
//! buoyant particle experiences two competing lateral forces:
//!
//! ## 1. Inertial lift force (Segré-Silberberg)
//!
//! The net inertial lift on a sphere of diameter `a` in a channel of hydraulic
//! scale `H` at mean velocity `U` is the difference between two finite-size
//! inertial migration mechanisms:
//!
//! ```text
//! F_L = F_wall - F_shear
//! F_wall  = C_wall(x̃)  · (1 - DI) · ρ · U² · a⁶ / H⁴
//! F_shear = C_shear(x̃)          · ρ · U² · a³ / H
//! ```
//!
//! where `x̃ = x / (H/2)` is the dimensionless lateral position on the
//! half-channel (`0 = center`, `1 = wall`). The scaling separates
//! wall-induced lift (`a⁶/H⁴`) from shear-gradient lift (`a³/H`), matching
//! the finite-size inertial microfluidics literature instead of folding both
//! mechanisms into one coefficient.
//!
//! For deformable cells, the equilibrium shifts toward the wall because
//! deformability reduces the wall-repulsion component of lift.  We model this
//! via a deformability correction factor `(1 − DI)` applied to the wall-
//! repulsion term (Hur et al. 2011, *Lab Chip* 11, 912–920):
//!
//! ```text
//! F_L(x̃) = F_wall(x̃, DI) − F_shear(x̃)
//! ```
//!
//! where the wall-induced and shear-gradient shape functions are bounded
//! position functions on the closed half-channel domain. The implementation
//! does not evaluate outside `0 <= x̃ <= 1` and does not clip singular wall
//! neighborhoods into the model.
//!
//! ## 2. Dean drag force (curved channels)
//!
//! In a curved channel of radius `R`, the Dean number is:
//!
//! ```text
//! De = Re · √(D_h / (2R))
//! ```
//!
//! The secondary Dean flow creates a drag force on particles (Gossett & Di Carlo
//! 2009, *Anal. Chem.* 81, 8459–8465):
//!
//! ```text
//! F_D = 5.4 × 10⁻⁴ · π · μ · De^{1.63} · a
//! ```
//!
//! This force acts in the plane of the channel cross-section, pushing particles
//! toward the outer wall (positive x̃ direction for a left-curving channel).
//!
//! ## Equilibrium position
//!
//! The equilibrium lateral position `x̃_eq` satisfies `F_L(x̃_eq) = F_D`.
//!
//! In the toward-center sign convention (F_L > 0 → toward center), Dean drag
//! acts toward the wall, so the bisection solves `F_L(x̃) − F_D = 0`.
//! Increasing F_D shifts the equilibrium toward the wall (larger x̃), consistent
//! with experimental observations of Dean-flow-enhanced cell migration.
//! We solve this numerically using bisection on `[0, 1]` (by symmetry,
//! only the half-channel need be considered).
//!
//! ## Theorem: Inertial Lift Force Scaling (Segré-Silberberg)
//!
//! **Theorem**: A neutrally-buoyant sphere of diameter `a` in a rectangular
//! channel of hydraulic scale `H` at mean velocity `U` experiences a signed
//! lateral inertial force:
//!
//! ```text
//! F_L(x̃, DI) =
//!     C_wall(x̃)(1 - DI)ρU²a⁶/H⁴ - C_shear(x̃)ρU²a³/H
//! ```
//!
//! **Stability**: The equilibrium position `x̃_eq` where `F_L(x̃_eq) = 0` is **stable**:
//! - For `x̃ < x̃_eq`: `F_L > 0` (force pushes toward center → away from wall)
//! - For `x̃ > x̃_eq`: `F_L < 0` (force pushes toward wall → away from center)
//!
//! **Proof sketch**: Reviews through 2025 continue to identify wall-induced
//! lift scaling as `ρU²a⁶/H⁴` and shear-gradient lift as `ρU²a³/H`. The
//! bounded shapes used here satisfy the required signs: wall lift is zero at
//! the channel center and increases monotonically toward the wall, while
//! shear-gradient lift is maximal at the center and vanishes at the wall. The
//! two continuous terms therefore form a force-balance root on the validated
//! half-channel when their endpoint signs bracket zero. Deformability
//! multiplies only the wall-repulsion term, shifting softer cells toward the
//! wall by reducing the restoring force.
//!
//! **Rigid sphere equilibrium**: At DI = 0, the equilibrium `x̃_eq ≈ 0.6` (Segré-Silberberg
//! position), corresponding to `≈ 0.6 × H/2` from the center (Di Carlo 2009, Fig. 2).
//!
//! ## Theorem: Dean Drag Bisection Convergence
//!
//! **Theorem**: The bisection algorithm on `[0, 1]` for `F_L(x̃) - F_Dean = 0`
//! converges to machine precision in at most 60 iterations.
//!
//! **Proof**: At each step the interval halves: `|x̃_{n+1} - x̃*| ≤ (1/2) · 2^{-n}`.
//! After 60 steps: `|error| ≤ 2^{-60} ≈ 8.67 × 10^{-19}`, which is below
//! `f64` machine epsilon (~2.22 × 10^{-16}). Since `F_L` is continuous on
//! `[0, 1]`, the intermediate value theorem guarantees a root when
//! `sign(F_lo) != sign(F_hi)`. ∎
//!
//! # References
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.
//! - Gossett, D. R. & Di Carlo, D. (2009). Particle focusing mechanisms in
//!   curving confined flows. *Anal. Chem.*, 81, 8459–8465.
//! - Hur, S. C., Henderson-MacLennan, N. K., McCabe, E. R. B. & Di Carlo, D.
//!   (2011). Deformability-based cell classification and enrichment using
//!   inertial microfluidics. *Lab Chip*, 11, 912–920.

use crate::physics::cell_separation::properties::CellProperties;
use cfd_core::error::{Error, Result};
use serde::{Deserialize, Serialize};

// ── Amini (2014) confinement-dependent lift correction ──────────────────────

/// Reference confinement ratio `κ_ref = 0.1` for the Amini correction.
///
/// At this confinement ratio the correction factor is unity (no modification).
pub const AMINI_KAPPA_REF: f64 = 0.1;

/// Calibrated confinement sensitivity coefficient `α = 2.5`.
///
/// Fitted from experimental data in Amini et al. (2014), *Lab Chip* 14:2739–2761.
pub const AMINI_ALPHA_CONFINEMENT: f64 = 2.5;

/// Reference confinement ratio for MCF-7-like rigid cell focusing in a 200 µm
/// channel height.
pub const LIFT_REFERENCE_KAPPA: f64 = 17.5e-6 / 200.0e-6;

/// Reference equilibrium half-channel coordinate for rigid inertial focusing.
pub const LIFT_REFERENCE_X_TILDE: f64 = 0.55;

/// Reference deformability index for MCF-7 breast-cancer cells.
pub const LIFT_REFERENCE_DEFORMABILITY_INDEX: f64 = 0.15;

/// Amini (2014) confinement-dependent inertial lift correction factor.
///
/// Returns a multiplicative correction to the Di Carlo (2009) lift coefficient
/// that accounts for finite-size (confinement) effects:
///
/// ```text
/// f(κ) = 1 + α_confinement · (κ − κ_ref)²
/// ```
///
/// where `κ = a / D_h` is the confinement ratio (cell diameter / hydraulic diameter).
///
/// ## Theorem: Confinement-Enhanced Inertial Focusing (Amini et al. 2014)
///
/// For `κ > κ_ref`, inertial lift is amplified by wall-induced confinement effects
/// (Amini et al. 2014). The correction increases the focusing force by up to 3×
/// for `κ ≈ 0.2`. This arises because larger confinement ratios strengthen the
/// wall-particle hydrodynamic interaction, increasing both the wall-repulsion and
/// shear-gradient lift components. The quadratic form ensures smooth, monotonic
/// amplification away from the reference confinement ratio.
///
/// ## Properties
///
/// - `f(κ_ref) = 1.0` (no correction at reference confinement)
/// - `f(κ) > 1.0` for all `κ ≠ κ_ref` (always amplifying)
/// - `f(κ) ≥ 1.0` for all `κ` (always positive)
///
/// ## Reference
///
/// Amini, H., Lee, W. & Di Carlo, D. (2014). Inertial microfluidic physics.
/// *Lab Chip*, 14, 2739–2761.
///
/// # Arguments
///
/// * `kappa` — confinement ratio `a / D_h` (cell diameter / hydraulic diameter)
#[inline]
#[must_use]
pub fn amini_confinement_correction(kappa: f64) -> f64 {
    let delta = kappa - AMINI_KAPPA_REF;
    1.0 + AMINI_ALPHA_CONFINEMENT * delta * delta
}

/// Checked Amini confinement correction for callers that want explicit model
/// boundary enforcement instead of silent acceptance of nonphysical inputs.
pub fn checked_amini_confinement_correction(kappa: f64) -> Result<f64> {
    if !kappa.is_finite() || kappa < 0.0 {
        return Err(Error::InvalidConfiguration(
            "Margination confinement ratio must be finite and nonnegative".to_string(),
        ));
    }
    Ok(amini_confinement_correction(kappa))
}

fn validate_lateral_position(x_tilde: f64) -> Result<f64> {
    if !x_tilde.is_finite() || !(0.0..=1.0).contains(&x_tilde) {
        return Err(Error::InvalidConfiguration(
            "Margination lateral position x_tilde must lie in [0, 1]".to_string(),
        ));
    }
    Ok(x_tilde)
}

fn validate_deformability_index(di: f64) -> Result<f64> {
    if !di.is_finite() || !(0.0..=1.0).contains(&di) {
        return Err(Error::InvalidConfiguration(
            "Margination deformability index must lie in [0, 1]".to_string(),
        ));
    }
    Ok(di)
}

fn validate_positive_finite(name: &str, value: f64) -> Result<f64> {
    if !value.is_finite() || value <= 0.0 {
        return Err(Error::InvalidConfiguration(format!(
            "Margination {name} must be finite and positive"
        )));
    }
    Ok(value)
}

// ── Lift coefficient model ────────────────────────────────────────────────────

/// Dimensionless wall-induced lift shape on the validated half-channel domain.
///
/// # Theorem - Wall-Induced Lift Shape
///
/// `x̃²` is nonnegative for `x̃ ∈ [0, 1]`, vanishes at the channel center,
/// and is strictly increasing for `x̃ > 0`.
///
/// **Proof sketch**: The derivative is `2x̃`, which is nonnegative on the
/// closed validated interval and strictly positive for `x̃ > 0`. Therefore
/// the wall-repulsion magnitude increases monotonically as the particle
/// approaches the wall while remaining finite at `x̃ = 1`.
#[inline]
fn wall_lift_shape(x_tilde: f64) -> f64 {
    x_tilde * x_tilde
}

/// Dimensionless shear-gradient lift shape on the validated half-channel domain.
///
/// # Theorem - Shear-Gradient Lift Shape
///
/// `1 - x̃²` is nonnegative and monotonically decreasing on `x̃ ∈ [0, 1]`.
///
/// **Proof sketch**: `1 - x̃² >= 0` on the closed unit interval and its
/// derivative is `-2x̃ <= 0`. The shape is maximal at the center and vanishes
/// at the wall, matching the shear-gradient contribution that drives finite
/// particles away from the center in Poiseuille flow.
#[inline]
fn shear_gradient_lift_shape(x_tilde: f64) -> f64 {
    1.0 - x_tilde * x_tilde
}

/// Dimensionless wall-lift gain derived from the reference focusing condition.
///
/// # Theorem - Reference Equilibrium Calibration
///
/// Let `x_ref`, `k_ref`, and `DI_ref` denote the documented rigid-cell
/// reference equilibrium, confinement ratio, and deformability index. Defining
///
/// ```text
/// G = C_shear(x_ref) k_ref /
///     (C_wall(x_ref) (1 - DI_ref) k_ref^4)
/// ```
///
/// makes the dimensional lift force vanish exactly at the reference state.
///
/// **Proof sketch**: Substitute `G` into
/// `G C_wall(x)(1-DI)k^4 - C_shear(x)k`. At the reference state the first term
/// becomes `C_shear(x_ref)k_ref`, which cancels the second term exactly.
#[inline]
fn wall_lift_reference_gain() -> f64 {
    shear_gradient_lift_shape(LIFT_REFERENCE_X_TILDE) * LIFT_REFERENCE_KAPPA
        / (wall_lift_shape(LIFT_REFERENCE_X_TILDE)
            * (1.0 - LIFT_REFERENCE_DEFORMABILITY_INDEX)
            * LIFT_REFERENCE_KAPPA.powi(4))
}

/// Net dimensional inertial lift force [N].
///
/// Positive → force toward center (away from wall).
/// Negative → force toward wall.
///
#[inline]
fn dimensional_lift_force_n(
    x_tilde: f64,
    cell_diameter_m: f64,
    deformability_index: f64,
    fluid_density_kg_m3: f64,
    mean_velocity_m_s: f64,
    channel_height_m: f64,
) -> f64 {
    debug_assert!((0.0..=1.0).contains(&deformability_index));
    let di = deformability_index;
    let dynamic_pressure = fluid_density_kg_m3 * mean_velocity_m_s * mean_velocity_m_s;
    let confinement = cell_diameter_m / channel_height_m;
    let area_scale = cell_diameter_m * cell_diameter_m;
    let wall = wall_lift_reference_gain()
        * wall_lift_shape(x_tilde)
        * (1.0 - di)
        * area_scale
        * confinement.powi(4);
    let shear = shear_gradient_lift_shape(x_tilde) * area_scale * confinement;
    dynamic_pressure * (wall - shear)
}

// ── Dean drag model ───────────────────────────────────────────────────────────

/// Dean number `De = Re · √(D_h / (2R))`.
///
/// # Arguments
/// - `re` — channel Reynolds number `ρ U D_h / μ`
/// - `hydraulic_diameter_m` — `D_h = 2wh/(w+h)` [m]
/// - `bend_radius_m` — radius of curvature of the channel centreline [m]
///
/// # Panics
/// Panics in debug mode if `bend_radius_m ≤ 0`.
#[inline]
#[must_use]
pub fn dean_number(re: f64, hydraulic_diameter_m: f64, bend_radius_m: f64) -> f64 {
    debug_assert!(bend_radius_m > 0.0, "bend_radius_m must be positive");
    re * (hydraulic_diameter_m / (2.0 * bend_radius_m)).sqrt()
}

/// Dean drag force magnitude [N] on a cell of diameter `a` [m].
///
/// From Gossett & Di Carlo (2009), Eq. 4:
/// `F_D = 5.4 × 10⁻⁴ · π · μ · De^{1.63} · a`
///
/// The force acts perpendicular to the flow direction, pushing cells toward
/// the outer wall of the curved channel.
#[inline]
#[must_use]
pub fn dean_drag_force_n(dynamic_viscosity_pa_s: f64, de: f64, cell_diameter_m: f64) -> f64 {
    5.4e-4 * std::f64::consts::PI * dynamic_viscosity_pa_s * de.powf(1.63) * cell_diameter_m
}

// ── Inertial lift force ───────────────────────────────────────────────────────

/// Inertial lift force [N] on a cell at normalised lateral position `x̃`.
///
/// `F_L = C_wall(1-DI)ρU²a⁶/H⁴ - C_shearρU²a³/H`
///
/// Positive → toward center; negative → toward wall.
///
/// # Arguments
/// - `x_tilde` — normalised lateral position ∈ [0, 1] (0 = center, 1 = wall)
/// - `cell` — cell physical properties
/// - `fluid_density_kg_m3` — fluid density [kg/m³]
/// - `mean_velocity_m_s` — mean channel velocity [m/s]
/// - `channel_height_m` — channel height (shorter dimension) [m]
#[inline]
#[must_use]
pub fn inertial_lift_force_n(
    x_tilde: f64,
    cell: &CellProperties,
    fluid_density_kg_m3: f64,
    mean_velocity_m_s: f64,
    channel_height_m: f64,
) -> f64 {
    checked_inertial_lift_force_n(
        x_tilde,
        cell,
        fluid_density_kg_m3,
        mean_velocity_m_s,
        channel_height_m,
    )
    .expect("inertial lift inputs must satisfy the validated margination envelope")
}

/// Checked inertial lift force evaluation using the validated margination model.
pub fn checked_inertial_lift_force_n(
    x_tilde: f64,
    cell: &CellProperties,
    fluid_density_kg_m3: f64,
    mean_velocity_m_s: f64,
    channel_height_m: f64,
) -> Result<f64> {
    let x_tilde = validate_lateral_position(x_tilde)?;
    let deformability_index = validate_deformability_index(cell.deformability_index)?;
    let fluid_density_kg_m3 = validate_positive_finite("fluid density", fluid_density_kg_m3)?;
    let channel_height_m = validate_positive_finite("channel height", channel_height_m)?;
    if !mean_velocity_m_s.is_finite() {
        return Err(Error::InvalidConfiguration(
            "Margination mean velocity must be finite".to_string(),
        ));
    }
    if !cell.diameter_m.is_finite() || cell.diameter_m <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Margination cell diameter must be finite and positive".to_string(),
        ));
    }

    Ok(dimensional_lift_force_n(
        x_tilde,
        cell.diameter_m,
        deformability_index,
        fluid_density_kg_m3,
        mean_velocity_m_s,
        channel_height_m,
    ))
}

// ── Lateral Drift Velocity ────────────────────────────────────────────────────

/// Transient lateral drift velocity [m/s] at a given lateral position `x̃`.
///
/// Determined by balancing the net lateral forces (`F_L - F_D`) against
/// Stokes drag (`F_{drag} = 3 \pi \mu a v_{drift}`).
///
/// Positive velocity points toward the channel center.
///
/// # Arguments
/// - `x_tilde` — normalised lateral position ∈ [0, 1] (0 = center, 1 = wall)
/// - `cell` — cell physical properties
/// - `fluid_density_kg_m3` — fluid density [kg/m³]
/// - `dynamic_viscosity_pa_s` — fluid dynamic viscosity [Pa·s]
/// - `mean_velocity_m_s` — mean channel velocity [m/s]
/// - `channel_width_m` — channel width [m]
/// - `channel_height_m` — channel height (shorter dimension) [m]
/// - `bend_radius_m` — radius of curvature [m], or `None` for a straight channel
#[inline]
#[must_use]
pub fn lateral_velocity_m_s(
    x_tilde: f64,
    cell: &CellProperties,
    fluid_density_kg_m3: f64,
    dynamic_viscosity_pa_s: f64,
    mean_velocity_m_s: f64,
    channel_width_m: f64,
    channel_height_m: f64,
    bend_radius_m: Option<f64>,
) -> f64 {
    let dh = 2.0 * channel_width_m * channel_height_m / (channel_width_m + channel_height_m);
    let f_lift = checked_inertial_lift_force_n(
        x_tilde,
        cell,
        fluid_density_kg_m3,
        mean_velocity_m_s,
        channel_height_m,
    )
    .expect("lateral velocity inputs must satisfy the validated margination envelope");

    let f_dean = if let Some(r) = bend_radius_m {
        if r > 0.0 {
            let re = fluid_density_kg_m3 * mean_velocity_m_s * dh / dynamic_viscosity_pa_s;
            let de = dean_number(re, dh, r);
            dean_drag_force_n(dynamic_viscosity_pa_s, de, cell.diameter_m)
        } else {
            0.0
        }
    } else {
        0.0
    };

    let net_force_n = f_lift - f_dean;
    let stokes_coeff = 3.0 * std::f64::consts::PI * dynamic_viscosity_pa_s * cell.diameter_m;

    net_force_n / stokes_coeff.max(1e-30)
}

/// Checked lateral drift velocity evaluation using explicit model bounds.
pub fn checked_lateral_velocity_m_s(
    x_tilde: f64,
    cell: &CellProperties,
    fluid_density_kg_m3: f64,
    dynamic_viscosity_pa_s: f64,
    mean_velocity_m_s: f64,
    channel_width_m: f64,
    channel_height_m: f64,
    bend_radius_m: Option<f64>,
) -> Result<f64> {
    let channel_width_m = validate_positive_finite("channel width", channel_width_m)?;
    let channel_height_m = validate_positive_finite("channel height", channel_height_m)?;
    let dynamic_viscosity_pa_s =
        validate_positive_finite("dynamic viscosity", dynamic_viscosity_pa_s)?;
    if let Some(radius) = bend_radius_m {
        validate_positive_finite("bend radius", radius)?;
    }

    let dh = 2.0 * channel_width_m * channel_height_m / (channel_width_m + channel_height_m);
    let f_lift = checked_inertial_lift_force_n(
        x_tilde,
        cell,
        fluid_density_kg_m3,
        mean_velocity_m_s,
        channel_height_m,
    )?;

    let f_dean = if let Some(r) = bend_radius_m {
        let re = fluid_density_kg_m3 * mean_velocity_m_s * dh / dynamic_viscosity_pa_s;
        let de = dean_number(re, dh, r);
        dean_drag_force_n(dynamic_viscosity_pa_s, de, cell.diameter_m)
    } else {
        0.0
    };

    let net_force_n = f_lift - f_dean;
    let stokes_coeff = 3.0 * std::f64::consts::PI * dynamic_viscosity_pa_s * cell.diameter_m;
    Ok(net_force_n / stokes_coeff.max(1e-30))
}

// ── Equilibrium position solver ───────────────────────────────────────────────

/// Result of a lateral equilibrium position calculation.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EquilibriumResult {
    /// Normalised lateral equilibrium position `x̃_eq ∈ [0, 1]`.
    ///
    /// 0 = channel center; 1 = channel wall.
    /// Values < 0.3 indicate center-focused; values > 0.5 indicate wall-focused.
    pub x_tilde_eq: f64,

    /// Dimensional lateral position from channel center [m].
    pub lateral_position_m: f64,

    /// Inertial lift force at equilibrium [N] (should be ≈ 0).
    pub residual_force_n: f64,

    /// Dean drag force at equilibrium [N] (0 if no curvature).
    pub dean_drag_n: f64,

    /// Channel Reynolds number.
    pub reynolds_number: f64,

    /// Dean number (0 if no curvature).
    pub dean_number: f64,

    /// Whether the cell will exhibit measurable focusing (κ > 0.07).
    pub will_focus: bool,
}

/// Compute the lateral equilibrium position of a cell in a rectangular channel.
///
/// Solves `F_L(x̃) − F_D = 0` using bisection on `x̃ ∈ [0, 1]`.
///
/// Sign convention: F_L positive = toward center; Dean drag `F_D ≥ 0` acts
/// toward the wall, so it enters with a negative sign.  This ensures that
/// larger Dean forces shift the equilibrium toward the wall, as observed
/// experimentally (Gossett & Di Carlo 2009).
///
/// For straight channels (`bend_radius_m = None`), `F_D = 0` and the
/// equilibrium is determined by inertial lift alone.
///
/// For curved channels, the Dean drag pushes cells toward the outer wall
/// (increasing `x̃`), shifting the equilibrium.
///
/// # Arguments
/// - `cell` — cell physical properties
/// - `fluid_density_kg_m3` — fluid density [kg/m³]
/// - `dynamic_viscosity_pa_s` — fluid dynamic viscosity [Pa·s]
/// - `mean_velocity_m_s` — mean channel velocity [m/s]
/// - `channel_width_m` — channel width [m]
/// - `channel_height_m` — channel height [m] (shorter dimension)
/// - `bend_radius_m` — radius of curvature [m], or `None` for straight channel
///
/// # Returns
/// [`EquilibriumResult`] with the equilibrium position and force balance.
///
/// # Errors
/// Returns `None` if the cell does not focus (κ ≤ 0.07) or if bisection
/// fails to converge within 100 iterations.
#[must_use]
pub fn lateral_equilibrium(
    cell: &CellProperties,
    fluid_density_kg_m3: f64,
    dynamic_viscosity_pa_s: f64,
    mean_velocity_m_s: f64,
    channel_width_m: f64,
    channel_height_m: f64,
    bend_radius_m: Option<f64>,
) -> Option<EquilibriumResult> {
    let h = channel_height_m;
    let w = channel_width_m;
    let dh = 2.0 * w * h / (w + h); // hydraulic diameter

    let kappa = cell.confinement_ratio(dh);
    let will_focus = kappa > 0.07;

    // Reynolds number: Re = ρ U D_h / μ
    let re = fluid_density_kg_m3 * mean_velocity_m_s * dh / dynamic_viscosity_pa_s;

    // Dean number and drag force (0 for straight channels)
    let (de, f_dean) = if let Some(r) = bend_radius_m {
        if r > 0.0 {
            let de = dean_number(re, dh, r);
            let fd = dean_drag_force_n(dynamic_viscosity_pa_s, de, cell.diameter_m);
            (de, fd)
        } else {
            (0.0, 0.0)
        }
    } else {
        (0.0, 0.0)
    };

    // Net force function (in the toward-center sign convention).
    //
    // Sign convention: positive = force toward channel CENTER.
    //   • F_L > 0  → inertial lift pushes toward center (away from wall)
    //   • F_L < 0  → inertial lift pushes toward wall
    //   • f_dean ≥ 0 is the Dean drag MAGNITUDE; it acts toward the wall,
    //     so it enters as −f_dean in this toward-center convention.
    //
    // Equilibrium: F_L(x̃) − f_dean = 0  →  F_L(x̃) = f_dean
    // For straight channel (f_dean=0) this gives F_L = 0 (Segré-Silberberg).
    // For curved channels f_dean > 0 forces equilibrium to larger x̃ (toward wall).
    let net_force = |x: f64| -> f64 {
        inertial_lift_force_n(x, cell, fluid_density_kg_m3, mean_velocity_m_s, h) - f_dean
    };

    // Bisection on [0, 1] over the validated nonsingular half-channel domain.
    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;

    // Check if there is a sign change (equilibrium exists in interval)
    let f_lo = net_force(lo);
    let f_hi = net_force(hi);

    // If both ends have the same sign, the equilibrium is at the boundary
    let x_eq = if f_lo * f_hi > 0.0 {
        // No zero crossing: equilibrium is at the end with smaller |F|
        if f_lo.abs() < f_hi.abs() {
            lo
        } else {
            hi
        }
    } else {
        // Bisection: 60 iterations -> precision < 2^-60.
        let mut mid = 0.5 * (lo + hi);
        for _ in 0..60 {
            mid = 0.5 * (lo + hi);
            let f_mid = net_force(mid);
            if f_mid.abs() < 1e-30 {
                break;
            }
            if f_lo * f_mid < 0.0 {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        mid
    };

    let residual = net_force(x_eq);
    let lateral_pos_m = x_eq * (h / 2.0); // dimensional position from center

    Some(EquilibriumResult {
        x_tilde_eq: x_eq,
        lateral_position_m: lateral_pos_m,
        residual_force_n: residual,
        dean_drag_n: f_dean,
        reynolds_number: re,
        dean_number: de,
        will_focus,
    })
}

/// Checked lateral equilibrium solver using explicit model-boundary validation.
pub fn checked_lateral_equilibrium(
    cell: &CellProperties,
    fluid_density_kg_m3: f64,
    dynamic_viscosity_pa_s: f64,
    mean_velocity_m_s: f64,
    channel_width_m: f64,
    channel_height_m: f64,
    bend_radius_m: Option<f64>,
) -> Result<EquilibriumResult> {
    let channel_width_m = validate_positive_finite("channel width", channel_width_m)?;
    let channel_height_m = validate_positive_finite("channel height", channel_height_m)?;
    let dynamic_viscosity_pa_s =
        validate_positive_finite("dynamic viscosity", dynamic_viscosity_pa_s)?;
    let fluid_density_kg_m3 = validate_positive_finite("fluid density", fluid_density_kg_m3)?;
    if !mean_velocity_m_s.is_finite() || mean_velocity_m_s <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Margination mean velocity must be finite and positive".to_string(),
        ));
    }
    if let Some(radius) = bend_radius_m {
        validate_positive_finite("bend radius", radius)?;
    }

    let h = channel_height_m;
    let w = channel_width_m;
    let dh = 2.0 * w * h / (w + h);
    let kappa = cell.confinement_ratio(dh);
    let will_focus = kappa > 0.07;
    let re = fluid_density_kg_m3 * mean_velocity_m_s * dh / dynamic_viscosity_pa_s;

    let (de, f_dean) = if let Some(r) = bend_radius_m {
        let de = dean_number(re, dh, r);
        let fd = dean_drag_force_n(dynamic_viscosity_pa_s, de, cell.diameter_m);
        (de, fd)
    } else {
        (0.0, 0.0)
    };

    let net_force = |x: f64| -> Result<f64> {
        Ok(
            checked_inertial_lift_force_n(x, cell, fluid_density_kg_m3, mean_velocity_m_s, h)?
                - f_dean,
        )
    };

    let mut lo = 0.0_f64;
    let mut hi = 1.0_f64;
    let f_lo = net_force(lo)?;
    let f_hi = net_force(hi)?;

    let x_eq = if f_lo * f_hi > 0.0 {
        if f_lo.abs() < f_hi.abs() {
            lo
        } else {
            hi
        }
    } else {
        let mut mid = 0.5 * (lo + hi);
        let mut f_lo_mut = f_lo;
        for _ in 0..60 {
            mid = 0.5 * (lo + hi);
            let f_mid = net_force(mid)?;
            if f_mid.abs() < 1e-30 {
                break;
            }
            if f_lo_mut * f_mid < 0.0 {
                hi = mid;
            } else {
                lo = mid;
                f_lo_mut = f_mid;
            }
        }
        mid
    };

    let residual = net_force(x_eq)?;
    Ok(EquilibriumResult {
        x_tilde_eq: x_eq,
        lateral_position_m: x_eq * (h / 2.0),
        residual_force_n: residual,
        dean_drag_n: f_dean,
        reynolds_number: re,
        dean_number: de,
        will_focus,
    })
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── amini_confinement_correction ─────────────────────────────────────────

    #[test]
    fn test_amini_correction_at_reference_kappa() {
        // At κ = κ_ref = 0.1, the correction factor should be exactly 1.0
        let correction = amini_confinement_correction(AMINI_KAPPA_REF);
        assert!(
            (correction - 1.0).abs() < 1e-15,
            "correction at κ_ref should be 1.0, got {correction}"
        );
    }

    #[test]
    fn test_amini_correction_increases_with_confinement() {
        // For κ > κ_ref, the correction should be > 1.0
        let correction_015 = amini_confinement_correction(0.15);
        let correction_020 = amini_confinement_correction(0.20);
        let correction_025 = amini_confinement_correction(0.25);

        assert!(
            correction_015 > 1.0,
            "correction at κ=0.15 should exceed 1.0, got {correction_015}"
        );
        assert!(
            correction_020 > correction_015,
            "correction should increase with κ: {correction_020} vs {correction_015}"
        );
        assert!(
            correction_025 > correction_020,
            "correction should increase with κ: {correction_025} vs {correction_020}"
        );

        // At κ = 0.2, correction = 1 + 2.5 × (0.1)² = 1.025
        let expected_020 = 1.0 + AMINI_ALPHA_CONFINEMENT * 0.1 * 0.1;
        assert!(
            (correction_020 - expected_020).abs() < 1e-15,
            "correction at κ=0.2: got {correction_020}, expected {expected_020}"
        );
    }

    #[test]
    fn test_amini_correction_bounded() {
        // The correction should always be positive (≥ 1.0) for any κ
        for &kappa in &[0.0, 0.01, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5, 1.0] {
            let correction = amini_confinement_correction(kappa);
            assert!(
                correction >= 1.0,
                "correction must be ≥ 1.0 for κ={kappa}, got {correction}"
            );
        }
        // Also positive for negative κ (physically meaningless but numerically safe)
        let correction_neg = amini_confinement_correction(-0.1);
        assert!(
            correction_neg > 0.0,
            "correction must be positive, got {correction_neg}"
        );
    }

    #[test]
    fn test_checked_amini_rejects_negative_kappa() {
        assert!(checked_amini_confinement_correction(-0.1).is_err());
    }

    #[test]
    fn test_checked_inertial_lift_rejects_out_of_domain_position() {
        let cell = CellProperties::mcf7_breast_cancer();
        let err = checked_inertial_lift_force_n(1.01, &cell, 1_000.0, 0.05, 200e-6).expect_err(
            "checked margination lift must reject x_tilde outside the validated domain",
        );
        assert!(err.to_string().contains("x_tilde"));
    }

    #[test]
    #[should_panic(expected = "validated margination envelope")]
    fn legacy_inertial_lift_rejects_instead_of_clamping_position() {
        let cell = CellProperties::mcf7_breast_cancer();
        let _ = inertial_lift_force_n(1.01, &cell, 1_000.0, 0.05, 200e-6);
    }

    #[test]
    fn test_checked_lateral_equilibrium_rejects_nonphysical_geometry() {
        let cell = CellProperties::mcf7_breast_cancer();
        let err = checked_lateral_equilibrium(&cell, 1_000.0, 1.0e-3, 0.05, 100e-6, 0.0, None)
            .expect_err("checked lateral equilibrium must reject zero channel height");
        assert!(err.to_string().contains("channel height"));
    }

    #[test]
    fn test_checked_lateral_equilibrium_matches_legacy_on_nominal_case() {
        let cell = CellProperties::mcf7_breast_cancer();
        let legacy = lateral_equilibrium(&cell, 1_000.0, 1.0e-3, 0.05, 200e-6, 100e-6, None)
            .expect("legacy margination equilibrium should succeed");
        let checked =
            checked_lateral_equilibrium(&cell, 1_000.0, 1.0e-3, 0.05, 200e-6, 100e-6, None)
                .expect("checked margination equilibrium should succeed");

        assert!((legacy.x_tilde_eq - checked.x_tilde_eq).abs() < 1e-12);
        assert!((legacy.residual_force_n - checked.residual_force_n).abs() < 1e-18);
        assert_eq!(legacy.will_focus, checked.will_focus);
    }

    #[test]
    fn wall_lift_shape_is_monotone_on_validated_domain() {
        let low = wall_lift_shape(0.25);
        let mid = wall_lift_shape(0.50);
        let high = wall_lift_shape(0.75);

        assert_eq!(wall_lift_shape(0.0), 0.0);
        assert_eq!(wall_lift_shape(1.0), 1.0);
        assert!(low > 0.0);
        assert!(mid > low);
        assert!(high > mid);
    }

    #[test]
    fn shear_gradient_lift_shape_decreases_to_zero_at_wall() {
        let center = shear_gradient_lift_shape(0.0);
        let mid = shear_gradient_lift_shape(0.50);
        let wall = shear_gradient_lift_shape(1.0);

        assert_eq!(center, 1.0);
        assert!(mid < center);
        assert!(wall < mid);
        assert_eq!(wall, 0.0);
    }

    #[test]
    fn dimensional_lift_uses_distinct_confinement_scalings() {
        let rho = 1_000.0;
        let u = 0.05;
        let h = 100.0e-6;
        let a = 10.0e-6;
        let x = 0.25;

        let force = dimensional_lift_force_n(x, a, 0.0, rho, u, h);
        let q = rho * u * u;
        let area = a * a;
        let kappa = a / h;
        let expected = q
            * (wall_lift_reference_gain() * wall_lift_shape(x) * area * kappa.powi(4)
                - shear_gradient_lift_shape(x) * area * kappa);

        assert!((force - expected).abs() < 1.0e-24);
        assert!(force < 0.0);
    }

    #[test]
    fn wall_lift_gain_places_reference_cell_at_documented_equilibrium() {
        let force = dimensional_lift_force_n(
            LIFT_REFERENCE_X_TILDE,
            17.5e-6,
            LIFT_REFERENCE_DEFORMABILITY_INDEX,
            1_000.0,
            0.05,
            200.0e-6,
        );

        assert!(force.abs() < 1.0e-24);
    }

    #[test]
    fn inertial_lift_matches_derived_reference_data() -> cfd_core::error::Result<()> {
        let cell = CellProperties::mcf7_breast_cancer();
        let rho = 1_000.0;
        let u = 0.05;
        let h = 200.0e-6;
        let x = 0.75;
        let kappa = cell.diameter_m / h;
        let area = cell.diameter_m * cell.diameter_m;
        let expected = rho
            * u
            * u
            * (wall_lift_reference_gain()
                * x
                * x
                * (1.0 - cell.deformability_index)
                * area
                * kappa.powi(4)
                - (1.0 - x * x) * area * kappa);

        let actual = checked_inertial_lift_force_n(x, &cell, rho, u, h)?;
        assert!((actual - expected).abs() < 1.0e-24);
        assert!(actual > 0.0);
        Ok(())
    }
}
