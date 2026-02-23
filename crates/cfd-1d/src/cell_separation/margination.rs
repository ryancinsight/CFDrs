//! Inertial lift and Dean drag forces for cell margination in microchannels.
//!
//! # Physical model
//!
//! In a straight rectangular channel at finite Reynolds number, a neutrally
//! buoyant particle experiences two competing lateral forces:
//!
//! ## 1. Inertial lift force (Segré-Silberberg)
//!
//! The net inertial lift on a sphere of diameter `a` in a channel of height `H`
//! at mean velocity `U` is (Di Carlo 2009, Eq. 1):
//!
//! ```text
//! F_L = C_L(x̃) · ρ · U² · a⁴ / H²
//! ```
//!
//! where `x̃ = x / (H/2)` is the dimensionless lateral position (−1 = wall,
//! 0 = center, +1 = opposite wall) and `C_L(x̃)` is the position-dependent
//! lift coefficient.  `C_L` is positive (toward center) near the wall and
//! negative (toward wall) near the center, creating a stable equilibrium at
//! `x̃_eq ≈ ±0.6` for rigid spheres in a square channel.
//!
//! For deformable cells, the equilibrium shifts toward the wall because
//! deformability reduces the wall-repulsion component of lift.  We model this
//! via a deformability correction factor `(1 − DI)` applied to the wall-
//! repulsion term (Hur et al. 2011, *Lab Chip* 11, 912–920):
//!
//! ```text
//! C_L(x̃) = C_wall(x̃) · (1 − DI) − C_center(x̃)
//! ```
//!
//! where `C_wall` is the wall-repulsion coefficient and `C_center` is the
//! Saffman-type shear-gradient lift toward the wall.
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
//! We solve this numerically using bisection on `[0, 0.95]` (by symmetry,
//! only the half-channel need be considered).
//!
//! # References
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.
//! - Gossett, D. R. & Di Carlo, D. (2009). Particle focusing mechanisms in
//!   curving confined flows. *Anal. Chem.*, 81, 8459–8465.
//! - Hur, S. C., Henderson-MacLennan, N. K., McCabe, E. R. B. & Di Carlo, D.
//!   (2011). Deformability-based cell classification and enrichment using
//!   inertial microfluidics. *Lab Chip*, 11, 912–920.

use crate::cell_separation::properties::CellProperties;
use serde::{Deserialize, Serialize};

// ── Lift coefficient model ────────────────────────────────────────────────────

/// Dimensionless wall-repulsion lift coefficient as a function of normalised
/// lateral position `x̃ ∈ [0, 1]` (0 = center, 1 = wall).
///
/// Fitted to the numerical data of Di Carlo et al. (2009), Fig. 2:
/// `C_wall(x̃) = 0.5 · x̃² · (1 − x̃)⁻¹`
///
/// This is a simplified but physically motivated fit: the wall repulsion
/// diverges as the particle approaches the wall (`x̃ → 1`) and vanishes at
/// the center (`x̃ = 0`).
#[inline]
fn c_wall(x_tilde: f64) -> f64 {
    let x = x_tilde.clamp(0.0, 0.95); // avoid singularity at wall
    0.5 * x * x / (1.0 - x).max(0.05)
}

/// Dimensionless shear-gradient (Saffman) lift coefficient as a function of
/// normalised lateral position `x̃ ∈ [0, 1]`.
///
/// Fitted to Di Carlo et al. (2009), Fig. 2:
/// `C_center(x̃) = 0.3 · (1 − x̃²)`
///
/// This term drives particles toward the wall (negative lift relative to center)
/// and is maximum at the center, vanishing at the wall.
#[inline]
fn c_center(x_tilde: f64) -> f64 {
    let x = x_tilde.clamp(0.0, 1.0);
    0.3 * (1.0 - x * x)
}

/// Net dimensionless lift coefficient `C_L(x̃, DI)`.
///
/// Positive → force toward center (away from wall).
/// Negative → force toward wall.
///
/// Deformability index `DI ∈ [0, 1]` reduces the wall-repulsion component,
/// shifting the equilibrium toward the wall for deformable cells (Hur 2011).
#[inline]
fn c_lift(x_tilde: f64, deformability_index: f64) -> f64 {
    let di = deformability_index.clamp(0.0, 1.0);
    c_wall(x_tilde) * (1.0 - di) - c_center(x_tilde)
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
pub fn dean_drag_force_n(
    dynamic_viscosity_pa_s: f64,
    de: f64,
    cell_diameter_m: f64,
) -> f64 {
    5.4e-4 * std::f64::consts::PI * dynamic_viscosity_pa_s * de.powf(1.63) * cell_diameter_m
}

// ── Inertial lift force ───────────────────────────────────────────────────────

/// Inertial lift force [N] on a cell at normalised lateral position `x̃`.
///
/// `F_L = C_L(x̃, DI) · ρ · U² · a⁴ / H²`
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
    let cl = c_lift(x_tilde, cell.deformability_index);
    let a = cell.diameter_m;
    let h = channel_height_m;
    cl * fluid_density_kg_m3 * mean_velocity_m_s * mean_velocity_m_s * a.powi(4) / (h * h)
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
/// Solves `F_L(x̃) − F_D = 0` using bisection on `x̃ ∈ [0, 0.95]`.
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

    // Bisection on [0, 0.95] (avoid wall singularity)
    let mut lo = 0.0_f64;
    let mut hi = 0.95_f64;

    // Check if there is a sign change (equilibrium exists in interval)
    let f_lo = net_force(lo);
    let f_hi = net_force(hi);

    // If both ends have the same sign, the equilibrium is at the boundary
    let x_eq = if f_lo * f_hi > 0.0 {
        // No zero crossing: equilibrium is at the end with smaller |F|
        if f_lo.abs() < f_hi.abs() { lo } else { hi }
    } else {
        // Bisection: 60 iterations → precision < 0.95 / 2^60 ≈ 8×10⁻¹⁹
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
