//! Cell-cell interaction correction for lateral equilibrium positions.
//!
//! Wraps [`cfd_core::physics::cell_interaction`] to apply the empirical
//! cell-free layer (CFL) model to the 1D inertial equilibrium solver.
//!
//! # When to use
//!
//! Use [`enhanced_lateral_equilibrium`] instead of
//! [`crate::physics::cell_separation::margination::lateral_equilibrium`] when:
//! - Feed hematocrit `≥ 0.01` (≥ 1% RBC volume fraction)
//! - Channel `D_h ≤ 300 µm` (CFL model valid range)
//!
//! At very dilute HCT (< 1%), the cell-cell interaction correction is
//! negligible (Γ ≈ 1.0) and the plain inertial model is sufficient.
//!
//! # Physical basis
//!
//! At finite hematocrit, RBCs form a dense core in the channel center,
//! creating a cell-free plasma layer near the walls.  WBCs larger than the
//! CFL thickness are excluded from the CFL by the RBC core boundary →
//! hydrodynamic collisions push them closer to the wall, enhancing separation
//! efficiency.
//!
//! The correction applies the enhancement factor Γ:
//!
//! ```text
//! x̃_eff = 1 − (1 − x̃_inertial) / Γ
//! ```
//!
//! See `cfd_core::physics::cell_interaction` for the physics and theorems.

use cfd_core::physics::cell_interaction as ci;
use cfd_core::error::{Error, Result};

use super::{margination, margination::EquilibriumResult, properties::CellProperties};

/// Compute the lateral equilibrium position of a WBC corrected for
/// RBC-WBC cell-cell interaction at the given feed hematocrit.
///
/// Applies the empirical CFL model (Fedosov 2012) on top of the single-particle
/// inertial equilibrium solver ([`margination::lateral_equilibrium`]).
///
/// # Returns
/// - `Some(result)` — equilibrium position with CFL correction applied.
/// - `None` — cell does not exhibit measurable inertial focusing (κ ≤ 0.07),
///   even after the CCI correction (since Γ only shifts x̃ further toward the
///   wall and cannot create focusing where none exists).
///
/// # Arguments
/// - `wbc` — WBC physical properties (diameter, deformability index, density)
/// - `fluid_density` — fluid density [kg/m³]
/// - `viscosity` — dynamic viscosity [Pa·s]
/// - `velocity` — mean channel velocity [m/s]
/// - `width` — channel width (wider dimension) [m]
/// - `height` — channel height (shorter dimension) [m]
/// - `bend_radius` — bend radius [m], or `None` for straight channel
/// - `hematocrit` — volumetric RBC fraction at device inlet (0.0–0.45)
#[must_use]
pub fn enhanced_lateral_equilibrium(
    wbc: &CellProperties,
    fluid_density: f64,
    viscosity: f64,
    velocity: f64,
    width: f64,
    height: f64,
    bend_radius: Option<f64>,
    hematocrit: f64,
) -> Option<EquilibriumResult> {
    let result = checked_enhanced_lateral_equilibrium(
        wbc,
        fluid_density,
        viscosity,
        velocity,
        width,
        height,
        bend_radius,
        hematocrit,
    )
    .ok()?;

    result.will_focus.then_some(result)
}

/// Checked WBC equilibrium including the hematocrit-dependent cell-interaction correction.
pub fn checked_enhanced_lateral_equilibrium(
    wbc: &CellProperties,
    fluid_density: f64,
    viscosity: f64,
    velocity: f64,
    width: f64,
    height: f64,
    bend_radius: Option<f64>,
    hematocrit: f64,
) -> Result<EquilibriumResult> {
    if !hematocrit.is_finite() || !(0.0..=1.0).contains(&hematocrit) {
        return Err(Error::InvalidConfiguration(
            "Cell-interaction hematocrit must lie in [0, 1]".to_string(),
        ));
    }

    // 1. Pure inertial equilibrium (existing single-particle model)
    let mut result = margination::checked_lateral_equilibrium(
        wbc,
        fluid_density,
        viscosity,
        velocity,
        width,
        height,
        bend_radius,
    )?;

    if !result.will_focus {
        return Ok(result);
    }

    // Skip correction at negligible HCT (< 0.5%) — Γ ≈ 1 and correction is
    // below floating-point noise in the inertial model.
    if hematocrit < 0.005 {
        return Ok(result);
    }

    // 2. Hydraulic diameter D_h = 2wh / (w + h)
    let dh = 2.0 * width * height / (width + height);

    // 3. Apply CFL + margination enhancement factor
    let x_tilde_corrected =
        ci::apply_cell_interaction(result.x_tilde_eq, wbc.diameter_m, hematocrit, dh);

    // 4. Update result with corrected position
    result.x_tilde_eq = x_tilde_corrected;
    result.lateral_position_m = x_tilde_corrected * (height * 0.5);

    Ok(result)
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::physics::cell_separation::properties::CellProperties;

    #[test]
    fn checked_enhanced_rejects_invalid_hematocrit() {
        let wbc = CellProperties::white_blood_cell();
        let err = checked_enhanced_lateral_equilibrium(
            &wbc,
            1060.0,
            3.5e-3,
            0.05,
            400e-6,
            80e-6,
            None,
            1.2,
        )
        .expect_err("checked enhanced equilibrium must reject hematocrit above unity");
        assert!(err.to_string().contains("hematocrit"));
    }

    /// At zero hematocrit, enhanced and plain equilibrium should agree
    #[test]
    fn enhanced_matches_plain_at_zero_hct() {
        let wbc = CellProperties::white_blood_cell();
        let w = 400e-6_f64;
        let h = 80e-6_f64;
        let v = 0.05_f64; // 5 cm/s
        let rho = 1060.0_f64;
        let mu = 3.5e-3_f64;

        let plain = margination::lateral_equilibrium(&wbc, rho, mu, v, w, h, None);
        let enhanced = enhanced_lateral_equilibrium(&wbc, rho, mu, v, w, h, None, 0.0);

        match (plain, enhanced) {
            (Some(p), Some(e)) => {
                assert!(
                    (p.x_tilde_eq - e.x_tilde_eq).abs() < 1e-9,
                    "At HCT=0 plain={:.6} enhanced={:.6}",
                    p.x_tilde_eq,
                    e.x_tilde_eq
                );
            }
            (None, None) => {} // Both agree cell doesn't focus
            _ => panic!("plain and enhanced should agree at HCT=0"),
        }
    }

    /// At physiological HCT, WBC should be pushed further toward wall
    #[test]
    fn enhanced_shifts_toward_wall_at_physiological_hct() {
        let wbc = CellProperties::white_blood_cell();
        let w = 200e-6_f64;
        let h = 80e-6_f64;
        let v = 0.10_f64;
        let rho = 1060.0_f64;
        let mu = 3.5e-3_f64;

        let plain = margination::lateral_equilibrium(&wbc, rho, mu, v, w, h, None);
        let enhanced = enhanced_lateral_equilibrium(&wbc, rho, mu, v, w, h, None, 0.40);

        if let (Some(p), Some(e)) = (plain, enhanced) {
            assert!(
                e.x_tilde_eq >= p.x_tilde_eq,
                "Enhanced should be ≥ plain at HCT=40%: plain={:.4} enhanced={:.4}",
                p.x_tilde_eq,
                e.x_tilde_eq
            );
        }
        // If neither focuses (κ ≤ 0.07) both are None — test is vacuously satisfied
    }

    /// lateral_position_m should be consistent with x_tilde_eq × (height/2)
    #[test]
    fn lateral_position_consistent_with_x_tilde() {
        let wbc = CellProperties::white_blood_cell();
        let w = 300e-6_f64;
        let h = 100e-6_f64;
        let v = 0.08_f64;
        let rho = 1060.0_f64;
        let mu = 3.5e-3_f64;

        if let Some(result) = enhanced_lateral_equilibrium(&wbc, rho, mu, v, w, h, None, 0.20) {
            let expected_pos = result.x_tilde_eq * (h * 0.5);
            assert!(
                (result.lateral_position_m - expected_pos).abs() < 1e-12,
                "lateral_position_m={:.3e} expected={:.3e}",
                result.lateral_position_m,
                expected_pos
            );
        }
    }

    #[test]
    fn checked_enhanced_matches_legacy_nominal_case() {
        let wbc = CellProperties::white_blood_cell();
        let legacy = enhanced_lateral_equilibrium(&wbc, 1060.0, 3.5e-3, 0.10, 200e-6, 80e-6, None, 0.20)
            .expect("legacy enhanced equilibrium should succeed");
        let checked = checked_enhanced_lateral_equilibrium(
            &wbc,
            1060.0,
            3.5e-3,
            0.10,
            200e-6,
            80e-6,
            None,
            0.20,
        )
        .expect("checked enhanced equilibrium should succeed");

        assert!((legacy.x_tilde_eq - checked.x_tilde_eq).abs() < 1e-12);
        assert_eq!(legacy.will_focus, checked.will_focus);
    }
}
