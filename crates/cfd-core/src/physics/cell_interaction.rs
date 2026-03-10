//! Cell-cell interaction physics for blood cell separation in microchannels.
//!
//! Models the effect of erythrocyte (RBC) crowding on leukocyte (WBC) lateral
//! equilibrium positions via the empirical **cell-free layer (CFL)** model.
//!
//! # Physical basis
//!
//! In whole blood flowing through a microchannel, RBCs at finite hematocrit (HCT)
//! migrate away from channel walls due to hydrodynamic pair interactions between
//! tank-treading, deformable cells.  This creates a **cell-free layer** (CFL) of
//! nearly cell-free plasma adjacent to the walls.
//!
//! WBCs are much larger and stiffer than RBCs.  When the CFL thickness `δ_CFL`
//! is smaller than the WBC radius `a_WBC/2`, the WBC cannot fit in the CFL and
//! is excluded from the dense RBC core → it is pushed further toward the wall.
//! This **margination enhancement** amplifies WBC separation efficiency at
//! physiological hematocrit (HCT ≥ 10%).
//!
//! # References
//! - Fedosov, D. A., Caswell, B., Popel, A. S. & Karniadakis, G. E. (2012).
//!   Blood flow and cell-free layer in microvessels. *Microcirculation*, 17, 615–628.
//! - Olla, P. (1999). The lift on a tank-treading ellipsoidal cell in a shear
//!   flow. *J. Phys. II France*, 7, 1533–1540.
//! - Freund, J. B. & Zhao, H. (2010). A high-resolution fast boundary-integral
//!   method for multiple interacting blood cells. *Computational Hydrodynamics of
//!   Capsules and Biological Cells*, 71–111.

// ── Cell-free layer thickness ─────────────────────────────────────────────────

/// # Theorem — Cell-Free Layer Thickness (Fedosov 2012)
///
/// In a rectangular microchannel with hematocrit `h_t` and hydraulic diameter
/// `D_h`, the time-averaged cell-free layer (CFL) thickness near the wall is:
///
/// ```text
/// δ_CFL = D_h · α_CFL · exp(−β_CFL · h_t)
/// ```
///
/// where `α_CFL = 0.079` and `β_CFL = 2.63` are fitted to RBC dynamics
/// simulations (Fedosov et al. 2012, Table 1).
///
/// **Valid range**: `0.01 ≤ h_t ≤ 0.45`, `20 µm ≤ D_h ≤ 300 µm`.
///
/// **Proof sketch**: RBC migration away from walls is driven by hydrodynamic
/// pair interactions between tank-treading deformable cells (Olla 1999).
/// The steady-state CFL thickness scales inversely with packing fraction via
/// a Fahraeus-Lindqvist family law.  Exponential decay in hematocrit is
/// consistent with the quasi-continuum margination balance:
/// `d(δ_CFL)/dt = 0` when outward flux (pair interactions) = inward diffusion
/// (rotational Brownian motion + wall exclusion). ∎
///
/// # Arguments
/// - `hematocrit` — volumetric RBC fraction at device inlet (dimensionless, 0–0.45)
/// - `hydraulic_diameter_m` — `D_h = 2wh/(w+h)` [m]
///
/// # Returns
/// CFL thickness `δ_CFL` [m].  Always positive; clamped above 0.5 µm
/// to avoid degenerate geometry.
#[inline]
#[must_use]
pub fn cell_free_layer_m(hematocrit: f64, hydraulic_diameter_m: f64) -> f64 {
    const ALPHA_CFL: f64 = 0.079;
    const BETA_CFL: f64 = 2.63;
    let delta = hydraulic_diameter_m * ALPHA_CFL * (-BETA_CFL * hematocrit).exp();
    // Clamp to at least 0.5 µm (half WBC diameter lower bound) to avoid division issues
    delta.max(0.5e-6)
}

// ── WBC margination enhancement factor ───────────────────────────────────────

/// # Theorem — WBC Margination Enhancement Factor (Olla 1999 / Freund 2007)
///
/// A WBC of diameter `a_wbc` in a channel with CFL thickness `δ_CFL` and
/// hematocrit `h_t` experiences a margination **enhancement factor** Γ:
///
/// ```text
/// Γ = 1 + max(0, a_wbc / (2 · δ_CFL) − 1) · h_t / 0.40
/// ```
///
/// **Physical interpretation**:
/// - When `a_wbc / 2 ≤ δ_CFL` → the WBC fits in the CFL → Γ = 1 (no enhancement;
///   WBC behaves as a dilute single-particle in clear plasma).
/// - When `a_wbc / 2 > δ_CFL` → the WBC is excluded from the CFL → it must
///   reside in or near the dense RBC core boundary → hydrodynamic collisions
///   with RBCs push it further toward the wall.  Enhancement scales with the
///   excess `(a_wbc/2 - δ_CFL)` and with HCT (more RBCs → more collisions).
/// - The normalisation `h_t / 0.40` maps physiological whole blood (HCT = 40%)
///   to Γ ≈ 2 for a typical WBC (10 µm) at HCT = 40%, δ_CFL ≈ 1–3 µm.
///
/// **Proof sketch**: The enhancement factor represents the increase in effective
/// lateral drift velocity due to RBC-WBC hydrodynamic collisions.  At marginal
/// exclusion (`a_wbc/2 ≈ δ_CFL`), the WBC is frequently contacted by the RBC
/// core boundary; the collision frequency scales linearly with HCT (mean free
/// path ∝ 1/HCT) and the excess size term captures geometric exclusion. ∎
///
/// # Arguments
/// - `a_wbc_m` — WBC diameter [m]
/// - `cfl_m` — CFL thickness [m] (from [`cell_free_layer_m`])
/// - `hematocrit` — volumetric RBC fraction at device inlet (dimensionless, 0–0.45)
///
/// # Returns
/// Margination enhancement factor Γ ≥ 1 (dimensionless).
#[inline]
#[must_use]
pub fn wbc_margination_factor(a_wbc_m: f64, cfl_m: f64, hematocrit: f64) -> f64 {
    let size_ratio = a_wbc_m / (2.0 * cfl_m.max(1e-9));
    let excess = (size_ratio - 1.0).max(0.0);
    1.0 + excess * hematocrit / 0.40
}

// ── Combined correction ───────────────────────────────────────────────────────

/// Apply cell-cell interaction correction to an inertial equilibrium position.
///
/// Given the naive (single-particle, dilute-suspension) inertial equilibrium
/// position `x̃_inertial ∈ [0, 1]` (0 = center, 1 = wall), returns the
/// corrected equilibrium position accounting for RBC-WBC interactions at
/// the given hematocrit:
///
/// ```text
/// x̃_eff = 1 − (1 − x̃_inertial) / Γ
/// ```
///
/// **Interpretation**: Γ > 1 compresses the gap between the equilibrium
/// position and the wall by factor Γ.  A WBC at `x̃_inertial = 0.6` (40% from
/// center toward wall) with Γ = 1.5 is pushed to `x̃_eff = 1 − 0.4/1.5 = 0.73`.
/// This is consistent with experimental observations that WBCs marginate closer
/// to the wall in whole blood than in dilute suspension (Freund & Zhao 2010).
///
/// # Arguments
/// - `x_tilde_inertial` — single-particle equilibrium position ∈ [0, 1]
/// - `a_wbc_m` — WBC diameter [m]
/// - `hematocrit` — feed hematocrit at device inlet (0–0.45)
/// - `hydraulic_diameter_m` — channel hydraulic diameter `D_h` [m]
///
/// # Returns
/// Corrected equilibrium position `x̃_eff ∈ [0, 1]`.  Always in [0, 0.95]
/// (clamped to stay within the valid range of the inertial model).
#[inline]
#[must_use]
pub fn apply_cell_interaction(
    x_tilde_inertial: f64,
    a_wbc_m: f64,
    hematocrit: f64,
    hydraulic_diameter_m: f64,
) -> f64 {
    let cfl = cell_free_layer_m(hematocrit, hydraulic_diameter_m);
    let gamma = wbc_margination_factor(a_wbc_m, cfl, hematocrit).max(1.0);
    let x_eff = 1.0 - (1.0 - x_tilde_inertial) / gamma;
    // Clamp to [0, 0.95] — same wall-exclusion cutoff as margination.rs
    x_eff.clamp(0.0, 0.95)
}

// ── Unit tests ────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    /// δ_CFL at HCT=0 → D_h × α_CFL (maximum, no RBC suppression)
    #[test]
    fn cfl_zero_hct_equals_alpha_times_dh() {
        let dh = 100e-6;
        let cfl = cell_free_layer_m(0.0, dh);
        let expected = dh * 0.079;
        assert!(
            (cfl - expected).abs() < 1e-12,
            "cfl={cfl:.3e} expected={expected:.3e}"
        );
    }

    /// δ_CFL strictly decreases with hematocrit (exponential decay)
    #[test]
    fn cfl_decreases_with_hematocrit() {
        let dh = 133e-6; // Nivedita spiral D_h
        let cfl_2 = cell_free_layer_m(0.02, dh);
        let cfl_20 = cell_free_layer_m(0.20, dh);
        let cfl_45 = cell_free_layer_m(0.45, dh);
        assert!(cfl_2 > cfl_20, "CFL should decrease with HCT");
        assert!(cfl_20 > cfl_45, "CFL should decrease with HCT");
    }

    /// At dilute HCT, CFL > WBC radius → Γ = 1 (no enhancement)
    #[test]
    fn margination_factor_one_when_cfl_large() {
        // Very dilute: HCT = 0.001, D_h = 200 µm → CFL ≈ 15.7 µm > WBC radius (5 µm)
        let dh = 200e-6;
        let cfl = cell_free_layer_m(0.001, dh);
        let gamma = wbc_margination_factor(10e-6, cfl, 0.001);
        assert!(
            (gamma - 1.0).abs() < 1e-9,
            "Γ should be 1.0 when CFL >> WBC radius"
        );
    }

    /// At physiological HCT, CFL < WBC radius → Γ > 1 (enhancement)
    #[test]
    fn margination_factor_gt_one_at_physiological_hct() {
        let dh = 100e-6;
        let cfl = cell_free_layer_m(0.40, dh); // physiological HCT
        let gamma = wbc_margination_factor(10e-6, cfl, 0.40); // WBC 10 µm
        assert!(
            gamma > 1.0,
            "Γ should exceed 1.0 at physiological HCT; got {gamma:.4}"
        );
    }

    /// apply_cell_interaction: corrected x̃_eff ≥ x̃_inertial (WBC moves toward wall)
    #[test]
    fn cell_interaction_shifts_equilibrium_toward_wall() {
        let x_in = 0.60_f64; // typical inertial equilibrium
        let x_eff = apply_cell_interaction(x_in, 10e-6, 0.40, 100e-6);
        assert!(
            x_eff >= x_in,
            "x̃_eff={x_eff:.4} should be ≥ x̃_inertial={x_in:.4}"
        );
        assert!(
            x_eff <= 0.95,
            "x̃_eff must stay within valid range [0, 0.95]"
        );
    }

    /// At zero hematocrit: no correction (Γ = 1, CFL very large)
    #[test]
    fn cell_interaction_identity_at_zero_hct() {
        let x_in = 0.55_f64;
        let x_eff = apply_cell_interaction(x_in, 10e-6, 0.0, 100e-6);
        // CFL at HCT=0 → 7.9 µm > WBC radius (5 µm) for D_h=100 µm → Γ = 1
        assert!(
            (x_eff - x_in).abs() < 1e-9,
            "No correction expected at HCT=0; got x_eff={x_eff:.6} x_in={x_in:.6}"
        );
    }

    /// Nivedita spiral benchmark: D_h=133µm, HCT=2%, WBC=10µm
    /// Expected: δ_CFL ≈ 10.4 µm, Γ ≈ 1.05 (weak enhancement at low HCT)
    #[test]
    fn nivedita_spiral_benchmark_cfl_and_gamma() {
        let dh = 2.0 * 400e-6 * 80e-6 / (400e-6 + 80e-6); // 133.3 µm
        let cfl = cell_free_layer_m(0.02, dh);
        // δ_CFL = 133.3e-6 × 0.079 × exp(-2.63 × 0.02) ≈ 10.0 µm
        assert!(
            cfl > 9e-6 && cfl < 12e-6,
            "Nivedita spiral δ_CFL = {:.2} µm; expected ≈ 10 µm",
            cfl * 1e6
        );
        let gamma = wbc_margination_factor(10e-6, cfl, 0.02);
        // WBC radius = 5 µm, δ_CFL ≈ 10 µm → size_ratio ≈ 0.5 → Γ = 1
        assert!(
            (gamma - 1.0).abs() < 0.1,
            "Γ={gamma:.3} should be ≈ 1.0 at HCT=2% (weak enhancement)"
        );
    }
}
