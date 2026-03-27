//! Cell-free layer (CFL) width correlations for microchannels.
//!
//! # Theorem — Cell-Free Layer Formation (Poiseuille flow + RBC margination)
//!
//! In microchannels with diameter $D < 300\,\mu m$, red blood cells migrate
//! toward the tube axis under the action of the Magnus-like lift force,
//! creating a plasma-only marginal layer of thickness $\delta$ near the wall.
//!
//! The CFL width governs:
//! 1. Apparent viscosity via the Fahraeus-Lindqvist effect (two-layer model)
//! 2. Nutrient/gas transport to the vessel wall
//! 3. Hematocrit partitioning at bifurcations (Zweifach-Fung)
//!
//! ## Sharan-Popel (2001) CFL width model
//!
//! ```text
//! δ/R = 1 − √(H_T / H_F)
//! ```
//!
//! where $H_T$ is the tube hematocrit and $H_F$ is the feed (discharge)
//! hematocrit. Since $H_T < H_F$ in microtubes (Fahraeus effect), the CFL
//! width is always positive.
//!
//! ## Fedosov (2010) computational correlation
//!
//! From dissipative particle dynamics (DPD) simulations:
//!
//! ```text
//! δ/R = 0.29 · (1 − H_t)^{0.84}       for D ∈ [20, 100] µm
//! ```
//!
//! # Two-layer viscosity model
//!
//! The apparent viscosity from a core-annulus (CFL) structure is:
//!
//! ```text
//! µ_app = µ_plasma / [1 − (1 − δ/R)^4 · (1 − µ_plasma/µ_core)]
//! ```
//!
//! This is exact for concentric Poiseuille flow with two immiscible fluids
//! of different viscosities.
//!
//! # References
//! - Sharan, M. & Popel, A.S. (2001). A two-phase model for flow of blood
//!   in narrow tubes with increased effective viscosity near the wall.
//!   *Biorheology* 38:415-428.
//! - Fedosov, D.A. et al. (2010). A multiscale red blood cell model with
//!   accurate mechanics. *Biophys. J.* 98:2215-2225.
//! - Balogh, P. & Bagchi, P. (2017). A computational approach to modeling
//!   cellular-scale blood flow in complex geometry. *J. Comput. Phys.* 334:280-307.

/// Fedosov (2010) DPD-based CFL width correlation.
///
/// ```text
/// δ/R = 0.29 · (1 − H_t)^{0.84}
/// ```
///
/// Valid for tube diameters $D \in [20, 100]\,\mu m$ and hematocrit $H_t \in [0.1, 0.5]$.
///
/// # Arguments
/// * `diameter_um` — Channel diameter [µm]
/// * `hematocrit` — Feed hematocrit [0, 1]
///
/// # Returns
/// CFL width [µm]
#[inline]
#[must_use]
pub fn cfl_width_fedosov(diameter_um: f64, hematocrit: f64) -> f64 {
    let ht = hematocrit.clamp(0.0, 1.0);
    let r = diameter_um.max(1.0) / 2.0;
    let delta_over_r = 0.29 * (1.0 - ht).powf(0.84);
    r * delta_over_r
}

/// Sharan-Popel (2001) CFL width from the Fahraeus ratio.
///
/// ```text
/// δ/R = 1 − √(H_T / H_F)
/// ```
///
/// The tube hematocrit $H_T$ is computed from the Pries (1990) Fahraeus
/// ratio model (see [`super::fahraeus_effect::tube_hematocrit_ratio`]).
///
/// # Arguments
/// * `diameter_um` — Channel diameter [µm]
/// * `hematocrit` — Feed hematocrit [0, 1]
///
/// # Returns
/// CFL width [µm]
#[must_use]
pub fn cfl_width_sharan_popel(diameter_um: f64, hematocrit: f64) -> f64 {
    let ht = hematocrit.clamp(0.0, 1.0);
    if ht < 1e-15 {
        // Pure plasma: CFL = full radius
        return diameter_um.max(1.0) / 2.0;
    }
    let r = diameter_um.max(1.0) / 2.0;
    let ht_ratio = super::fahraeus_effect::tube_hematocrit_ratio(diameter_um);
    // H_T / H_F = ht_ratio, so sqrt(ht_ratio) gives the core radius fraction
    let delta_over_r = (1.0 - ht_ratio.sqrt()).max(0.0);
    r * delta_over_r
}

/// Two-layer (core + CFL) apparent viscosity model.
///
/// # Theorem — Two-Layer Viscosity (Sharan & Popel 2001)
///
/// For concentric Poiseuille flow with an inner core of viscosity $\mu_c$
/// and an outer annulus (CFL) of viscosity $\mu_p$, with CFL thickness $\delta$
/// and tube radius $R$:
///
/// ```text
/// Q = (π R⁴ ΔP) / (8 µ_app L)
/// ```
///
/// where the apparent viscosity is:
///
/// ```text
/// µ_app = µ_p / [1 − (1 − δ/R)⁴ · (1 − µ_p/µ_c)]
/// ```
///
/// **Proof**: Integrate the Navier-Stokes equations in cylindrical coordinates
/// for two concentric regions with matched velocity and shear at the interface
/// $r = R - \delta$. The total flow rate is the sum of the core and annular
/// contributions, yielding an effective viscosity that reduces to $\mu_p$
/// when $\delta = R$ (pure plasma) and to $\mu_c$ when $\delta = 0$ (no CFL).
///
/// # Arguments
/// * `diameter_um` — Channel diameter [µm]
/// * `hematocrit` — Feed hematocrit [0, 1]
/// * `mu_plasma_pa_s` — Plasma viscosity [Pa·s] (typically 0.0012)
/// * `mu_core_pa_s` — Core (RBC-rich) viscosity [Pa·s] (typically 3.5-5× plasma)
///
/// # Returns
/// Apparent viscosity [Pa·s]
#[must_use]
pub fn two_layer_viscosity(
    diameter_um: f64,
    hematocrit: f64,
    mu_plasma_pa_s: f64,
    mu_core_pa_s: f64,
) -> f64 {
    if hematocrit < 1e-15 {
        return mu_plasma_pa_s;
    }

    let r = diameter_um.max(1.0) / 2.0;
    let delta = cfl_width_fedosov(diameter_um, hematocrit);
    let delta_over_r = (delta / r).clamp(0.0, 1.0);
    let core_ratio = 1.0 - delta_over_r; // (R-δ)/R

    let denom = 1.0 - core_ratio.powi(4) * (1.0 - mu_plasma_pa_s / mu_core_pa_s.max(1e-15));
    if denom.abs() < 1e-15 {
        return mu_core_pa_s;
    }

    mu_plasma_pa_s / denom
}

#[cfg(test)]
mod tests {
    use super::*;

    const MU_PLASMA: f64 = 0.0012;

    /// CFL width increases with diameter (larger tubes → larger absolute CFL).
    #[test]
    fn cfl_width_increases_with_diameter() {
        let d_small = cfl_width_fedosov(20.0, 0.45);
        let d_medium = cfl_width_fedosov(50.0, 0.45);
        let d_large = cfl_width_fedosov(100.0, 0.45);
        assert!(
            d_small < d_medium && d_medium < d_large,
            "CFL: {d_small:.3} < {d_medium:.3} < {d_large:.3}"
        );
    }

    /// CFL width decreases with hematocrit (more RBCs → thinner CFL).
    #[test]
    fn cfl_width_decreases_with_hematocrit() {
        let d = 50.0;
        let cfl_low = cfl_width_fedosov(d, 0.20);
        let cfl_mid = cfl_width_fedosov(d, 0.35);
        let cfl_high = cfl_width_fedosov(d, 0.50);
        assert!(
            cfl_low > cfl_mid && cfl_mid > cfl_high,
            "CFL: {cfl_low:.3} > {cfl_mid:.3} > {cfl_high:.3}"
        );
    }

    /// CFL width is always positive and less than the radius.
    #[test]
    fn cfl_width_bounded() {
        for d in [10.0, 30.0, 50.0, 100.0, 200.0] {
            for ht in [0.1, 0.3, 0.45, 0.6] {
                let cfl = cfl_width_fedosov(d, ht);
                assert!(cfl > 0.0, "CFL must be positive: d={d}, ht={ht}");
                assert!(cfl < d / 2.0, "CFL must be < radius: cfl={cfl:.3}, R={}", d / 2.0);
            }
        }
    }

    /// Two-layer viscosity: pure plasma (Ht=0) returns plasma viscosity.
    #[test]
    fn two_layer_pure_plasma() {
        let mu = two_layer_viscosity(50.0, 0.0, MU_PLASMA, 0.005);
        assert!(
            (mu - MU_PLASMA).abs() < 1e-10,
            "Pure plasma viscosity: {mu:.6} should be {MU_PLASMA:.6}"
        );
    }

    /// Two-layer viscosity increases with hematocrit.
    #[test]
    fn two_layer_viscosity_increases_with_hematocrit() {
        let mu_core = 0.005;
        let mu_low = two_layer_viscosity(50.0, 0.20, MU_PLASMA, mu_core);
        let mu_mid = two_layer_viscosity(50.0, 0.35, MU_PLASMA, mu_core);
        let mu_high = two_layer_viscosity(50.0, 0.50, MU_PLASMA, mu_core);
        assert!(
            mu_low < mu_mid && mu_mid < mu_high,
            "Viscosity: {mu_low:.6} < {mu_mid:.6} < {mu_high:.6}"
        );
    }

    /// Two-layer viscosity is between plasma and core viscosity.
    #[test]
    fn two_layer_viscosity_between_bounds() {
        let mu_core = 0.005;
        let mu = two_layer_viscosity(50.0, 0.45, MU_PLASMA, mu_core);
        assert!(
            mu > MU_PLASMA && mu < mu_core,
            "Apparent viscosity {mu:.6} should be between {MU_PLASMA} and {mu_core}"
        );
    }
}
