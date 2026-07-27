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
/// * `diameter` — Channel diameter
/// * `hematocrit` — Feed hematocrit [0, 1]
///
/// # Returns
/// CFL width
#[inline]
pub fn cfl_width_fedosov(
    diameter: aequitas::systems::si::quantities::Length,
    hematocrit: f64,
) -> cfd_core::error::Result<aequitas::systems::si::quantities::Length> {
    let diameter_um = diameter.into_base() * 1.0e6;
    if !diameter_um.is_finite() || diameter_um <= 0.0 || !hematocrit.is_finite() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "Fedosov CFL inputs must be finite and diameter must be positive".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&hematocrit) {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "Fedosov CFL hematocrit must lie in [0, 1]".to_string(),
        ));
    }

    let r = diameter_um / 2.0;
    let delta_over_r = 0.29 * (1.0 - hematocrit).powf(0.84);
    Ok(aequitas::systems::si::quantities::Length::from_base(
        r * delta_over_r * 1.0e-6,
    ))
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
/// * `diameter` — Channel diameter
/// * `hematocrit` — Feed hematocrit [0, 1]
///
/// # Returns
/// CFL width
pub fn cfl_width_sharan_popel(
    diameter: aequitas::systems::si::quantities::Length,
    hematocrit: f64,
) -> cfd_core::error::Result<aequitas::systems::si::quantities::Length> {
    let diameter_um = diameter.into_base() * 1.0e6;
    if !diameter_um.is_finite() || diameter_um <= 0.0 || !hematocrit.is_finite() {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "Sharan-Popel CFL inputs must be finite and diameter must be positive".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&hematocrit) {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "Sharan-Popel CFL hematocrit must lie in [0, 1]".to_string(),
        ));
    }

    if hematocrit < 1e-15 {
        // Pure plasma: CFL = full radius
        return Ok(aequitas::systems::si::quantities::Length::from_base(
            diameter.into_base() / 2.0,
        ));
    }
    let r = diameter_um / 2.0;
    let ht_ratio = super::fahraeus_effect::tube_hematocrit_ratio(diameter)?;
    // H_T / H_F = ht_ratio, so sqrt(ht_ratio) gives the core radius fraction
    let delta_over_r = (1.0 - ht_ratio.sqrt()).max(0.0);
    Ok(aequitas::systems::si::quantities::Length::from_base(
        r * delta_over_r * 1.0e-6,
    ))
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
/// * `diameter` — Channel diameter
/// * `hematocrit` — Feed hematocrit [0, 1]
/// * `mu_plasma` — Plasma dynamic viscosity
/// * `mu_core` — Core (RBC-rich) dynamic viscosity
///
/// # Returns
/// Apparent dynamic viscosity
pub fn two_layer_viscosity(
    diameter: aequitas::systems::si::quantities::Length,
    hematocrit: f64,
    mu_plasma: aequitas::systems::si::quantities::DynamicViscosity,
    mu_core: aequitas::systems::si::quantities::DynamicViscosity,
) -> cfd_core::error::Result<aequitas::systems::si::quantities::DynamicViscosity> {
    let diameter_um = diameter.into_base() * 1.0e6;
    let mu_plasma_pa_s = mu_plasma.into_base();
    let mu_core_pa_s = mu_core.into_base();
    if !diameter_um.is_finite()
        || diameter_um <= 0.0
        || !hematocrit.is_finite()
        || !mu_plasma_pa_s.is_finite()
        || !mu_core_pa_s.is_finite()
    {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "Two-layer viscosity inputs must be finite and diameter/viscosities must be positive"
                .to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&hematocrit) {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "Two-layer viscosity hematocrit must lie in [0, 1]".to_string(),
        ));
    }
    if mu_plasma_pa_s <= 0.0 || mu_core_pa_s <= 0.0 {
        return Err(cfd_core::error::Error::InvalidConfiguration(
            "Two-layer viscosity values must be positive".to_string(),
        ));
    }

    if hematocrit < 1e-15 {
        return Ok(mu_plasma);
    }

    let r = diameter_um / 2.0;
    let delta = cfl_width_fedosov(diameter, hematocrit)?;
    let delta_over_r = (delta.into_base() * 1.0e6 / r).clamp(0.0, 1.0);
    let core_ratio = 1.0 - delta_over_r; // (R-δ)/R

    let denom = 1.0 - core_ratio.powi(4) * (1.0 - mu_plasma_pa_s / mu_core_pa_s);
    if denom.abs() < 1e-15 {
        return Ok(mu_core);
    }

    Ok(aequitas::systems::si::quantities::DynamicViscosity::from_base(mu_plasma_pa_s / denom))
}

#[cfg(test)]
mod tests {
    use super::*;
    use aequitas::systems::si::quantities::{DynamicViscosity, Length};

    const MU_PLASMA: DynamicViscosity = DynamicViscosity::from_base(0.0012);

    fn diameter_um(value: f64) -> Length {
        Length::from_base(value * 1.0e-6)
    }

    #[test]
    fn cfl_width_fedosov_rejects_nonpositive_diameter() {
        assert!(cfl_width_fedosov(diameter_um(0.0), 0.45).is_err());
    }

    /// CFL width increases with diameter (larger tubes → larger absolute CFL).
    #[test]
    fn cfl_width_increases_with_diameter() -> cfd_core::error::Result<()> {
        let d_small = cfl_width_fedosov(diameter_um(20.0), 0.45)?.into_base();
        let d_medium = cfl_width_fedosov(diameter_um(50.0), 0.45)?.into_base();
        let d_large = cfl_width_fedosov(diameter_um(100.0), 0.45)?.into_base();
        assert!(
            d_small < d_medium && d_medium < d_large,
            "CFL: {d_small:.3} < {d_medium:.3} < {d_large:.3}"
        );
        Ok(())
    }

    /// CFL width decreases with hematocrit (more RBCs → thinner CFL).
    #[test]
    fn cfl_width_decreases_with_hematocrit() -> cfd_core::error::Result<()> {
        let d = diameter_um(50.0);
        let cfl_low = cfl_width_fedosov(d, 0.20)?.into_base();
        let cfl_mid = cfl_width_fedosov(d, 0.35)?.into_base();
        let cfl_high = cfl_width_fedosov(d, 0.50)?.into_base();
        assert!(
            cfl_low > cfl_mid && cfl_mid > cfl_high,
            "CFL: {cfl_low:.3} > {cfl_mid:.3} > {cfl_high:.3}"
        );
        Ok(())
    }

    /// CFL width is always positive and less than the radius.
    #[test]
    fn cfl_width_bounded() -> cfd_core::error::Result<()> {
        for d in [10.0, 30.0, 50.0, 100.0, 200.0] {
            for ht in [0.1, 0.3, 0.45, 0.6] {
                let cfl = cfl_width_fedosov(diameter_um(d), ht)?.into_base();
                assert!(cfl > 0.0, "CFL must be positive: d={d}, ht={ht}");
                assert!(
                    cfl < d / 2.0,
                    "CFL must be < radius: cfl={cfl:.3}, R={}",
                    d / 2.0
                );
            }
        }
        Ok(())
    }

    /// Two-layer viscosity: pure plasma (Ht=0) returns plasma viscosity.
    #[test]
    fn two_layer_pure_plasma() -> cfd_core::error::Result<()> {
        let mu = two_layer_viscosity(
            diameter_um(50.0),
            0.0,
            MU_PLASMA,
            DynamicViscosity::from_base(0.005),
        )?;
        assert!(
            (mu.into_base() - MU_PLASMA.into_base()).abs() < 1e-10,
            "Pure plasma viscosity: {:.6} should be {:.6}",
            mu.into_base(),
            MU_PLASMA.into_base()
        );
        Ok(())
    }

    /// Two-layer viscosity increases with hematocrit.
    #[test]
    fn two_layer_viscosity_increases_with_hematocrit() -> cfd_core::error::Result<()> {
        let mu_core = DynamicViscosity::from_base(0.005);
        let mu_low = two_layer_viscosity(diameter_um(50.0), 0.20, MU_PLASMA, mu_core)?;
        let mu_mid = two_layer_viscosity(diameter_um(50.0), 0.35, MU_PLASMA, mu_core)?;
        let mu_high = two_layer_viscosity(diameter_um(50.0), 0.50, MU_PLASMA, mu_core)?;
        assert!(
            mu_low.into_base() < mu_mid.into_base() && mu_mid.into_base() < mu_high.into_base(),
            "Viscosity: {:.6} < {:.6} < {:.6}",
            mu_low.into_base(),
            mu_mid.into_base(),
            mu_high.into_base()
        );
        Ok(())
    }

    /// Two-layer viscosity is between plasma and core viscosity.
    #[test]
    fn two_layer_viscosity_between_bounds() -> cfd_core::error::Result<()> {
        let mu_core = DynamicViscosity::from_base(0.005);
        let mu = two_layer_viscosity(diameter_um(50.0), 0.45, MU_PLASMA, mu_core)?;
        assert!(
            mu.into_base() > MU_PLASMA.into_base() && mu.into_base() < mu_core.into_base(),
            "Apparent viscosity {:.6} should be between {:.6} and {:.6}",
            mu.into_base(),
            MU_PLASMA.into_base(),
            mu_core.into_base()
        );
        Ok(())
    }
}
