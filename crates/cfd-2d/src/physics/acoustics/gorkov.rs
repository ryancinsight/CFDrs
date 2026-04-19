//! Primary Acoustic Radiation Force (ARF) and Gor'kov Potential.
//!
//! # Theorem — Gor'kov Potential (Gor'kov 1962)
//!
//! The acoustic radiation force on a small spherical particle ($a \ll \lambda$)
//! in an inviscid fluid is the negative gradient of the acoustic potential $U$:
//!
//! ```text
//! U = 2 \pi a^3 \rho_0 \left( \frac{\langle p^2 \rangle}{3 \rho_0^2 c_0^2} f_1 - \frac{\langle v^2 \rangle}{2} f_2 \right)
//! ```
//! And $\vec{F}_{ARF} = -\nabla U$.
//!
//! For a 1D standing wave $p(x) = p_a \cos(kx)$, the local time-averaged force simplifies to:
//!
//! ```text
//! F_x = 4 \pi \Phi a^3 k E_{ac} \sin(2kx)
//! ```
//!
//! where $E_{ac} = p_a^2 / (4 \rho_0 c_0^2)$ is the acoustic energy density, and
//! $\Phi$ is the acoustic contrast factor.
//!
//! **Proof of migration limits**:
//! By definition, $\vec{F}_{ARF} = -\nabla U$.
//! Particles with $\Phi > 0$ (e.g., solid cells in water, $f_1 > 0, f_2 > 0$) have minima
//! of $U$ at the pressure nodes ($\cos(kx)=0$). Thus, they migrate to the pressure nodes.
//! Particles with $\Phi < 0$ (e.g., lipid droplets or microbubbles) have minima at
//! pressure antinodes, migrating there.
//!
//! # References
//! - Gor'kov, L. P. (1962). On the forces acting on a small particle in an acoustical field
//!   in an ideal fluid. *Soviet Physics Doklady*, 6(9), 773-775.
//! - Bruus, H. (2012). Acoustofluidics 7: The acoustic radiation force on small particles.
//!   *Lab on a Chip*, 12(6), 1014-1021.

use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// Core parameters for calculating Acoustic Radiation Force fields
#[derive(Debug, Clone, Copy)]
pub struct GorkovPotential<T: RealField + Copy> {
    /// Fluid density $\rho_0$ [kg/m³]
    pub fluid_density: T,
    /// Fluid speed of sound $c_0$ [m/s]
    pub fluid_sound_speed: T,
    /// Particle density $\rho_p$ [kg/m³]
    pub particle_density: T,
    /// Particle speed of sound $c_p$ [m/s]
    pub particle_sound_speed: T,
    /// Particle radius $a$ [m]
    pub particle_radius: T,
}

impl<T: RealField + Copy + Float + FromPrimitive> GorkovPotential<T> {
    /// Instantiate a typical configuration for standard human Red Blood Cells in water/plasma.
    ///
    /// Values derived from Bruus (2012) standard physiological tables:
    /// - Plasma: $\rho_0 = 1000$ kg/m³, $c_0 = 1500$ m/s
    /// - RBC: $\rho_p = 1090$ kg/m³, $c_p = 1639$ m/s
    /// - Radius: $a = 2.5$ µm (spherical equivalent volume)
    #[must_use]
    pub fn typical_rbc() -> Self {
        Self {
            fluid_density: T::from_f64(1000.0).unwrap(),
            fluid_sound_speed: T::from_f64(1500.0).unwrap(),
            particle_density: T::from_f64(1090.0).unwrap(),
            particle_sound_speed: T::from_f64(1639.0).unwrap(),
            particle_radius: T::from_f64(2.5e-6).unwrap(),
        }
    }

    /// Compute the fluid compressibility $\kappa_0 = 1 / (\rho_0 c_0^2)$.
    #[inline]
    #[must_use]
    pub fn fluid_compressibility(&self) -> T {
        let one = T::one();
        one / (self.fluid_density * self.fluid_sound_speed * self.fluid_sound_speed)
    }

    /// Compute the particle compressibility $\kappa_p = 1 / (\rho_p c_p^2)$.
    #[inline]
    #[must_use]
    pub fn particle_compressibility(&self) -> T {
        let one = T::one();
        one / (self.particle_density * self.particle_sound_speed * self.particle_sound_speed)
    }

    /// Compute the monopole scattering coefficient $f_1$ (compressibility contrast).
    ///
    /// ```text
    /// f_1 = 1 - \frac{\kappa_p}{\kappa_0}
    /// ```
    #[inline]
    #[must_use]
    pub fn f1_monopole(&self) -> T {
        let one = T::one();
        one - (self.particle_compressibility() / self.fluid_compressibility())
    }

    /// Compute the dipole scattering coefficient $f_2$ (density contrast).
    ///
    /// ```text
    /// f_2 = \frac{2(\rho_p - \rho_0)}{2\rho_p + \rho_0}
    /// ```
    #[inline]
    #[must_use]
    pub fn f2_dipole(&self) -> T {
        let two = T::from_f64(2.0).unwrap();
        let num = two * (self.particle_density - self.fluid_density);
        let den = two * self.particle_density + self.fluid_density;
        num / den
    }

    /// Compute the 1D acoustic contrast factor $\Phi$.
    ///
    /// ```text
    /// \Phi = \frac{f_1}{3} + \frac{f_2}{2}
    /// ```
    #[inline]
    #[must_use]
    pub fn contrast_factor(&self) -> T {
        let three = T::from_f64(3.0).unwrap();
        let two = T::from_f64(2.0).unwrap();
        (self.f1_monopole() / three) + (self.f2_dipole() / two)
    }

    /// Compute the Primary Acoustic Radiation Force (ARF) in a 1D standing wave.
    ///
    /// Evaluates the analytical gradient of the potential well:
    /// $F_x(x) = 4\pi \Phi a^3 k E_{ac} \sin(2kx)$
    ///
    /// # Arguments
    /// * `pressure_amplitude` - $p_a$ Peak acoustic pressure [Pa]
    /// * `frequency` - $f$ driving frequency [Hz]
    /// * `x` - Spatial coordinate along wave axis [m]
    #[inline]
    #[must_use]
    pub fn standing_wave_force_1d(&self, pressure_amplitude: T, frequency: T, x: T) -> T {
        let pi = T::from_f64(std::f64::consts::PI).unwrap();
        let two = T::from_f64(2.0).unwrap();
        let four = T::from_f64(4.0).unwrap();

        let wave_number = (two * pi * frequency) / self.fluid_sound_speed;
        let e_ac = (pressure_amplitude * pressure_amplitude)
            / (four * self.fluid_density * self.fluid_sound_speed * self.fluid_sound_speed);

        let a3 = self.particle_radius * self.particle_radius * self.particle_radius;
        let phi = self.contrast_factor();

        let sin_2kx = Float::sin(two * wave_number * x);

        four * pi * phi * a3 * wave_number * e_ac * sin_2kx
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    /// Theorem: Solid particles (like cells) in water move to pressure nodes,
    /// meaning $\Phi > 0$. Microbubbles ($\rho_p \ll \rho_0, \kappa_p \gg \kappa_0$)
    /// move to antinodes, meaning $\Phi < 0$ and strongly negative.
    #[test]
    fn thermodynamic_contrast_factor_signs() {
        let rbc = GorkovPotential::<f64>::typical_rbc();
        assert!(
            rbc.contrast_factor() > 0.0,
            "RBC contrast factor MUST be positive to move to nodes"
        );

        let microbubble = GorkovPotential::<f64> {
            fluid_density: 1000.0,
            fluid_sound_speed: 1500.0,
            particle_density: 1.2, // air
            particle_sound_speed: 343.0,
            particle_radius: 2.0e-6,
        };
        assert!(
            microbubble.contrast_factor() < 0.0,
            "Microbubble contrast factor MUST be negative to move to antinodes"
        );
        // Air in water is exceptionally responsive
        assert!(
            microbubble.contrast_factor() < -1000.0,
            "Microbubble contrast factor should be massively negative"
        );
    }

    /// Force at pressure node x=0: F = 4πΦa³kE_ac·sin(0) = 0 exactly.
    #[test]
    fn force_at_pressure_node_is_zero() {
        let rbc = GorkovPotential::<f64>::typical_rbc();
        let f = rbc.standing_wave_force_1d(1e6, 1e6, 0.0);
        assert_relative_eq!(f, 0.0, epsilon = 1e-30);
    }

    /// Analytical f1 and f2 for RBC in plasma.
    /// f1 = 1 − κ_p/κ_0 = 1 − (ρ₀c₀²)/(ρ_pc_p²)
    /// f2 = 2(ρ_p − ρ₀)/(2ρ_p + ρ₀)
    #[test]
    fn f1_f2_analytical_values() {
        let rbc = GorkovPotential::<f64>::typical_rbc();
        let kappa_0 = 1.0 / (1000.0 * 1500.0_f64.powi(2));
        let kappa_p = 1.0 / (1090.0 * 1639.0_f64.powi(2));
        let f1_expected = 1.0 - kappa_p / kappa_0;
        let f2_expected = 2.0 * (1090.0 - 1000.0) / (2.0 * 1090.0 + 1000.0);

        assert_relative_eq!(rbc.f1_monopole(), f1_expected, epsilon = 1e-10);
        assert_relative_eq!(rbc.f2_dipole(), f2_expected, epsilon = 1e-10);
    }

    /// Force sign at x = λ/8 for Φ > 0 particle: force pushes toward node (x=0).
    /// sin(2k·λ/8) = sin(π/2) = 1 > 0, so F > 0 for Φ > 0.
    /// At x = λ/8, the nearest node is x=0, and F > 0 means force is in the
    /// +x direction away from x=0; but that is the antinode-to-node direction
    /// since x=λ/4 is the antinode. Particles at x < λ/4 feel force toward
    /// the nearest node (x=0 or x=λ/2), confirming migration.
    #[test]
    fn force_sign_at_quarter_wavelength() {
        let rbc = GorkovPotential::<f64>::typical_rbc();
        let freq = 1e6_f64;
        let c0 = 1500.0_f64;
        let lambda = c0 / freq;
        // At x = λ/8, sin(2kx) = sin(2·(2π/λ)·(λ/8)) = sin(π/2) = 1
        let x = lambda / 8.0;
        let f = rbc.standing_wave_force_1d(1e6, freq, x);
        // Φ > 0 → F_x has same sign as sin(2kx) = +1 → F > 0
        assert!(f > 0.0, "Force at λ/8 must be positive for Φ>0 particle");
    }
}
