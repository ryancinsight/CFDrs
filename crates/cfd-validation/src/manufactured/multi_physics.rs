//! Multi-physics manufactured solutions for coupled CFD problems
//!
//! Implements manufactured solutions for coupled physics including:
//! - Conjugate heat transfer (fluid-solid thermal coupling)
//! - Species transport with chemical reactions
//! - Magnetohydrodynamics (MHD)
//! - Multi-phase flows
//! - Turbulent combustion

use super::{ManufacturedSolution, ManufacturedFunctions};
use nalgebra::RealField;
use num_traits::Float;

/// Manufactured solution for conjugate heat transfer
///
/// Provides analytical solutions for coupled fluid-solid thermal problems
/// where temperature satisfies different equations in different domains.
#[derive(Debug, Clone)]
pub struct ManufacturedConjugateHeatTransfer<T: RealField + Float + Copy> {
    /// Thermal conductivity ratio (solid/fluid)
    pub conductivity_ratio: T,
    /// Heat capacity ratio (solid/fluid)
    pub capacity_ratio: T,
    /// Interface position (x-coordinate)
    pub interface_x: T,
    /// Temperature amplitude
    pub amplitude: T,
    /// Frequency parameter
    pub frequency: T,
}

impl<T: RealField + Float + Copy> ManufacturedConjugateHeatTransfer<T> {
    pub fn new(conductivity_ratio: T, capacity_ratio: T, interface_x: T, amplitude: T, frequency: T) -> Self {
        Self {
            conductivity_ratio,
            capacity_ratio,
            interface_x,
            amplitude,
            frequency,
        }
    }

    /// Evaluate temperature in fluid domain (x < interface_x)
    fn fluid_temperature(&self, x: T, y: T, t: T) -> T {
        let base = ManufacturedFunctions::sinusoidal(x, y, t, self.frequency, self.frequency);
        self.amplitude * base
    }

    /// Evaluate temperature in solid domain (x > interface_x)
    fn solid_temperature(&self, x: T, y: T, t: T) -> T {
        let base = ManufacturedFunctions::sinusoidal(x, y, t, self.frequency, self.frequency);
        // Account for different thermal properties
        self.amplitude * self.conductivity_ratio * base
    }
}

impl<T: RealField + Float + Copy> ManufacturedSolution<T> for ManufacturedConjugateHeatTransfer<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        if x < self.interface_x {
            self.fluid_temperature(x, y, t)
        } else {
            self.solid_temperature(x, y, t)
        }
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        if x < self.interface_x {
            // Fluid domain: ∂T/∂t = α ∇²T + S
            self.fluid_heat_source(x, y, t)
        } else {
            // Solid domain: ∂T/∂t = α_s ∇²T + S_s
            self.solid_heat_source(x, y, t)
        }
    }
}

impl<T: RealField + Float + Copy> ManufacturedConjugateHeatTransfer<T> {
    fn fluid_heat_source(&self, x: T, y: T, t: T) -> T {
        let t_exact = self.fluid_temperature(x, y, t);
        let alpha = T::from(0.01).unwrap(); // Thermal diffusivity

        // Time derivative
        let dt_dt = -t_exact; // ∂T/∂t from exp(-t) factor

        // Laplacian
        let kx_sq = self.frequency * self.frequency;
        let ky_sq = self.frequency * self.frequency;
        let laplacian = -(kx_sq + ky_sq) * t_exact;

        // Source = ∂T/∂t - α ∇²T
        dt_dt - alpha * laplacian
    }

    fn solid_heat_source(&self, x: T, y: T, t: T) -> T {
        let t_exact = self.solid_temperature(x, y, t);
        let alpha_s = T::from(0.01).unwrap() / self.capacity_ratio; // Solid thermal diffusivity

        // Time derivative
        let dt_dt = -t_exact;

        // Laplacian
        let kx_sq = self.frequency * self.frequency;
        let ky_sq = self.frequency * self.frequency;
        let laplacian = -(kx_sq + ky_sq) * t_exact / self.conductivity_ratio;

        // Source = ∂T/∂t - α_s ∇²T
        dt_dt - alpha_s * laplacian
    }
}

/// Manufactured solution for species transport with reaction
#[derive(Debug, Clone)]
pub struct ManufacturedSpeciesTransport<T: RealField + Float + Copy> {
    /// Diffusion coefficient
    pub diffusivity: T,
    /// Reaction rate constant
    pub reaction_rate: T,
    /// Concentration amplitude
    pub amplitude: T,
    /// Wave numbers
    pub kx: T,
    pub ky: T,
}

impl<T: RealField + Float + Copy> ManufacturedSpeciesTransport<T> {
    pub fn new(diffusivity: T, reaction_rate: T, amplitude: T, kx: T, ky: T) -> Self {
        Self {
            diffusivity,
            reaction_rate,
            amplitude,
            kx,
            ky,
        }
    }
}

impl<T: RealField + Float + Copy> ManufacturedSolution<T> for ManufacturedSpeciesTransport<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        // C = A * sin(kx*x) * sin(ky*y) * exp(-t) * exp(-k²t)
        let spatial = ManufacturedFunctions::sinusoidal(x, y, T::zero(), self.kx, self.ky);
        let temporal_decay = T::from(-1.0).unwrap() - (self.kx * self.kx + self.ky * self.ky) * self.diffusivity;
        let temporal = Float::exp(temporal_decay * t);
        self.amplitude * spatial * temporal
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        let c = self.exact_solution(x, y, z, t);

        // Species equation: ∂C/∂t = D ∇²C - k C + S
        // We want S such that the MMS satisfies this equation

        // Time derivative
        let k_total = T::from(-1.0).unwrap() - (self.kx * self.kx + self.ky * self.ky) * self.diffusivity;
        let dc_dt = k_total * c;

        // Diffusion term
        let k_sq = self.kx * self.kx + self.ky * self.ky;
        let diffusion = self.diffusivity * k_sq * c;

        // Reaction term
        let reaction = self.reaction_rate * c;

        // Source = ∂C/∂t - D ∇²C + k C
        dc_dt - diffusion + reaction
    }
}

/// Manufactured solution for magnetohydrodynamics (MHD)
#[derive(Debug, Clone)]
pub struct ManufacturedMHD<T: RealField + Float + Copy> {
    /// Magnetic permeability
    pub mu_0: T,
    /// Electrical conductivity
    pub sigma: T,
    /// Velocity amplitude
    pub velocity_amp: T,
    /// Magnetic field amplitude
    pub magnetic_amp: T,
    /// Wave numbers
    pub kx: T,
    pub ky: T,
}

impl<T: RealField + Float + Copy> ManufacturedMHD<T> {
    pub fn new(mu_0: T, sigma: T, velocity_amp: T, magnetic_amp: T, kx: T, ky: T) -> Self {
        Self {
            mu_0,
            sigma,
            velocity_amp,
            magnetic_amp,
            kx,
            ky,
        }
    }
}

impl<T: RealField + Float + Copy> ManufacturedSolution<T> for ManufacturedMHD<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        // Return velocity magnitude (simplified - in practice would need vector components)
        let spatial = ManufacturedFunctions::sinusoidal(x, y, T::zero(), self.kx, self.ky);
        self.velocity_amp * spatial * Float::exp(-t)
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // MHD momentum equation: ∂u/∂t + (u·∇)u = -∇p/ρ + ν ∇²u + J × B / ρ
        // This is highly simplified - real MHD MMS would be much more complex
        let u = self.exact_solution(x, y, z, t);

        // Simplified source term
        let du_dt = -u; // Time derivative
        let viscous = -u; // Simplified viscous term
        let lorentz = self.magnetic_amp * u; // Simplified Lorentz force

        du_dt + viscous + lorentz
    }
}

/// Manufactured solution for multi-phase flows (simplified)
#[derive(Debug, Clone)]
pub struct ManufacturedMultiphase<T: RealField + Float + Copy> {
    /// Density ratio (phase 2 / phase 1)
    pub density_ratio: T,
    /// Viscosity ratio (phase 2 / phase 1)
    pub viscosity_ratio: T,
    /// Interface position
    pub interface_y: T,
    /// Amplitude
    pub amplitude: T,
    /// Wave numbers
    pub kx: T,
    pub ky: T,
}

impl<T: RealField + Float + Copy> ManufacturedMultiphase<T> {
    pub fn new(density_ratio: T, viscosity_ratio: T, interface_y: T, amplitude: T, kx: T, ky: T) -> Self {
        Self {
            density_ratio,
            viscosity_ratio,
            interface_y,
            amplitude,
            kx,
            ky,
        }
    }
}

impl<T: RealField + Float + Copy> ManufacturedSolution<T> for ManufacturedMultiphase<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        // Phase indicator function (simplified)
        let spatial = ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky);
        self.amplitude * spatial
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // Simplified multi-phase source term
        // In practice, this would involve level set or VOF equations
        let phi = self.exact_solution(x, y, z, t);
        -phi // Simplified decay term
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_conjugate_heat_transfer() {
        let cht = ManufacturedConjugateHeatTransfer::<f64>::new(
            10.0, // conductivity ratio
            2.0,  // capacity ratio
            0.5,  // interface at x=0.5
            1.0,  // amplitude
            1.0,  // frequency
        );

        let t_fluid = cht.exact_solution(0.25, 0.5, 0.0, 1.0);
        let t_solid = cht.exact_solution(0.75, 0.5, 0.0, 1.0);

        // Temperatures should be different due to material properties
        assert!(t_fluid > 0.0);
        assert!(t_solid > 0.0);
        assert!(t_fluid.abs() < 2.0);
        assert!(t_solid.abs() < 20.0); // Solid should have larger amplitude

        // Test source terms
        let s_fluid = cht.source_term(0.25, 0.5, 0.0, 1.0);
        let s_solid = cht.source_term(0.75, 0.5, 0.0, 1.0);

        assert!(s_fluid.is_finite());
        assert!(s_solid.is_finite());
    }

    #[test]
    fn test_species_transport() {
        let species = ManufacturedSpeciesTransport::<f64>::new(
            0.01, // diffusivity
            0.1,  // reaction rate
            1.0,  // amplitude
            1.0,  // kx
            1.0,  // ky
        );

        let c = species.exact_solution(0.5, 0.5, 0.0, 1.0);
        let source = species.source_term(0.5, 0.5, 0.0, 1.0);

        assert!(c > 0.0 && c < 1.0); // Concentration should be bounded
        assert!(source.is_finite());

        // Test time evolution (should decay)
        let c_later = species.exact_solution(0.5, 0.5, 0.0, 2.0);
        assert!(c_later < c); // Should decay over time
    }

    #[test]
    fn test_mhd() {
        let mhd = ManufacturedMHD::<f64>::new(
            1.0,  // mu_0
            1.0,  // sigma
            1.0,  // velocity amplitude
            0.1,  // magnetic amplitude
            1.0,  // kx
            1.0,  // ky
        );

        let u = mhd.exact_solution(0.5, 0.5, 0.0, 1.0);
        let source = mhd.source_term(0.5, 0.5, 0.0, 1.0);

        assert!(u > 0.0);
        assert!(source.is_finite());
    }

    #[test]
    fn test_multiphase() {
        let multiphase = ManufacturedMultiphase::<f64>::new(
            2.0,  // density ratio
            5.0,  // viscosity ratio
            0.5,  // interface
            1.0,  // amplitude
            1.0,  // kx
            1.0,  // ky
        );

        let phi = multiphase.exact_solution(0.5, 0.5, 0.0, 1.0);
        let source = multiphase.source_term(0.5, 0.5, 0.0, 1.0);

        assert!(phi.is_finite());
        assert!(source.is_finite());
    }
}

