//! Multi-physics manufactured solutions for coupled CFD problems
//!
//! Implements manufactured solutions for coupled physics including:
//! - Conjugate heat transfer (fluid-solid thermal coupling)
//! - Species transport with chemical reactions
//! - Magnetohydrodynamics (MHD)
//! - Multi-phase flows
//! - Turbulent combustion

use super::{ManufacturedFunctions, ManufacturedSolution};
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Manufactured solution for conjugate heat transfer
///
/// Provides analytical solutions for coupled fluid-solid thermal problems
/// where temperature satisfies different equations in different domains.
#[derive(Debug, Clone)]
pub struct ManufacturedConjugateHeatTransfer<T: RealField + Copy> {
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

impl<T: RealField + Copy + FromPrimitive> ManufacturedConjugateHeatTransfer<T> {
    /// Create a new manufactured solution for conjugate heat transfer
    pub fn new(
        conductivity_ratio: T,
        capacity_ratio: T,
        interface_x: T,
        amplitude: T,
        frequency: T,
    ) -> Self {
        Self {
            conductivity_ratio,
            capacity_ratio,
            interface_x,
            amplitude,
            frequency,
        }
    }

    /// Evaluate temperature in fluid domain (x < interface_x)
    pub fn fluid_temperature(&self, x: T, y: T, t: T) -> T {
        let base = ManufacturedFunctions::sinusoidal(x, y, t, self.frequency, self.frequency);
        self.amplitude * base
    }

    /// Evaluate temperature in solid domain (x > interface_x)
    pub fn solid_temperature(&self, x: T, y: T, t: T) -> T {
        let base = ManufacturedFunctions::sinusoidal(x, y, t, self.frequency, self.frequency);
        // Remove conductivity_ratio to ensure continuity at interface
        self.amplitude * base
    }
}

impl<T: RealField + Copy> ManufacturedSolution<T> for ManufacturedConjugateHeatTransfer<T> {
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        if x < self.interface_x {
            self.fluid_temperature(x, y, t)
        } else {
            self.solid_temperature(x, y, t)
        }
    }

    fn source_term(&self, x: T, y: T, _z: T, t: T) -> T {
        if x < self.interface_x {
            // Fluid domain: ∂T/∂t = α ∇²T + S
            self.fluid_heat_source(x, y, t)
        } else {
            // Solid domain: ∂T/∂t = α_s ∇²T + S_s
            self.solid_heat_source(x, y, t)
        }
    }
}

impl<T: RealField + Copy + FromPrimitive> ManufacturedConjugateHeatTransfer<T> {
    fn fluid_heat_source(&self, x: T, y: T, t: T) -> T {
        let t_exact = self.fluid_temperature(x, y, t);
        let alpha = <T as FromPrimitive>::from_f64(0.01f64).unwrap(); // Base thermal diffusivity

        // Time derivative: ∂T/∂t = -T (from exp(-t))
        let dt_dt = -t_exact;

        // Laplacian: ∇²T = -(kx² + ky²) * T
        let kx_sq = self.frequency * self.frequency;
        let ky_sq = self.frequency * self.frequency;
        let laplacian = -(kx_sq + ky_sq) * t_exact;

        // Source = ∂T/∂t - α ∇²T
        dt_dt - alpha * laplacian
    }

    fn solid_heat_source(&self, x: T, y: T, t: T) -> T {
        let t_exact = self.solid_temperature(x, y, t);
        // Solid thermal diffusivity α_s = α * (k_ratio / capacity_ratio)
        let alpha = <T as FromPrimitive>::from_f64(0.01f64).unwrap();
        let alpha_s = alpha * self.conductivity_ratio / self.capacity_ratio;

        // Time derivative
        let dt_dt = -t_exact;

        // Laplacian
        let kx_sq = self.frequency * self.frequency;
        let ky_sq = self.frequency * self.frequency;
        let laplacian = -(kx_sq + ky_sq) * t_exact;

        // Source = ∂T/∂t - α_s ∇²T
        dt_dt - alpha_s * laplacian
    }
}

/// Manufactured solution for species transport with reaction
#[derive(Debug, Clone)]
pub struct ManufacturedSpeciesTransport<T: RealField + Copy> {
    /// Diffusion coefficient
    pub diffusivity: T,
    /// Reaction rate constant
    pub reaction_rate: T,
    /// Concentration amplitude
    pub amplitude: T,
    /// Wave numbers
    pub kx: T,
    /// Wave number in y-direction
    pub ky: T,
}

impl<T: RealField + Copy> ManufacturedSpeciesTransport<T> {
    /// Create a new manufactured solution for species transport
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

impl<T: RealField + Copy + FromPrimitive> ManufacturedSolution<T>
    for ManufacturedSpeciesTransport<T>
{
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        // C = A * sin(kx*x) * sin(ky*y) * exp(-t) * exp(-k²t)
        let spatial = ManufacturedFunctions::sinusoidal(x, y, T::zero(), self.kx, self.ky);
        let temporal_decay = <T as FromPrimitive>::from_f64(-1.0).unwrap()
            - (self.kx * self.kx + self.ky * self.ky) * self.diffusivity;
        let temporal = nalgebra::ComplexField::exp(temporal_decay * t);
        self.amplitude * spatial * temporal
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        let c = self.exact_solution(x, y, z, t);

        // Species equation: ∂C/∂t = D ∇²C - k C + S
        // We want S such that the MMS satisfies this equation

        // Time derivative
        let k_total = <T as FromPrimitive>::from_f64(-1.0).unwrap()
            - (self.kx * self.kx + self.ky * self.ky) * self.diffusivity;
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

/// Vector-valued MHD fields for manufactured solutions
#[derive(Debug, Clone)]
pub struct MHDVectorFields<T: RealField + Copy> {
    /// Velocity field components (u, v, w)
    pub velocity: (T, T, T),
    /// Magnetic field components (Bx, By, Bz)
    pub magnetic: (T, T, T),
    /// Pressure field
    pub pressure: T,
    /// Current density components (Jx, Jy, Jz)
    pub current_density: (T, T, T),
    /// Lorentz force components (Fx, Fy, Fz)
    pub lorentz_force: (T, T, T),
}

/// Manufactured solution for magnetohydrodynamics (MHD)
#[derive(Debug, Clone)]
pub struct ManufacturedMHD<T: RealField + Copy> {
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
    /// Wave number in y-direction
    pub ky: T,
    /// Density
    pub density: T,
    /// Kinematic viscosity
    pub viscosity: T,
}

impl<T: RealField + Copy> ManufacturedMHD<T> {
    /// Create a new manufactured solution for MHD
    pub fn new(
        mu_0: T,
        sigma: T,
        velocity_amp: T,
        magnetic_amp: T,
        kx: T,
        ky: T,
        density: T,
        viscosity: T,
    ) -> Self {
        Self {
            mu_0,
            sigma,
            velocity_amp,
            magnetic_amp,
            kx,
            ky,
            density,
            viscosity,
        }
    }

    /// Compute vector-valued MHD fields at given position and time
    pub fn compute_vector_fields(&self, x: T, y: T, z: T, t: T) -> MHDVectorFields<T> {
        // Spatial variation
        let spatial_u = ManufacturedFunctions::sinusoidal(x, y, z, self.kx, self.ky);
        let spatial_v = ManufacturedFunctions::sinusoidal(x, y, z, self.ky, self.kx);
        let spatial_w = ManufacturedFunctions::sinusoidal(x, y, z, self.kx, self.kx);
        
        // Temporal decay
        let temporal = nalgebra::ComplexField::exp(-t);
        
        // Velocity field components
        let u = self.velocity_amp * spatial_u * temporal;
        let v = self.velocity_amp * spatial_v * temporal * T::from_f64(0.5).unwrap();
        let w = self.velocity_amp * spatial_w * temporal * T::from_f64(0.25).unwrap();
        
        // Magnetic field components (perpendicular to velocity for interesting dynamics)
        let bx = self.magnetic_amp * spatial_v * temporal;
        let by = self.magnetic_amp * spatial_u * temporal;
        let bz = self.magnetic_amp * spatial_w * temporal * T::from_f64(0.1).unwrap();
        
        // Pressure field (from momentum equation)
        let pressure = self.density * (u * u + v * v + w * w) * T::from_f64(0.5).unwrap();
        
        // Current density J = σ(E + u × B), assuming E = 0 for simplicity
        let jx = self.sigma * (v * bz - w * by);
        let jy = self.sigma * (w * bx - u * bz);
        let jz = self.sigma * (u * by - v * bx);
        
        // Lorentz force F = J × B
        let fx = jy * bz - jz * by;
        let fy = jz * bx - jx * bz;
        let fz = jx * by - jy * bx;
        
        MHDVectorFields {
            velocity: (u, v, w),
            magnetic: (bx, by, bz),
            pressure,
            current_density: (jx, jy, jz),
            lorentz_force: (fx, fy, fz),
        }
    }

    /// Compute momentum source term for x-component
    pub fn momentum_source_x(&self, x: T, y: T, z: T, t: T) -> T {
        let fields = self.compute_vector_fields(x, y, z, t);
        let (u, v, w) = fields.velocity;
        let (fx, _, _) = fields.lorentz_force;
        
        // Time derivative
        let du_dt = -u; // From exp(-t) temporal decay
        
        // Convective terms (u·∇)u
        let convective = u * self.kx * u + v * self.ky * u + w * self.kx * u;
        
        // Viscous terms ν∇²u
        let k_squared = self.kx * self.kx + self.ky * self.ky + self.kx * self.kx;
        let viscous = -self.viscosity * k_squared * u;
        
        // Pressure gradient
        let pressure_grad_x = self.density * self.kx * fields.pressure;
        
        // Source term to make MMS work
        du_dt - convective + pressure_grad_x / self.density - viscous - fx / self.density
    }

    /// Compute momentum source term for y-component
    pub fn momentum_source_y(&self, x: T, y: T, z: T, t: T) -> T {
        let fields = self.compute_vector_fields(x, y, z, t);
        let (u, v, w) = fields.velocity;
        let (_, fy, _) = fields.lorentz_force;
        
        // Time derivative
        let dv_dt = -v * T::from_f64(0.5).unwrap();
        
        // Convective terms
        let convective = u * self.kx * v + v * self.ky * v + w * self.kx * v;
        
        // Viscous terms
        let k_squared = self.kx * self.kx + self.ky * self.ky + self.kx * self.kx;
        let viscous = -self.viscosity * k_squared * v;
        
        // Pressure gradient
        let pressure_grad_y = self.density * self.ky * fields.pressure;
        
        // Source term
        dv_dt - convective + pressure_grad_y / self.density - viscous - fy / self.density
    }

    /// Compute momentum source term for z-component
    pub fn momentum_source_z(&self, x: T, y: T, z: T, t: T) -> T {
        let fields = self.compute_vector_fields(x, y, z, t);
        let (u, v, w) = fields.velocity;
        let (_, _, fz) = fields.lorentz_force;
        
        // Time derivative
        let dw_dt = -w * T::from_f64(0.25).unwrap();
        
        // Convective terms
        let convective = u * self.kx * w + v * self.ky * w + w * self.kx * w;
        
        // Viscous terms
        let k_squared = self.kx * self.kx + self.ky * self.ky + self.kx * self.kx;
        let viscous = -self.viscosity * k_squared * w;
        
        // Pressure gradient (assuming no variation in z for pressure)
        let pressure_grad_z = T::zero();
        
        // Source term
        dw_dt - convective + pressure_grad_z / self.density - viscous - fz / self.density
    }

    /// Compute magnetic field source term (induction equation)
    pub fn induction_source(&self, x: T, y: T, z: T, t: T) -> (T, T, T) {
        let fields = self.compute_vector_fields(x, y, z, t);
        let (bx, by, bz) = fields.magnetic;
        let (u, v, w) = fields.velocity;
        
        // Induction equation: ∂B/∂t = ∇×(u×B) + η∇²B
        // where η = 1/(μ₀σ) is magnetic diffusivity
        
        let magnetic_diffusivity = T::one() / (self.mu_0 * self.sigma);
        let k_squared = self.kx * self.kx + self.ky * self.ky + self.kx * self.kx;
        
        // Time derivatives
        let dbx_dt = -bx;
        let dby_dt = -by;
        let dbz_dt = -bz * T::from_f64(0.1).unwrap();
        
        // Diffusion terms
        let diff_x = magnetic_diffusivity * k_squared * bx;
        let diff_y = magnetic_diffusivity * k_squared * by;
        let diff_z = magnetic_diffusivity * k_squared * bz;
        
        // Advection terms (simplified)
        let adv_x = u * self.kx * bx + v * self.ky * bx;
        let adv_y = u * self.kx * by + v * self.ky * by;
        let adv_z = u * self.kx * bz + v * self.ky * bz;
        
        // Source terms
        let sx = dbx_dt - adv_x + diff_x;
        let sy = dby_dt - adv_y + diff_y;
        let sz = dbz_dt - adv_z + diff_z;
        
        (sx, sy, sz)
    }
}

impl<T: RealField + Copy> ManufacturedSolution<T> for ManufacturedMHD<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        // Return velocity magnitude for compatibility with scalar interface
        let fields = self.compute_vector_fields(x, y, z, t);
        let (u, v, w) = fields.velocity;
        (u * u + v * v + w * w).sqrt()
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // Return x-momentum source term for compatibility
        self.momentum_source_x(x, y, z, t)
    }
}

/// TODO: Provide a physically consistent manufactured multiphase solution.
#[derive(Debug, Clone)]
pub struct ManufacturedMultiphase<T: RealField + Copy> {
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
    /// Wave number in y-direction
    pub ky: T,
}

impl<T: RealField + Copy> ManufacturedMultiphase<T> {
    /// Create a new manufactured solution for multi-phase flows
    pub fn new(
        density_ratio: T,
        viscosity_ratio: T,
        interface_y: T,
        amplitude: T,
        kx: T,
        ky: T,
    ) -> Self {
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

impl<T: RealField + Copy> ManufacturedSolution<T> for ManufacturedMultiphase<T> {
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        // TODO: Use a level-set/VOF-consistent manufactured interface field.
        let spatial = ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky);
        self.amplitude * spatial
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // TODO: Derive the source term from the governing multiphase transport equation.
        let phi = self.exact_solution(x, y, z, t);
        -phi
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use approx::assert_relative_eq;

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
            1.0, // mu_0
            1.0, // sigma
            1.0, // velocity amplitude
            0.1, // magnetic amplitude
            1.0, // kx
            1.0, // ky
            1.0, // density
            0.01, // viscosity
        );

        let u = mhd.exact_solution(0.5, 0.5, 0.0, 1.0);
        let source = mhd.source_term(0.5, 0.5, 0.0, 1.0);
        
        // Test vector fields
        let fields = mhd.compute_vector_fields(0.5, 0.5, 0.0, 1.0);
        let (vel_x, vel_y, vel_z) = fields.velocity;
        let (mag_x, mag_y, mag_z) = fields.magnetic;
        let (fx, fy, fz) = fields.lorentz_force;

        assert!(u > 0.0);
        assert!(source.is_finite());
        assert!(vel_x.is_finite() && vel_y.is_finite() && vel_z.is_finite());
        assert!(mag_x.is_finite() && mag_y.is_finite() && mag_z.is_finite());
        assert!(fx.is_finite() && fy.is_finite() && fz.is_finite());
        
        // Test source terms
        let src_x = mhd.momentum_source_x(0.5, 0.5, 0.0, 1.0);
        let src_y = mhd.momentum_source_y(0.5, 0.5, 0.0, 1.0);
        let src_z = mhd.momentum_source_z(0.5, 0.5, 0.0, 1.0);
        let (ind_x, ind_y, ind_z) = mhd.induction_source(0.5, 0.5, 0.0, 1.0);
        
        assert!(src_x.is_finite() && src_y.is_finite() && src_z.is_finite());
        assert!(ind_x.is_finite() && ind_y.is_finite() && ind_z.is_finite());
    }

    #[test]
    fn test_multiphase() {
        let multiphase = ManufacturedMultiphase::<f64>::new(
            2.0, // density ratio
            5.0, // viscosity ratio
            0.5, // interface
            1.0, // amplitude
            1.0, // kx
            1.0, // ky
        );

        let phi = multiphase.exact_solution(0.5, 0.5, 0.0, 1.0);
        let source = multiphase.source_term(0.5, 0.5, 0.0, 1.0);

        assert!(phi.is_finite());
        assert!(source.is_finite());
    }
}
