//! Manufactured solutions for turbulent flow validation
//!
//! Implements manufactured solutions for turbulent Navier-Stokes equations
//! including k-ε, k-ω, and Spalart-Allmaras turbulence models.

use super::{ManufacturedSolution, ManufacturedFunctions};
use nalgebra::RealField;
use num_traits::Float;

/// Manufactured solution for k-ε turbulence model
///
/// Provides analytical solutions for turbulent kinetic energy (k)
/// and dissipation rate (ε) that satisfy the k-ε transport equations.
#[derive(Debug, Clone)]
pub struct ManufacturedKEpsilon<T: RealField + Float + Copy> {
    /// Wave number in x-direction
    pub kx: T,
    /// Wave number in y-direction
    pub ky: T,
    /// Amplitude of the solution
    pub amplitude: T,
    /// Turbulent viscosity
    pub nu_t: T,
    /// Production scaling
    pub production_scale: T,
}

impl<T: RealField + Float + Copy> ManufacturedKEpsilon<T> {
    /// Create a new manufactured k-ε solution
    pub fn new(kx: T, ky: T, amplitude: T, nu_t: T) -> Self {
        Self {
            kx,
            ky,
            amplitude,
            nu_t,
            production_scale: T::from(0.1).unwrap(), // Reasonable production scaling
        }
    }
}

impl<T: RealField + Float + Copy> ManufacturedSolution<T> for ManufacturedKEpsilon<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        // Turbulent kinetic energy: k = A * sin(kx*x) * sin(ky*y) * exp(-t)
        ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky) * self.amplitude
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // Source term for k-equation: P_k - ε + ∇·(ν_t ∇k) - ∂k/∂t
        let k = self.exact_solution(x, y, z, t);

        // Time derivative: ∂k/∂t (exact from manufactured solution)
        let dk_dt = -k; // Since k ∝ exp(-t)

        // Diffusion term: ∇·(ν_t ∇k) (exact Laplacian)
        let kx_sq = self.kx * self.kx;
        let ky_sq = self.ky * self.ky;
        let laplacian_coeff = self.nu_t * (kx_sq + ky_sq);
        let diffusion = laplacian_coeff * k;

        // Exact production term: P_k = -⟨u_i'u_j'⟩ ∂U_i/∂x_j
        // For MMS, we need the exact strain rate tensor from the manufactured velocity field
        let arg = self.kx * x + self.ky * y + self.omega * t;

        // Velocity gradients (exact derivatives of manufactured solution)
        let du_dx = self.kx * T::cos(arg);
        let du_dy = -self.kx * T::sin(arg);
        let dv_dx = self.ky * T::cos(arg);
        let dv_dy = -self.ky * T::sin(arg);

        // Strain rate tensor components: S_ij = (1/2)(∂U_i/∂x_j + ∂U_j/∂x_i)
        let s_xx = du_dx;
        let s_xy = T::from_f64(0.5).unwrap() * (du_dy + dv_dx);
        let s_yy = dv_dy;

        // Production term: P_k = 2ν_t * S_ij * S_ij (trace of S·S)
        let strain_rate_magnitude_sq = s_xx * s_xx + T::from_f64(2.0).unwrap() * s_xy * s_xy + s_yy * s_yy;
        let production = T::from_f64(2.0).unwrap() * self.nu_t * strain_rate_magnitude_sq;

        // Exact dissipation rate from k-ε model equilibrium
        // In equilibrium: ε = P_k * C_ε2 / (C_ε1 - 1 + C_ε2)
        // But for MMS validation, we need to compute exact ε
        let c_mu = T::from_f64(0.09).unwrap();
        let c_eps1 = T::from_f64(1.44).unwrap();
        let c_eps2 = T::from_f64(1.92).unwrap();

        // For manufactured solution, ε should satisfy the exact transport equation
        // ε = P_k - ∇·(ν_t ∇k) + ∂k/∂t (from rearranging k-equation)
        let epsilon = production - diffusion + dk_dt;

        // Ensure ε is positive and physically reasonable
        let epsilon = epsilon.max(T::from_f64(1e-10).unwrap());

        // Source = P_k - ε + ∇·(ν_t ∇k) - ∂k/∂t
        production - epsilon + diffusion - dk_dt
    }
}

/// Manufactured solution for k-ω turbulence model
#[derive(Debug, Clone)]
pub struct ManufacturedKOmega<T: RealField + Float + Copy> {
    /// Wave number in x-direction
    pub kx: T,
    /// Wave number in y-direction
    pub ky: T,
    /// Amplitude for k
    pub k_amplitude: T,
    /// Amplitude for ω
    pub omega_amplitude: T,
    /// Turbulent viscosity
    pub nu_t: T,
}

impl<T: RealField + Float + Copy> ManufacturedKOmega<T> {
    pub fn new(kx: T, ky: T, k_amplitude: T, omega_amplitude: T, nu_t: T) -> Self {
        Self {
            kx,
            ky,
            k_amplitude,
            omega_amplitude,
            nu_t,
        }
    }
}

impl<T: RealField + Float + Copy> ManufacturedSolution<T> for ManufacturedKOmega<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        // Return k for now - in practice would need to specify which quantity
        ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky) * self.k_amplitude
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // Source term for k-equation in k-ω model
        let k = ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky) * self.k_amplitude;
        let omega = ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky) * self.omega_amplitude;

        // Time derivative
        let dk_dt = -k;

        // Production: P_k = ν_t |∇U|²
        let strain_rate_sq = self.kx * self.kx + self.ky * self.ky;
        let production = self.nu_t * strain_rate_sq * T::from(0.1).unwrap();

        // Dissipation: β* k ω
        let beta = T::from(0.09).unwrap(); // Standard k-ω constant
        let dissipation = beta * k * omega;

        // Diffusion: ∇·(ν_t ∇k)
        let kx_sq = self.kx * self.kx;
        let ky_sq = self.ky * self.ky;
        let diffusion = self.nu_t * (kx_sq + ky_sq) * k;

        production - dissipation + diffusion - dk_dt
    }
}

/// Manufactured solution for Spalart-Allmaras turbulence model
#[derive(Debug, Clone)]
pub struct ManufacturedSpalartAllmaras<T: RealField + Float + Copy> {
    /// Wave number in x-direction
    pub kx: T,
    /// Wave number in y-direction
    pub ky: T,
    /// Amplitude for ν̃
    pub amplitude: T,
    /// Wall distance scaling
    pub wall_distance: T,
}

impl<T: RealField + Float + Copy> ManufacturedSpalartAllmaras<T> {
    pub fn new(kx: T, ky: T, amplitude: T, wall_distance: T) -> Self {
        Self {
            kx,
            ky,
            amplitude,
            wall_distance,
        }
    }
}

impl<T: RealField + Float + Copy> ManufacturedSolution<T> for ManufacturedSpalartAllmaras<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        // Modified vorticity ν̃ = ν + ν_t
        let base = ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky);
        let wall_factor = if self.wall_distance < T::from(10.0).unwrap() { self.wall_distance } else { T::from(10.0).unwrap() }; // Damping near wall
        base * self.amplitude * wall_factor
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // Source term for Spalart-Allmaras equation
        let nu_tilde = self.exact_solution(x, y, z, t);

        // Time derivative
        let dnu_dt = -nu_tilde;

        // Production term: C_b1 (1 - f_t2) ν̃ |Ω̃|
        let cb1 = T::from(0.1355).unwrap();
        let vorticity = self.kx * self.kx + self.ky * self.ky; // |Ω̃| approximation
        let production = cb1 * nu_tilde * num_traits::Float::sqrt(vorticity);

        // Destruction term: C_w1 f_w (ν̃/d)²
        let cw1 = T::from(3.239067816775729).unwrap();
        let kappa = T::from(0.41).unwrap();
        let d = if self.wall_distance > T::from(1e-10).unwrap() { self.wall_distance } else { T::from(1e-10).unwrap() };
        let chi = nu_tilde / d;
        let chi3 = chi * chi * chi;
        let denom = chi3 + T::from(7.1).unwrap() * T::from(7.1).unwrap() * T::from(7.1).unwrap();
        let fv1 = num_traits::Float::cbrt(chi3 / denom);
        let fv2 = T::one() - chi / (T::one() + chi * fv1);
        let s = vorticity + nu_tilde * fv2 / (kappa * kappa * d * d);
        let r = if nu_tilde / (s * kappa * kappa * d * d) < T::from(10.0).unwrap() {
            nu_tilde / (s * kappa * kappa * d * d)
        } else {
            T::from(10.0).unwrap()
        };
        let g = r + T::from(1.0).unwrap() / (T::from(9.0).unwrap() + r * r);
        let fw_inner = T::one() + T::from(6.0).unwrap() * T::from(6.0).unwrap();
        let fw_denom = fw_inner + g * g;
        let fw = g * num_traits::Float::sqrt(fw_inner / fw_denom);

        let destruction = cw1 * fw * nu_tilde * nu_tilde / (d * d);

        // Diffusion term: ∇·((ν + ν̃) ∇ν̃)
        let kx_sq = self.kx * self.kx;
        let ky_sq = self.ky * self.ky;
        let nu = T::from(1e-5).unwrap(); // Molecular viscosity
        let diffusion_coeff = nu + nu_tilde;
        let diffusion = diffusion_coeff * (kx_sq + ky_sq) * nu_tilde;

        production - destruction + diffusion - dnu_dt
    }
}

/// Manufactured solution for Reynolds stress components
#[derive(Debug, Clone)]
pub struct ManufacturedReynoldsStress<T: RealField + Float + Copy> {
    /// Wave numbers
    pub kx: T,
    pub ky: T,
    /// Amplitudes for different stress components
    pub uu_amp: T,
    pub vv_amp: T,
    pub uv_amp: T,
    /// Mean strain rate
    pub strain_rate: T,
}

impl<T: RealField + Float + Copy> ManufacturedReynoldsStress<T> {
    pub fn new(kx: T, ky: T, uu_amp: T, vv_amp: T, uv_amp: T, strain_rate: T) -> Self {
        Self {
            kx,
            ky,
            uu_amp,
            vv_amp,
            uv_amp,
            strain_rate,
        }
    }
}

impl<T: RealField + Float + Copy> ManufacturedSolution<T> for ManufacturedReynoldsStress<T> {
    fn exact_solution(&self, x: T, y: T, z: T, t: T) -> T {
        // Return -uv (Reynolds shear stress) as primary quantity
        ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky) * self.uv_amp
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // Source term for Reynolds stress transport equation
        // Simplified: Production - Dissipation + Diffusion - Time derivative

        let uv = self.exact_solution(x, y, z, t);
        let duv_dt = -uv; // Time derivative

        // Production: -2ν_t S_ij (simplified for 2D)
        let production = T::from(-2.0).unwrap() * self.strain_rate * uv;

        // Dissipation: ε/3 * 2k/3 (simplified)
        let dissipation = uv * T::from(0.1).unwrap(); // ε ≈ 0.1k

        // Diffusion: ∇·(ν_t ∇(-uv))
        let kx_sq = self.kx * self.kx;
        let ky_sq = self.ky * self.ky;
        let nu_t = T::from(0.01).unwrap(); // Turbulent viscosity
        let diffusion = nu_t * (kx_sq + ky_sq) * uv;

        production - dissipation + diffusion - duv_dt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_k_epsilon_manufactured_solution() {
        let mms = ManufacturedKEpsilon::<f64>::new(1.0, 1.0, 1.0, 0.01);

        let x = 0.5;
        let y = 0.5;
        let t = 1.0;

        let k = mms.exact_solution(x, y, 0.0, t);
        let source = mms.source_term(x, y, 0.0, t);

        // Verify solution is reasonable
        assert!(k > 0.0);
        assert!(source.is_finite());

        println!("k-ε MMS: k={:.6}, source={:.6}", k, source);
    }

    #[test]
    fn test_k_omega_manufactured_solution() {
        let mms = ManufacturedKOmega::<f64>::new(1.0, 1.0, 1.0, 10.0, 0.01);

        let x = 0.5;
        let y = 0.5;
        let t = 1.0;

        let solution = mms.exact_solution(x, y, 0.0, t);
        let source = mms.source_term(x, y, 0.0, t);

        assert!(solution > 0.0);
        assert!(source.is_finite());

        println!("k-ω MMS: solution={:.6}, source={:.6}", solution, source);
    }

    #[test]
    fn test_spalart_allmaras_manufactured_solution() {
        let mms = ManufacturedSpalartAllmaras::<f64>::new(1.0, 1.0, 0.1, 0.01);

        let x = 0.5;
        let y = 0.5;
        let t = 1.0;

        let nu_tilde = mms.exact_solution(x, y, 0.0, t);
        let source = mms.source_term(x, y, 0.0, t);

        assert!(nu_tilde >= 0.0); // ν̃ should be non-negative
        assert!(source.is_finite());

        println!("SA MMS: ν̃={:.6}, source={:.6}", nu_tilde, source);
    }

    #[test]
    fn test_reynolds_stress_manufactured_solution() {
        let mms = ManufacturedReynoldsStress::<f64>::new(1.0, 1.0, 1.0, 1.0, 0.5, 1.0);

        let x = 0.5;
        let y = 0.5;
        let t = 1.0;

        let uv = mms.exact_solution(x, y, 0.0, t);
        let source = mms.source_term(x, y, 0.0, t);

        assert!(uv.is_finite());
        assert!(source.is_finite());

        println!("Reynolds stress MMS: -uv={:.6}, source={:.6}", uv, source);
    }
}
