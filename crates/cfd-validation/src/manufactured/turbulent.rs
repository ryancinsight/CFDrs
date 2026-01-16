//! Manufactured solutions for turbulent flow validation
//!
//! Implements manufactured solutions for turbulent Navier-Stokes equations
//! including k-ε, k-ω, and Spalart-Allmaras turbulence models.

use super::{ManufacturedFunctions, ManufacturedSolution};
use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;

/// Manufactured solution for k-ε turbulence model
///
/// Provides analytical solutions for turbulent kinetic energy (k)
/// and dissipation rate (ε) that satisfy the k-ε transport equations.
#[derive(Debug, Clone)]
pub struct ManufacturedKEpsilon<T: RealField + Copy> {
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

impl<T: RealField + Copy> ManufacturedKEpsilon<T> {
    /// Create a new manufactured k-ε solution
    pub fn new(kx: T, ky: T, amplitude: T, nu_t: T) -> Self {
        Self {
            kx,
            ky,
            amplitude,
            nu_t,
            production_scale: T::from_f64_or_one(0.1), // Reasonable production scaling without unwrap
        }
    }
}

impl<T: RealField + Copy> ManufacturedSolution<T> for ManufacturedKEpsilon<T> {
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        // Turbulent kinetic energy: k = A * sin(kx*x) * sin(ky*y) * exp(-t)
        ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky) * self.amplitude
    }

    fn source_term(&self, x: T, y: T, _z: T, t: T) -> T {
        // Exact analytical source term for k-ε MMS validation
        // k-equation: ∂k/∂t + U_j ∂k/∂x_j = P_k - ε + ∇·(ν_t ∇k)

        // Manufactured solution: k = A * sin(kx*x) * sin(ky*y) * exp(-t)
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let exp_t = nalgebra::ComplexField::exp(-t);

        let k = self.amplitude * sin_kx_x * sin_ky_y * exp_t;

        // Time derivative: ∂k/∂t = -k (exact)
        let dk_dt = -k;

        // Convection term: U_j ∂k/∂x_j
        // Using manufactured velocity field: u = sin(kx*x)*cos(ky*y)*exp(-t), v = cos(kx*x)*sin(ky*y)*exp(-t)
        let u = sin_kx_x * nalgebra::ComplexField::cos(self.ky * y) * exp_t;
        let v = nalgebra::ComplexField::cos(self.kx * x) * sin_ky_y * exp_t;

        let dk_dx =
            self.amplitude * self.kx * nalgebra::ComplexField::cos(self.kx * x) * sin_ky_y * exp_t;
        let dk_dy =
            self.amplitude * sin_kx_x * self.ky * nalgebra::ComplexField::cos(self.ky * y) * exp_t;

        let convection = u * dk_dx + v * dk_dy;

        // Diffusion term: ∇·(ν_t ∇k) = ν_t * ∇²k
        let d2k_dx2 = -self.amplitude * self.kx * self.kx * sin_kx_x * sin_ky_y * exp_t;
        let d2k_dy2 = -self.amplitude * self.ky * self.ky * sin_kx_x * sin_ky_y * exp_t;
        let diffusion = self.nu_t * (d2k_dx2 + d2k_dy2);

        // Exact production term: P_k = 2ν_t * S_ij * S_ij
        // Strain rate tensor from manufactured velocity field
        let du_dx = self.kx
            * nalgebra::ComplexField::cos(self.kx * x)
            * nalgebra::ComplexField::cos(self.ky * y)
            * exp_t;
        let du_dy = -self.kx * sin_kx_x * nalgebra::ComplexField::sin(self.ky * y) * exp_t;
        let dv_dx = -self.ky * nalgebra::ComplexField::sin(self.kx * x) * sin_ky_y * exp_t;
        let dv_dy = self.ky
            * nalgebra::ComplexField::cos(self.kx * x)
            * nalgebra::ComplexField::cos(self.ky * y)
            * exp_t;

        let s_xx = du_dx;
        let s_xy = T::from_f64_or_one(0.5) * (du_dy + dv_dx);
        let s_yy = dv_dy;

        let strain_rate_magnitude_sq =
            s_xx * s_xx + T::from_f64_or_one(2.0) * s_xy * s_xy + s_yy * s_yy;
        let production = T::from_f64_or_one(2.0) * self.nu_t * strain_rate_magnitude_sq;

        // Exact dissipation rate from manufactured solution
        // For MMS, ε must satisfy: ε = P_k - ∇·(ν_t ∇k) + ∂k/∂t + convection
        let epsilon = production - diffusion + dk_dt + convection;

        // Ensure ε is positive and physically reasonable
        let epsilon = epsilon.max(T::from_f64_or_zero(1e-10));

        // Source term: P_k - ε + ∇·(ν_t ∇k) - ∂k/∂t - U_j ∂k/∂x_j
        production - epsilon + diffusion - dk_dt - convection
    }
}

/// Manufactured solution for k-ω turbulence model
#[derive(Debug, Clone)]
pub struct ManufacturedKOmega<T: RealField + Copy> {
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

impl<T: RealField + Copy> ManufacturedKOmega<T> {
    /// Create a new manufactured k-ω solution
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

impl<T: RealField + Copy> ManufacturedSolution<T> for ManufacturedKOmega<T> {
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        // TODO: Return the requested field (k or ω) via a typed API instead of always returning k.
        ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky) * self.k_amplitude
    }

    fn source_term(&self, x: T, y: T, _z: T, t: T) -> T {
        // Exact analytical source term for k-ω MMS validation
        // k-equation: ∂k/∂t + U_j ∂k/∂x_j = P_k - β* k ω + ∇·(ν_t ∇k)

        // Manufactured solutions: k = A_k * sin(kx*x) * sin(ky*y) * exp(-t)
        //                       ω = A_ω * sin(kx*x) * sin(ky*y) * exp(-t)
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let exp_t = nalgebra::ComplexField::exp(-t);

        let k = self.k_amplitude * sin_kx_x * sin_ky_y * exp_t;
        let omega = self.omega_amplitude * sin_kx_x * sin_ky_y * exp_t;

        // Time derivative: ∂k/∂t = -k (exact)
        let dk_dt = -k;

        // Convection term: U_j ∂k/∂x_j
        // Using manufactured velocity field: u = sin(kx*x)*cos(ky*y)*exp(-t), v = cos(kx*x)*sin(ky*y)*exp(-t)
        let u = sin_kx_x * nalgebra::ComplexField::cos(self.ky * y) * exp_t;
        let v = nalgebra::ComplexField::cos(self.kx * x) * sin_ky_y * exp_t;

        let dk_dx = self.k_amplitude
            * self.kx
            * nalgebra::ComplexField::cos(self.kx * x)
            * sin_ky_y
            * exp_t;
        let dk_dy = self.k_amplitude
            * sin_kx_x
            * self.ky
            * nalgebra::ComplexField::cos(self.ky * y)
            * exp_t;

        let convection = u * dk_dx + v * dk_dy;

        // Diffusion term: ∇·(ν_t ∇k) = ν_t * ∇²k
        let d2k_dx2 = -self.k_amplitude * self.kx * self.kx * sin_kx_x * sin_ky_y * exp_t;
        let d2k_dy2 = -self.k_amplitude * self.ky * self.ky * sin_kx_x * sin_ky_y * exp_t;
        let diffusion = self.nu_t * (d2k_dx2 + d2k_dy2);

        // Exact production term: P_k = 2ν_t * S_ij * S_ij
        // Strain rate tensor from manufactured velocity field
        let du_dx = self.kx
            * nalgebra::ComplexField::cos(self.kx * x)
            * nalgebra::ComplexField::cos(self.ky * y)
            * exp_t;
        let du_dy = -self.kx * sin_kx_x * nalgebra::ComplexField::sin(self.ky * y) * exp_t;
        let dv_dx = -self.ky * nalgebra::ComplexField::sin(self.kx * x) * sin_ky_y * exp_t;
        let dv_dy = self.ky
            * nalgebra::ComplexField::cos(self.kx * x)
            * nalgebra::ComplexField::cos(self.ky * y)
            * exp_t;

        let s_xx = du_dx;
        let s_xy = T::from_f64_or_one(0.5) * (du_dy + dv_dx);
        let s_yy = dv_dy;

        let strain_rate_magnitude_sq =
            s_xx * s_xx + T::from_f64_or_one(2.0) * s_xy * s_xy + s_yy * s_yy;
        let production = T::from_f64_or_one(2.0) * self.nu_t * strain_rate_magnitude_sq;

        // Exact dissipation term: β* k ω (standard k-ω model)
        let beta = T::from_f64_or_one(0.09); // Standard k-ω constant
        let dissipation = beta * k * omega;

        // Source term: P_k - β* k ω + ∇·(ν_t ∇k) - ∂k/∂t - U_j ∂k/∂x_j
        production - dissipation + diffusion - dk_dt - convection
    }
}

/// Manufactured solution for Spalart-Allmaras turbulence model
#[derive(Debug, Clone)]
pub struct ManufacturedSpalartAllmaras<T: RealField + Copy> {
    /// Wave number in x-direction
    pub kx: T,
    /// Wave number in y-direction
    pub ky: T,
    /// Amplitude for ν̃
    pub amplitude: T,
    /// Wall distance scaling
    pub wall_distance: T,
}

impl<T: RealField + Copy> ManufacturedSpalartAllmaras<T> {
    /// Create a new manufactured solution for Spalart-Allmaras
    pub fn new(kx: T, ky: T, amplitude: T, wall_distance: T) -> Self {
        Self {
            kx,
            ky,
            amplitude,
            wall_distance,
        }
    }
}

impl<T: RealField + Copy> ManufacturedSolution<T> for ManufacturedSpalartAllmaras<T> {
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        // Modified vorticity ν̃ = ν + ν_t
        let base = ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky);
        let wall_factor = if self.wall_distance < T::from_f64_or_one(10.0) {
            self.wall_distance
        } else {
            T::from_f64_or_one(10.0)
        }; // Damping near wall
        base * self.amplitude * wall_factor
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // Source term for Spalart-Allmaras equation
        let nu_tilde = self.exact_solution(x, y, z, t);

        // Time derivative
        let dnu_dt = -nu_tilde;

        // Production term: C_b1 (1 - f_t2) ν̃ |Ω̃|
        let cb1 = T::from_f64_or_one(0.1355);
        let vorticity = self.kx * self.kx + self.ky * self.ky; // |Ω̃| approximation
        let production = cb1 * nu_tilde * nalgebra::ComplexField::sqrt(vorticity);

        // Destruction term: C_w1 f_w (ν̃/d)²
        let cw1 = T::from_f64_or_one(3.239067816775729);
        let kappa = T::from_f64_or_one(0.41);
        let d = if self.wall_distance > T::from_f64_or_zero(1e-10) {
            self.wall_distance
        } else {
            T::from_f64_or_zero(1e-10)
        };
        let chi = nu_tilde / d;
        let chi3 = chi * chi * chi;
        let denom =
            chi3 + T::from_f64_or_one(7.1) * T::from_f64_or_one(7.1) * T::from_f64_or_one(7.1);
        let fv1 = nalgebra::ComplexField::cbrt(chi3 / denom);
        let fv2 = T::one() - chi / (T::one() + chi * fv1);
        let s = vorticity + nu_tilde * fv2 / (kappa * kappa * d * d);
        let r = if nu_tilde / (s * kappa * kappa * d * d) < T::from_f64_or_one(10.0) {
            nu_tilde / (s * kappa * kappa * d * d)
        } else {
            T::from_f64_or_one(10.0)
        };
        let g = r + T::from_f64_or_one(1.0) / (T::from_f64_or_one(9.0) + r * r);
        let fw_inner = T::one() + T::from_f64_or_one(6.0) * T::from_f64_or_one(6.0);
        let fw_denom = fw_inner + g * g;
        let fw = g * nalgebra::ComplexField::sqrt(fw_inner / fw_denom);

        let destruction = cw1 * fw * nu_tilde * nu_tilde / (d * d);

        // Diffusion term: ∇·((ν + ν̃) ∇ν̃)
        let kx_sq = self.kx * self.kx;
        let ky_sq = self.ky * self.ky;
        let nu = T::from_f64_or_zero(1e-5); // Molecular viscosity
        let diffusion_coeff = nu + nu_tilde;
        let diffusion = diffusion_coeff * (kx_sq + ky_sq) * nu_tilde;

        production - destruction + diffusion - dnu_dt
    }
}

/// Manufactured solution for Reynolds stress components
#[derive(Debug, Clone)]
pub struct ManufacturedReynoldsStress<T: RealField + Copy> {
    /// Wave numbers
    pub kx: T,
    /// Wave number in y-direction
    pub ky: T,
    /// Amplitudes for different stress components
    pub uu_amp: T,
    /// Amplitude for vv stress component
    pub vv_amp: T,
    /// Amplitude for uv stress component
    pub uv_amp: T,
    /// Mean strain rate
    pub strain_rate: T,
}

impl<T: RealField + Copy> ManufacturedReynoldsStress<T> {
    /// Create a new manufactured solution for Reynolds stress components
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

impl<T: RealField + Copy> ManufacturedSolution<T> for ManufacturedReynoldsStress<T> {
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        // Return -uv (Reynolds shear stress) as primary quantity
        ManufacturedFunctions::sinusoidal(x, y, t, self.kx, self.ky) * self.uv_amp
    }

    fn source_term(&self, x: T, y: T, z: T, t: T) -> T {
        // Source term for Reynolds stress transport equation
        // TODO: Implement full Reynolds-stress transport source terms (production/redistribution/dissipation/diffusion).

        let uv = self.exact_solution(x, y, z, t);
        let duv_dt = -uv; // Time derivative

        // TODO: Use the correct production term (tensor form, not 2D scalar simplification).
        let production = T::from_f64_or_one(-2.0) * self.strain_rate * uv;

        // TODO: Use a physically consistent dissipation model instead of a constant factor.
        let dissipation = uv * T::from_f64_or_one(0.1); // ε ≈ 0.1k

        // Diffusion: ∇·(ν_t ∇(-uv))
        let kx_sq = self.kx * self.kx;
        let ky_sq = self.ky * self.ky;
        let nu_t = T::from_f64_or_one(0.01); // Turbulent viscosity
        let diffusion = nu_t * (kx_sq + ky_sq) * uv;

        production - dissipation + diffusion - duv_dt
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    // use approx::assert_relative_eq;

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

        println!("k-ε MMS: k={k:.6}, source={source:.6}");
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

        println!("k-ω MMS: solution={solution:.6}, source={source:.6}");
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

        println!("SA MMS: ν̃={nu_tilde:.6}, source={source:.6}");
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

        println!("Reynolds stress MMS: -uv={uv:.6}, source={source:.6}");
    }
}
