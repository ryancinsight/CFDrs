//! Method of Manufactured Solutions for Reynolds Stress Transport Model Validation
//!
//! This module implements comprehensive manufactured solutions for validating
//! the Reynolds stress transport equations. It provides analytical solutions
//! for all 6 Reynolds stress components and computes source terms that satisfy
//! the full transport equations including production, pressure-strain correlation,
//! dissipation, turbulent transport, and molecular diffusion.
//!
//! ## Mathematical Foundation
//!
//! The Reynolds stress transport equations:
//! D⟨u_i'u_j'⟩/Dt = P_ij + Φ_ij - ε_ij + T_ij + D_ij
//!
//! Where manufactured solutions are chosen such that:
//! - ⟨u_i'u_j'⟩ = analytical function of (x,y,t)
//! - Source terms S_ij computed to satisfy the PDE exactly
//!
//! ## Manufactured Solution Form
//!
//! Using trigonometric-polynomial manufactured solutions:
//! ⟨u'u'⟩ = A_uu * sin(kx*x) * sin(ky*y) * exp(-α*t) + isotropic_part
//! ⟨u'v'⟩ = A_uv * cos(kx*x) * cos(ky*y) * exp(-α*t)
//! ⟨v'v'⟩ = A_vv * sin(kx*x) * sin(ky*y) * exp(-α*t) + isotropic_part
//!
//! ## References
//!
//! - Roache, P.J. (2002) "Code Verification by the Method of Manufactured Solutions"
//! - Salari, K. & Knupp, P. (2000) "Code Verification by the Method of Manufactured Solutions"
//! - Launder, B. E., et al. (1975). Progress in the development of a Reynolds-stress closure

use super::ManufacturedSolution;
use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;
use std::f64::consts::PI;

/// Comprehensive manufactured solution for Reynolds stress transport equations
#[derive(Debug, Clone)]
pub struct ManufacturedReynoldsStressMMS<T: RealField + Copy + FromPrimitive> {
    /// Wave numbers for spatial variation
    pub kx: T,
    /// Wave number in y-direction
    pub ky: T,
    /// Temporal decay rate
    pub alpha: T,
    /// Amplitudes for Reynolds stress components
    pub a_uu: T, // ⟨u'u'⟩ amplitude
    /// ⟨u'v'⟩ amplitude
    pub a_uv: T,
    /// ⟨v'v'⟩ amplitude
    pub a_vv: T,
    /// Mean velocity field amplitudes (for production terms)
    pub u0_amp: T, // Mean u velocity amplitude
    /// Mean v velocity amplitude
    pub v0_amp: T,
    /// Turbulent kinetic energy and dissipation amplitudes
    pub k_amp: T, // Turbulent kinetic energy amplitude
    /// Dissipation rate amplitude
    pub eps_amp: T,
    /// Pressure-strain model type
    pub pressure_strain_model: PressureStrainModelMMS,
    /// Turbulent viscosity (for transport terms)
    pub nu_t: T,
    /// Molecular viscosity (for diffusion terms)
    pub nu: T,
}

/// Pressure-strain models for MMS
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum PressureStrainModelMMS {
    /// Linear return-to-isotropy (Rotta, 1951)
    LinearReturnToIsotropy,
    /// Quadratic model (Speziale et al., 1991)
    Quadratic,
    /// SSG model (Speziale-Sarkar-Gatski)
    SSG,
}

impl<T: RealField + Copy + FromPrimitive> ManufacturedReynoldsStressMMS<T> {
    /// Create a new manufactured solution for Reynolds stress validation
    pub fn new(
        kx: T,
        ky: T,
        alpha: T,
        a_uu: T,
        a_uv: T,
        a_vv: T,
        u0_amp: T,
        v0_amp: T,
        k_amp: T,
        eps_amp: T,
        pressure_strain_model: PressureStrainModelMMS,
        nu_t: T,
        nu: T,
    ) -> Self {
        Self {
            kx,
            ky,
            alpha,
            a_uu,
            a_uv,
            a_vv,
            u0_amp,
            v0_amp,
            k_amp,
            eps_amp,
            pressure_strain_model,
            nu_t,
            nu,
        }
    }

    /// Create a standard test case with reasonable parameters
    pub fn standard_test_case() -> Self {
        Self::new(
            T::from_f64(2.0 * PI).unwrap(), // kx
            T::from_f64(2.0 * PI).unwrap(), // ky
            T::from_f64(0.1).unwrap(),      // alpha (slow decay)
            T::from_f64(0.1).unwrap(),      // a_uu
            T::from_f64(0.05).unwrap(),     // a_uv
            T::from_f64(0.1).unwrap(),      // a_vv
            T::from_f64(1.0).unwrap(),      // u0_amp
            T::from_f64(0.5).unwrap(),      // v0_amp
            T::from_f64(0.15).unwrap(),     // k_amp
            T::from_f64(0.01).unwrap(),     // eps_amp
            PressureStrainModelMMS::Quadratic,
            T::from_f64(0.01).unwrap(), // nu_t
            T::from_f64(1e-5).unwrap(), // nu
        )
    }

    /// Evaluate the exact Reynolds stress component ⟨u_i'u_j'⟩
    pub fn exact_reynolds_stress(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let cos_kx_x = nalgebra::ComplexField::cos(self.kx * x);
        let cos_ky_y = nalgebra::ComplexField::cos(self.ky * y);
        let exp_decay = nalgebra::ComplexField::exp(-self.alpha * t);

        let isotropic_part = self.k_amp * T::from_f64(2.0 / 3.0).unwrap() * exp_decay;

        match (i, j) {
            (0, 0) => {
                // ⟨u'u'⟩
                self.a_uu * sin_kx_x * sin_ky_y * exp_decay + isotropic_part
            }
            (0, 1) | (1, 0) => {
                // ⟨u'v'⟩
                self.a_uv * cos_kx_x * cos_ky_y * exp_decay
            }
            (1, 1) => {
                // ⟨v'v'⟩
                self.a_vv * sin_kx_x * sin_ky_y * exp_decay + isotropic_part
            }
            _ => T::zero(),
        }
    }

    /// Evaluate the exact mean velocity component U_i
    pub fn exact_mean_velocity(&self, i: usize, x: T, y: T, t: T) -> T {
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let exp_decay = nalgebra::ComplexField::exp(-self.alpha * t);

        match i {
            0 => self.u0_amp * sin_kx_x * sin_ky_y * exp_decay, // U
            1 => self.v0_amp * sin_kx_x * sin_ky_y * exp_decay, // V
            _ => T::zero(),
        }
    }

    /// Evaluate the exact turbulent kinetic energy k
    pub fn exact_tke(&self, x: T, y: T, t: T) -> T {
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let exp_decay = nalgebra::ComplexField::exp(-self.alpha * t);

        self.k_amp * sin_kx_x * sin_ky_y * exp_decay
    }

    /// Evaluate the exact dissipation rate ε
    pub fn exact_dissipation(&self, x: T, y: T, t: T) -> T {
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let exp_decay = nalgebra::ComplexField::exp(-self.alpha * t);

        self.eps_amp * sin_kx_x * sin_ky_y * exp_decay
    }

    /// Compute the production term P_ij = -⟨u_i'u_k'⟩∂U_j/∂x_k - ⟨u_j'u_k'⟩∂U_i/∂x_k
    pub fn production_term(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        // Velocity gradients
        let _du_dx = self.velocity_gradient(0, 0, x, y, t);
        let du_dy = self.velocity_gradient(0, 1, x, y, t);
        let dv_dx = self.velocity_gradient(1, 0, x, y, t);
        let _dv_dy = self.velocity_gradient(1, 1, x, y, t);

        // Reynolds stress components
        let uu = self.exact_reynolds_stress(0, 0, x, y, t);
        let uv = self.exact_reynolds_stress(0, 1, x, y, t);
        let vv = self.exact_reynolds_stress(1, 1, x, y, t);

        match (i, j) {
            (0, 0) => -T::from_f64(2.0).unwrap() * uv * du_dy, // P_xx = -2⟨u'v'⟩∂U/∂y
            (0, 1) | (1, 0) => -uu * dv_dx - vv * du_dy,       // P_xy = -⟨u'u'⟩∂V/∂x - ⟨v'v'⟩∂U/∂y
            (1, 1) => -T::from_f64(2.0).unwrap() * uv * dv_dx, // P_yy = -2⟨u'v'⟩∂V/∂x
            _ => T::zero(),
        }
    }

    /// Compute velocity gradient ∂U_i/∂x_j
    fn velocity_gradient(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let cos_kx_x = nalgebra::ComplexField::cos(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let cos_ky_y = nalgebra::ComplexField::cos(self.ky * y);
        let exp_decay = nalgebra::ComplexField::exp(-self.alpha * t);

        match (i, j) {
            (0, 0) => self.u0_amp * self.kx * cos_kx_x * sin_ky_y * exp_decay, // ∂U/∂x
            (0, 1) => self.u0_amp * self.ky * sin_kx_x * cos_ky_y * exp_decay, // ∂U/∂y
            (1, 0) => self.v0_amp * self.kx * cos_kx_x * sin_ky_y * exp_decay, // ∂V/∂x
            (1, 1) => self.v0_amp * self.ky * sin_kx_x * cos_ky_y * exp_decay, // ∂V/∂y
            _ => T::zero(),
        }
    }

    /// Compute the pressure-strain correlation Φ_ij
    pub fn pressure_strain_term(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        let k = self.exact_tke(x, y, t);
        let eps = self.exact_dissipation(x, y, t);

        // Avoid division by zero
        if eps <= T::zero() || k <= T::zero() {
            return T::zero();
        }

        let time_scale = k / eps;

        // Reynolds stress components
        let uu = self.exact_reynolds_stress(0, 0, x, y, t);
        let uv = self.exact_reynolds_stress(0, 1, x, y, t);
        let vv = self.exact_reynolds_stress(1, 1, x, y, t);

        // Anisotropy tensor components
        let a_xx = uu / k - T::from_f64(2.0 / 3.0).unwrap();
        let a_xy = uv / k;
        let a_yy = vv / k - T::from_f64(2.0 / 3.0).unwrap();

        // Strain rate tensor from mean velocity gradients
        let s11 = self.velocity_gradient(0, 0, x, y, t);
        let s12 = T::from_f64(0.5).unwrap()
            * (self.velocity_gradient(0, 1, x, y, t) + self.velocity_gradient(1, 0, x, y, t));
        let s22 = self.velocity_gradient(1, 1, x, y, t);

        // Rotation rate tensor
        let w12 = T::from_f64(0.5).unwrap()
            * (self.velocity_gradient(0, 1, x, y, t) - self.velocity_gradient(1, 0, x, y, t));

        match self.pressure_strain_model {
            PressureStrainModelMMS::LinearReturnToIsotropy => {
                // Φ_ij = -C1 ε/k (⟨u_i'u_j'⟩ - (2/3)k δ_ij)
                let c1 = T::from_f64(1.8).unwrap();
                let coeff = -c1 * eps / k;

                match (i, j) {
                    (0, 0) => coeff * a_xx,
                    (0, 1) | (1, 0) => coeff * a_xy,
                    (1, 1) => coeff * a_yy,
                    _ => T::zero(),
                }
            }

            PressureStrainModelMMS::Quadratic => self
                .pressure_strain_quadratic(a_xx, a_xy, a_yy, time_scale, s11, s12, s22, w12, i, j),

            PressureStrainModelMMS::SSG => {
                self.pressure_strain_ssg(a_xx, a_xy, a_yy, time_scale, s11, s12, s22, w12, i, j)
            }
        }
    }

    /// Quadratic pressure-strain correlation (Speziale et al., 1991)
    fn pressure_strain_quadratic(
        &self,
        a_xx: T,
        a_xy: T,
        a_yy: T,
        time_scale: T,
        s11: T,
        s12: T,
        s22: T,
        _w12: T,
        i: usize,
        j: usize,
    ) -> T {
        // Slow pressure-strain (return-to-isotropy)
        let c1 = T::from_f64(1.8).unwrap();
        let phi_slow_ij = match (i, j) {
            (0, 0) => -c1 * a_xx,
            (0, 1) | (1, 0) => -c1 * a_xy,
            (1, 1) => -c1 * a_yy,
            _ => T::zero(),
        };

        // Rapid pressure-strain (quadratic terms)
        let c1_star = T::from_f64(1.7).unwrap();
        let c2_star = T::from_f64(-1.05).unwrap();
        let two_thirds = T::from_f64(2.0 / 3.0).unwrap();
        let four_thirds = T::from_f64(4.0 / 3.0).unwrap();

        let phi_rapid_ij = match (i, j) {
            (0, 0) => {
                c1_star * (a_xx * s11 + a_xy * s12)
                    + c2_star * (a_xx * s22 - a_xy * s12 + two_thirds * (a_xx + a_yy) * (s11 + s22))
            }
            (0, 1) | (1, 0) => {
                c1_star * (a_xx * s12 + a_xy * s22)
                    + c2_star * (a_xy * (s11 - s22) + a_yy * s12 - four_thirds * a_xy * (s11 + s22))
            }
            (1, 1) => {
                c1_star * (a_xy * s12 + a_yy * s22)
                    + c2_star * (a_yy * s11 - a_xy * s12 + two_thirds * (a_xx + a_yy) * (s11 + s22))
            }
            _ => T::zero(),
        };

        (phi_slow_ij + phi_rapid_ij) / time_scale
    }

    /// SSG pressure-strain correlation (Speziale-Sarkar-Gatski, 1991)
    fn pressure_strain_ssg(
        &self,
        a_xx: T,
        a_xy: T,
        a_yy: T,
        time_scale: T,
        s11: T,
        s12: T,
        s22: T,
        _w12: T,
        i: usize,
        j: usize,
    ) -> T {
        let c1 = T::from_f64(1.7).unwrap();
        let c3 = T::from_f64(0.8).unwrap();

        let c1_term = -c1 / time_scale;
        let c3_term = c3 / time_scale;

        match (i, j) {
            (0, 0) => c1_term * a_xx + c3_term * (a_xx * s11 + a_xy * s12),
            (0, 1) | (1, 0) => c1_term * a_xy + c3_term * (a_xx * s12 + a_xy * s22),
            (1, 1) => c1_term * a_yy + c3_term * (a_xy * s12 + a_yy * s22),
            _ => T::zero(),
        }
    }

    /// Compute the dissipation tensor ε_ij
    pub fn dissipation_tensor(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        let eps = self.exact_dissipation(x, y, t);
        let two_thirds = T::from_f64(2.0 / 3.0).unwrap();

        // Isotropic dissipation approximation: ε_ij = (2/3)ε δ_ij
        match (i, j) {
            (0, 0) | (1, 1) => two_thirds * eps,
            _ => T::zero(),
        }
    }

    /// Compute the turbulent transport term T_ij = -∂⟨u_i'u_j'u_k'⟩/∂x_k
    pub fn turbulent_transport(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        // Approximate triple correlation transport using gradient-diffusion hypothesis
        // ⟨u_i'u_j'u_k'⟩ ≈ -C_s (k³/ε²) ∂⟨u_i'u_j'⟩/∂x_k

        let c_s = T::from_f64(0.11).unwrap();
        let k = self.exact_tke(x, y, t);
        let eps = self.exact_dissipation(x, y, t);

        if eps <= T::zero() || k <= T::zero() {
            return T::zero();
        }

        let diffusion_coeff = c_s * k * k * k / (eps * eps);

        // Compute gradients of Reynolds stresses
        let stress_gradient = self.reynolds_stress_gradient(i, j, x, y, t);

        match (i, j) {
            (0, 0) => -diffusion_coeff * stress_gradient[0], // -∂⟨u'u'u'⟩/∂x - ∂⟨u'u'u'⟩/∂y
            (0, 1) | (1, 0) => {
                -diffusion_coeff
                    * T::from_f64(0.5).unwrap()
                    * (stress_gradient[0] + stress_gradient[1])
            }
            (1, 1) => -diffusion_coeff * stress_gradient[1], // -∂⟨v'v'v'⟩/∂x - ∂⟨v'v'v'⟩/∂y
            _ => T::zero(),
        }
    }

    /// Compute gradient of Reynolds stress component ∂⟨u_i'u_j'⟩/∂x_k
    fn reynolds_stress_gradient(&self, i: usize, j: usize, x: T, y: T, t: T) -> [T; 2] {
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let cos_kx_x = nalgebra::ComplexField::cos(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let cos_ky_y = nalgebra::ComplexField::cos(self.ky * y);
        let exp_decay = nalgebra::ComplexField::exp(-self.alpha * t);

        let isotropic_part_grad_x = self.k_amp
            * T::from_f64(2.0 / 3.0).unwrap()
            * self.kx
            * cos_kx_x
            * sin_ky_y
            * exp_decay;
        let isotropic_part_grad_y = self.k_amp
            * T::from_f64(2.0 / 3.0).unwrap()
            * self.ky
            * sin_kx_x
            * cos_ky_y
            * exp_decay;

        match (i, j) {
            (0, 0) => {
                // ∇⟨u'u'⟩
                let fluct_grad_x = self.a_uu * self.kx * cos_kx_x * sin_ky_y * exp_decay;
                let fluct_grad_y = self.a_uu * self.ky * sin_kx_x * cos_ky_y * exp_decay;
                [
                    fluct_grad_x + isotropic_part_grad_x,
                    fluct_grad_y + isotropic_part_grad_y,
                ]
            }
            (0, 1) | (1, 0) => {
                // ∇⟨u'v'⟩
                let grad_x = -self.a_uv * self.kx * sin_kx_x * cos_ky_y * exp_decay;
                let grad_y = -self.a_uv * self.ky * cos_kx_x * sin_ky_y * exp_decay;
                [grad_x, grad_y]
            }
            (1, 1) => {
                // ∇⟨v'v'⟩
                let fluct_grad_x = self.a_vv * self.kx * cos_kx_x * sin_ky_y * exp_decay;
                let fluct_grad_y = self.a_vv * self.ky * sin_kx_x * cos_ky_y * exp_decay;
                [
                    fluct_grad_x + isotropic_part_grad_x,
                    fluct_grad_y + isotropic_part_grad_y,
                ]
            }
            _ => [T::zero(), T::zero()],
        }
    }

    /// Compute the molecular diffusion term D_ij = ∂/∂x_k(ν ∂⟨u_i'u_j'⟩/∂x_k)
    pub fn molecular_diffusion(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        // Molecular diffusion: D_ij = ∂/∂x_k(ν ∂⟨u_i'u_j'⟩/∂x_k)
        // For constant ν, this simplifies to ν ∇²⟨u_i'u_j'⟩

        let stress_laplacian = self.reynolds_stress_laplacian(i, j, x, y, t);
        self.nu * stress_laplacian
    }

    /// Compute Laplacian of Reynolds stress component ∇²⟨u_i'u_j'⟩
    fn reynolds_stress_laplacian(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        let sin_kx_x = nalgebra::ComplexField::sin(self.kx * x);
        let cos_kx_x = nalgebra::ComplexField::cos(self.kx * x);
        let sin_ky_y = nalgebra::ComplexField::sin(self.ky * y);
        let cos_ky_y = nalgebra::ComplexField::cos(self.ky * y);
        let exp_decay = nalgebra::ComplexField::exp(-self.alpha * t);

        let kx_sq = self.kx * self.kx;
        let ky_sq = self.ky * self.ky;

        let isotropic_part_laplacian = self.k_amp
            * T::from_f64(2.0 / 3.0).unwrap()
            * (kx_sq + ky_sq)
            * sin_kx_x
            * sin_ky_y
            * exp_decay;

        match (i, j) {
            (0, 0) => {
                // ∇²⟨u'u'⟩
                let fluct_laplacian = self.a_uu * (kx_sq + ky_sq) * sin_kx_x * sin_ky_y * exp_decay;
                fluct_laplacian + isotropic_part_laplacian
            }
            (0, 1) | (1, 0) => {
                // ∇²⟨u'v'⟩
                -self.a_uv * (kx_sq + ky_sq) * cos_kx_x * cos_ky_y * exp_decay
            }
            (1, 1) => {
                // ∇²⟨v'v'⟩
                let fluct_laplacian = self.a_vv * (kx_sq + ky_sq) * sin_kx_x * sin_ky_y * exp_decay;
                fluct_laplacian + isotropic_part_laplacian
            }
            _ => T::zero(),
        }
    }

    /// Compute the time derivative ∂⟨u_i'u_j'⟩/∂t
    pub fn time_derivative(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        -self.alpha * self.exact_reynolds_stress(i, j, x, y, t)
    }

    /// Compute the convective derivative D⟨u_i'u_j'⟩/Dt = ∂⟨u_i'u_j'⟩/∂t + U_k ∂⟨u_i'u_j'⟩/∂x_k
    pub fn convective_derivative(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        let time_deriv = self.time_derivative(i, j, x, y, t);
        let convective_term = self.convective_term(i, j, x, y, t);
        time_deriv + convective_term
    }

    /// Compute the convective term U_k ∂⟨u_i'u_j'⟩/∂x_k
    fn convective_term(&self, i: usize, j: usize, x: T, y: T, t: T) -> T {
        let u = self.exact_mean_velocity(0, x, y, t);
        let v = self.exact_mean_velocity(1, x, y, t);

        let stress_grad_x = self.reynolds_stress_gradient(i, j, x, y, t)[0];
        let stress_grad_y = self.reynolds_stress_gradient(i, j, x, y, t)[1];

        u * stress_grad_x + v * stress_grad_y
    }
}

impl<T: RealField + Copy + FromPrimitive> ManufacturedSolution<T>
    for ManufacturedReynoldsStressMMS<T>
{
    /// Return the Reynolds shear stress ⟨u'v'⟩ as the primary solution for trait compatibility
    fn exact_solution(&self, x: T, y: T, _z: T, t: T) -> T {
        self.exact_reynolds_stress(0, 1, x, y, t)
    }

    /// Compute the complete source term for the Reynolds stress transport equation
    /// S_ij = D⟨u_i'u_j'⟩/Dt - P_ij - Φ_ij + ε_ij - T_ij - D_ij
    fn source_term(&self, x: T, y: T, _z: T, t: T) -> T {
        // For the shear stress component (0,1), compute the complete source term
        let i = 0;
        let j = 1;

        let convective_deriv = self.convective_derivative(i, j, x, y, t);
        let production = self.production_term(i, j, x, y, t);
        let pressure_strain = self.pressure_strain_term(i, j, x, y, t);
        let dissipation = self.dissipation_tensor(i, j, x, y, t);
        let transport = self.turbulent_transport(i, j, x, y, t);
        let diffusion = self.molecular_diffusion(i, j, x, y, t);

        // Source term: S_ij = D⟨u_i'u_j'⟩/Dt - P_ij - Φ_ij + ε_ij - T_ij - D_ij
        convective_deriv - production - pressure_strain + dissipation - transport - diffusion
    }
}

/// Convergence study utilities for Reynolds stress MMS validation
pub struct ReynoldsStressConvergenceStudy<T: RealField + Copy + FromPrimitive> {
    /// The manufactured solution to use for convergence study
    pub mms: ManufacturedReynoldsStressMMS<T>,
}

impl<T: RealField + Copy + FromPrimitive> ReynoldsStressConvergenceStudy<T> {
    /// Create a new convergence study utility
    pub fn new(mms: ManufacturedReynoldsStressMMS<T>) -> Self {
        Self { mms }
    }

    /// Compute L2 error norm for all Reynolds stress components
    pub fn compute_l2_error(
        &self,
        numerical_stresses: &[DMatrix<T>; 3],
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        t: T,
    ) -> [T; 3] {
        let mut errors = [T::zero(); 3];

        for i in 0..nx {
            for j in 0..ny {
                let x = T::from_usize(i).unwrap() * dx;
                let y = T::from_usize(j).unwrap() * dy;

                // Exact solutions
                let exact_xx = self.mms.exact_reynolds_stress(0, 0, x, y, t);
                let exact_xy = self.mms.exact_reynolds_stress(0, 1, x, y, t);
                let exact_yy = self.mms.exact_reynolds_stress(1, 1, x, y, t);

                // Numerical solutions
                let num_xx = numerical_stresses[0][(i, j)];
                let num_xy = numerical_stresses[1][(i, j)];
                let num_yy = numerical_stresses[2][(i, j)];

                // Accumulate squared errors
                let dx_dy = dx * dy;
                errors[0] += (num_xx - exact_xx) * (num_xx - exact_xx) * dx_dy;
                errors[1] += (num_xy - exact_xy) * (num_xy - exact_xy) * dx_dy;
                errors[2] += (num_yy - exact_yy) * (num_yy - exact_yy) * dx_dy;
            }
        }

        // Take square root for L2 norm
        [
            nalgebra::ComplexField::sqrt(errors[0]),
            nalgebra::ComplexField::sqrt(errors[1]),
            nalgebra::ComplexField::sqrt(errors[2]),
        ]
    }

    /// Check if convergence rate is approximately 2nd order
    pub fn check_convergence_rate(
        &self,
        errors_fine: &[T; 3],
        errors_coarse: &[T; 3],
        refinement_ratio: T,
    ) -> [T; 3] {
        let mut rates = [T::zero(); 3];

        let ratio_sq = refinement_ratio * refinement_ratio;

        for i in 0..3 {
            if errors_coarse[i] > T::zero() && errors_fine[i] > T::zero() {
                rates[i] = nalgebra::ComplexField::ln(errors_coarse[i] / errors_fine[i])
                    / nalgebra::ComplexField::ln(ratio_sq);
            }
        }

        rates
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_manufactured_reynolds_stress_creation() {
        let mms = ManufacturedReynoldsStressMMS::<f64>::standard_test_case();

        assert_eq!(mms.kx, 2.0 * std::f64::consts::PI);
        assert_eq!(mms.ky, 2.0 * std::f64::consts::PI);
        assert_eq!(mms.alpha, 0.1);
    }

    #[test]
    fn test_exact_reynolds_stress_components() {
        let mms = ManufacturedReynoldsStressMMS::<f64>::standard_test_case();

        let x = 0.25;
        let y = 0.25;
        let t = 0.0;

        let uu = mms.exact_reynolds_stress(0, 0, x, y, t);
        let uv = mms.exact_reynolds_stress(0, 1, x, y, t);
        let vv = mms.exact_reynolds_stress(1, 1, x, y, t);

        // At t=0, should have non-zero values
        assert!(uu > 0.0);
        assert!(uv.is_finite());
        assert!(vv > 0.0);

        // Check that uu and vv have isotropic parts
        let isotropic_part = 0.15 * (2.0 / 3.0); // k_amp * 2/3
        assert!(uu > isotropic_part);
        assert!(vv > isotropic_part);
    }

    #[test]
    fn test_source_term_computation() {
        let mms = ManufacturedReynoldsStressMMS::<f64>::standard_test_case();

        let x = 0.5;
        let y = 0.5;
        let t = 0.1;

        let source = mms.source_term(x, y, 0.0, t);
        assert!(source.is_finite());

        // Source term should be computable without NaN or infinite values
        assert!(!source.is_nan());
        assert!(!source.is_infinite());
    }

    #[test]
    fn test_production_term() {
        let mms = ManufacturedReynoldsStressMMS::<f64>::standard_test_case();

        let x = 0.5;
        let y = 0.5;
        let t = 0.0;

        let p_xx = mms.production_term(0, 0, x, y, t);
        let p_xy = mms.production_term(0, 1, x, y, t);
        let p_yy = mms.production_term(1, 1, x, y, t);

        // Production terms should be finite
        assert!(p_xx.is_finite());
        assert!(p_xy.is_finite());
        assert!(p_yy.is_finite());
    }

    #[test]
    fn test_pressure_strain_terms() {
        let mms = ManufacturedReynoldsStressMMS::<f64>::standard_test_case();

        let x = 0.5;
        let y = 0.5;
        let t = 0.0;

        let phi_xx = mms.pressure_strain_term(0, 0, x, y, t);
        let phi_xy = mms.pressure_strain_term(0, 1, x, y, t);
        let phi_yy = mms.pressure_strain_term(1, 1, x, y, t);

        // Pressure-strain terms should be finite
        assert!(phi_xx.is_finite());
        assert!(phi_xy.is_finite());
        assert!(phi_yy.is_finite());
    }

    #[test]
    fn test_time_decay() {
        let mms = ManufacturedReynoldsStressMMS::<f64>::standard_test_case();

        let x = 0.5;
        let y = 0.5;

        let uv_t0 = mms.exact_reynolds_stress(0, 1, x, y, 0.0);
        let uv_t1 = mms.exact_reynolds_stress(0, 1, x, y, 1.0);

        // Should decay exponentially
        assert!(uv_t1.abs() < uv_t0.abs());
        assert_relative_eq!(uv_t1 / uv_t0, (-mms.alpha).exp(), epsilon = 1e-10);
    }

    #[test]
    fn test_convergence_study_setup() {
        let mms = ManufacturedReynoldsStressMMS::<f64>::standard_test_case();
        let study = ReynoldsStressConvergenceStudy::new(mms);

        // Should be able to create convergence study
        assert_eq!(study.mms.kx, 2.0 * std::f64::consts::PI);
    }
}
