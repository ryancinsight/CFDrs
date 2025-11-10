//! Angular momentum conservation checker for CFD simulations
//!
//! Validates angular momentum conservation: ∂(ρr×u)/∂t + ∇·(ρr×u⊗u + pI × r) = ∇·τ × r
//! For incompressible flows: ∇·(r × u) = 0 (angular momentum conservation)

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use cfd_core::error::Result;
use cfd_core::conversion::SafeFromF64;
use nalgebra::{DMatrix, RealField, Vector2};
use num_traits::FromPrimitive;

/// Angular momentum conservation checker
pub struct AngularMomentumChecker<T: RealField + Copy> {
    tolerance: T,
    nx: usize,
    ny: usize,
    center: Vector2<T>, // Center of rotation
}

impl<T: RealField + Copy + FromPrimitive> AngularMomentumChecker<T> {
    /// Create new angular momentum conservation checker
    pub fn new(tolerance: T, nx: usize, ny: usize, center: Vector2<T>) -> Self {
        Self { tolerance, nx, ny, center }
    }

    /// Create checker with origin as center
    pub fn new_centered(tolerance: T, nx: usize, ny: usize) -> Self {
        let center = Vector2::new(T::zero(), T::zero());
        Self::new(tolerance, nx, ny, center)
    }

    /// Check angular momentum conservation for 2D incompressible flow
    ///
    /// For incompressible flows, angular momentum conservation requires:
    /// ∇·(r × u) = 0, where r is position vector from center of rotation
    pub fn check_angular_momentum_2d(
        &self,
        u: &DMatrix<T>,
        v: &DMatrix<T>,
        dx: T,
        dy: T,
    ) -> Result<ConservationReport<T>> {
        assert_eq!(u.nrows(), self.nx);
        assert_eq!(u.ncols(), self.ny);
        assert_eq!(v.nrows(), self.nx);
        assert_eq!(v.ncols(), self.ny);

        let mut max_error = T::zero();
        let mut total_error = T::zero();
        let mut count = 0;

        // Check angular momentum conservation at interior points
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                // Position vector from center
                let x = T::from_usize(i).unwrap() * dx - self.center.x;
                let y = T::from_usize(j).unwrap() * dy - self.center.y;
                let r_squared = x * x + y * y;

                if r_squared > T::zero() {
                    // Angular momentum density: r × u = r_x * v - r_y * u
                    // Divergence: ∇·(r × u) = ∂(r_y * v)/∂x + ∂(-r_x * u)/∂y

                    // ∂(r_y * v)/∂x using central difference
                    let rv_right = y * v[(i + 1, j)];
                    let rv_left = y * v[(i - 1, j)];
                    let drv_dx = (rv_right - rv_left) / (<T as SafeFromF64>::from_f64_or_one(2.0) * dx);

                    // ∂(-r_x * u)/∂y using central difference
                    let minus_ru_top = -x * u[(i, j + 1)];
                    let minus_ru_bottom = -x * u[(i, j - 1)];
                    let d_minus_ru_dy = (minus_ru_top - minus_ru_bottom) / (<T as SafeFromF64>::from_f64_or_one(2.0) * dy);

                    // Total divergence
                    let divergence = drv_dx + d_minus_ru_dy;

                    let abs_divergence = divergence.abs();
                    if abs_divergence > max_error {
                        max_error = abs_divergence;
                    }
                    total_error = total_error + divergence * divergence;
                    count += 1;
                }
            }
        }

        let rms_error = if count > 0 {
            (total_error / T::from_usize(count).unwrap()).sqrt()
        } else {
            T::zero()
        };

        let is_conserved = max_error < self.tolerance;

        let mut details = std::collections::HashMap::new();
        details.insert("rms_error".to_string(), rms_error);

        Ok(ConservationReport {
            check_name: "Angular Momentum Conservation".to_string(),
            error: max_error,
            is_conserved,
            tolerance: self.tolerance,
            details,
        })
    }

    /// Check angular momentum conservation for axisymmetric flow
    ///
    /// In cylindrical coordinates: ∂(ρr²ω)/∂t + ∇·(ρr²ω u) = ∇·(μ r² ∂ω/∂z) + boundary terms
    pub fn check_axisymmetric_angular_momentum(
        &self,
        u_r: &DMatrix<T>, // Radial velocity
        u_z: &DMatrix<T>, // Axial velocity
        omega: &DMatrix<T>, // Angular velocity (vorticity in θ direction)
        viscosity: T,
        density: T,
        dr: T,
        dz: T,
    ) -> Result<ConservationReport<T>> {
        assert_eq!(u_r.nrows(), self.nx);
        assert_eq!(u_r.ncols(), self.ny);
        assert_eq!(u_z.nrows(), self.nx);
        assert_eq!(u_z.ncols(), self.ny);
        assert_eq!(omega.nrows(), self.nx);
        assert_eq!(omega.ncols(), self.ny);

        let mut max_error = T::zero();
        let mut total_error = T::zero();
        let mut count = 0;

        // Check angular momentum conservation at interior points
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let r = T::from_usize(i).unwrap() * dr;

                if r > T::zero() {
                    // Angular momentum per unit mass: r²ω
                    // Convective term: ∇·(r²ω u_r) in r-direction + ∇·(r²ω u_z) in z-direction

                    // ∂(r²ω u_r)/∂r term
                    let r2_omega = r * r * omega[(i, j)];
                    let conv_r_right = r2_omega * u_r[(i + 1, j)];
                    let conv_r_left = r2_omega * u_r[(i - 1, j)];
                    let d_conv_r_dr = (conv_r_right - conv_r_left) / (<T as SafeFromF64>::from_f64_or_one(2.0) * dr);

                    // ∂(r²ω u_z)/∂z term
                    let conv_z_top = r2_omega * u_z[(i, j + 1)];
                    let conv_z_bottom = r2_omega * u_z[(i, j - 1)];
                    let d_conv_z_dz = (conv_z_top - conv_z_bottom) / (<T as SafeFromF64>::from_f64_or_one(2.0) * dz);

                    // Viscous diffusion: ∇·(μ r² ∂ω/∂z) ≈ μ r² ∂²ω/∂z²
                    let d2omega_dz2 = (omega[(i, j + 1)] - <T as SafeFromF64>::from_f64_or_one(2.0) * omega[(i, j)] + omega[(i, j - 1)]) / (dz * dz);
                    let viscous = viscosity * r * r * d2omega_dz2;

                    // Angular momentum equation residual
                    let residual = d_conv_r_dr + d_conv_z_dz - viscous;

                    let abs_residual = residual.abs();
                    if abs_residual > max_error {
                        max_error = abs_residual;
                    }
                    total_error = total_error + residual * residual;
                    count += 1;
                }
            }
        }

        let rms_error = if count > 0 {
            (total_error / T::from_usize(count).unwrap()).sqrt()
        } else {
            T::zero()
        };

        let is_conserved = max_error < self.tolerance;

        let mut details = std::collections::HashMap::new();
        details.insert("rms_error".to_string(), rms_error);

        Ok(ConservationReport {
            check_name: "Axisymmetric Angular Momentum Conservation".to_string(),
            error: max_error,
            is_conserved,
            tolerance: self.tolerance,
            details,
        })
    }
}

impl<T: RealField + Copy + FromPrimitive> ConservationChecker<T> for AngularMomentumChecker<T> {
    type FlowField = Vec<DMatrix<T>>;

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        if field.len() < 2 {
            return Err(cfd_core::error::Error::InvalidInput(
                "Angular momentum check requires at least u,v velocity fields".to_string()
            ));
        }

        let u = &field[0];
        let v = &field[1];

        // Use 2D check as primary result
        self.check_angular_momentum_2d(u, v, <T as SafeFromF64>::from_f64_or_one(1.0), <T as SafeFromF64>::from_f64_or_one(1.0))
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }

    fn name(&self) -> &str {
        "Angular Momentum Conservation Checker"
    }
}

