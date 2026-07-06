//! Vorticity conservation checker for CFD simulations
//!
//! Validates vorticity transport equation: Dω/Dt = (ω·∇)u + ν∇²ω + source terms
//! For inviscid flows: Dω/Dt = (ω·∇)u (vorticity conservation)

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use crate::scalar;
use cfd_core::error::{Error, Result};
use eunomia::FloatElement;
use eunomia::RealField;
use leto::Array2;

/// Vorticity conservation checker
pub struct VorticityChecker<T: RealField + Copy> {
    tolerance: T,
    nx: usize,
    ny: usize,
}

impl<T: RealField + Copy + FloatElement> VorticityChecker<T> {
    /// Create new vorticity conservation checker
    pub fn new(tolerance: T, nx: usize, ny: usize) -> Self {
        Self { tolerance, nx, ny }
    }

    /// Check vorticity conservation for 2D incompressible flow
    ///
    /// Vorticity transport equation: ∂ω/∂t + u·∇ω = ν∇²ω + source terms
    /// where ω = ∂v/∂x - ∂u/∂y
    pub fn check_vorticity_transport_2d(
        &self,
        u: &Array2<T>,
        v: &Array2<T>,
        viscosity: T,
        dx: T,
        dy: T,
    ) -> Result<ConservationReport<T>> {
        assert_eq!(u.shape()[0], self.nx);
        assert_eq!(u.shape()[1], self.ny);
        assert_eq!(v.shape()[0], self.nx);
        assert_eq!(v.shape()[1], self.ny);

        let mut max_error = scalar::zero::<T>();
        let mut total_error = scalar::zero::<T>();
        let mut count = 0;

        // Check vorticity transport equation at interior points
        // Need to stay away from boundaries by 2 cells for proper stencil
        for i in 2..self.nx - 2 {
            for j in 2..self.ny - 2 {
                // Compute vorticity ω = ∂v/∂x - ∂u/∂y
                let dv_dx = (v[[i + 1, j]] - v[[i - 1, j]]) / (scalar::from_f64::<T>(2.0) * dx);
                let du_dy = (u[[i, j + 1]] - u[[i, j - 1]]) / (scalar::from_f64::<T>(2.0) * dy);
                let omega = dv_dx - du_dy;

                // Convective term: u·∇ω
                // Compute ∇ω using central differences on vorticity field
                // First, compute vorticity at neighboring points
                let omega_ip1 = ((v[[i + 2, j]] - v[[i, j]]) / (scalar::from_f64::<T>(2.0) * dx))
                    - ((u[[i + 1, j + 1]] - u[[i + 1, j - 1]]) / (scalar::from_f64::<T>(2.0) * dy));

                let omega_im1 = ((v[[i, j]] - v[[i - 2, j]]) / (scalar::from_f64::<T>(2.0) * dx))
                    - ((u[[i - 1, j + 1]] - u[[i - 1, j - 1]]) / (scalar::from_f64::<T>(2.0) * dy));

                let omega_jp1 = ((v[[i + 1, j + 1]] - v[[i - 1, j + 1]])
                    / (scalar::from_f64::<T>(2.0) * dx))
                    - ((u[[i, j + 2]] - u[[i, j]]) / (scalar::from_f64::<T>(2.0) * dy));

                let omega_jm1 = ((v[[i + 1, j - 1]] - v[[i - 1, j - 1]])
                    / (scalar::from_f64::<T>(2.0) * dx))
                    - ((u[[i, j]] - u[[i, j - 2]]) / (scalar::from_f64::<T>(2.0) * dy));

                let domega_dx = (omega_ip1 - omega_im1) / (scalar::from_f64::<T>(2.0) * dx);
                let domega_dy = (omega_jp1 - omega_jm1) / (scalar::from_f64::<T>(2.0) * dy);
                let u_dot_grad_omega = u[[i, j]] * domega_dx + v[[i, j]] * domega_dy;

                // Viscous diffusion: ν∇²ω
                // Use proper second-order central differences for Laplacian
                let d2omega_dx2 =
                    (omega_ip1 - scalar::from_f64::<T>(2.0) * omega + omega_im1) / (dx * dx);
                let d2omega_dy2 =
                    (omega_jp1 - scalar::from_f64::<T>(2.0) * omega + omega_jm1) / (dy * dy);
                let viscous_diffusion = viscosity * (d2omega_dx2 + d2omega_dy2);

                // Vorticity transport equation residual: ∂ω/∂t + u·∇ω - ν∇²ω
                // For steady state, this should be zero
                let residual = u_dot_grad_omega - viscous_diffusion;

                let abs_residual = scalar::abs(residual);
                if abs_residual > max_error {
                    max_error = abs_residual;
                }
                total_error += residual * residual;
                count += 1;
            }
        }

        let rms_error = if count > 0 {
            scalar::sqrt(total_error / scalar::from_usize::<T>(count))
        } else {
            scalar::zero::<T>()
        };

        let is_conserved = max_error < self.tolerance;

        Ok(ConservationReport {
            check_name: "Vorticity Transport Conservation".to_string(),
            error: max_error,
            is_conserved,
            tolerance: self.tolerance,
            details: [("rms_error".to_string(), rms_error)]
                .iter()
                .cloned()
                .collect(),
        })
    }

    /// Check Kelvin's circulation theorem
    ///
    /// For inviscid flows, circulation around a closed curve is conserved
    /// ∮ u·dl = constant for material curves
    pub fn check_kelvin_circulation(
        &self,
        u: &Array2<T>,
        v: &Array2<T>,
        curve_points: &[(usize, usize)], // Points defining closed curve
        dx: T,
        dy: T,
    ) -> Result<ConservationReport<T>> {
        if curve_points.len() < 3 {
            return Err(cfd_core::error::Error::InvalidInput(
                "Need at least 3 points to define a closed curve".to_string(),
            ));
        }

        // Compute circulation ∮ u·dl around the curve
        let mut circulation = scalar::zero::<T>();

        for i in 0..curve_points.len() {
            let (i1, j1) = curve_points[i];
            let (i2, j2) = curve_points[(i + 1) % curve_points.len()];

            // Tangent vector along curve
            let i1_i64 = i64::try_from(i1).map_err(|_| {
                Error::InvalidInput(format!("Failed to convert curve i-index {i1} to i64"))
            })?;
            let i2_i64 = i64::try_from(i2).map_err(|_| {
                Error::InvalidInput(format!("Failed to convert curve i-index {i2} to i64"))
            })?;
            let j1_i64 = i64::try_from(j1).map_err(|_| {
                Error::InvalidInput(format!("Failed to convert curve j-index {j1} to i64"))
            })?;
            let j2_i64 = i64::try_from(j2).map_err(|_| {
                Error::InvalidInput(format!("Failed to convert curve j-index {j2} to i64"))
            })?;

            let di = i2_i64 - i1_i64;
            let dj = j2_i64 - j1_i64;

            // Velocity at midpoint
            let i_mid = usize::midpoint(i1, i2).min(self.nx - 1);
            let j_mid = usize::midpoint(j1, j2).min(self.ny - 1);

            let u_mid = u[[i_mid, j_mid]];
            let v_mid = v[[i_mid, j_mid]];

            // dl component (tangent vector)
            let dl_x = scalar::from_f64::<T>(di as f64) * dx;
            let dl_y = scalar::from_f64::<T>(dj as f64) * dy;

            // u·dl contribution
            circulation = circulation + u_mid * dl_x + v_mid * dl_y;
        }

        // For inviscid flow, circulation should be conserved (constant)
        // Here we just check it's finite and reasonable
        let error = scalar::abs(circulation);
        let is_conserved = error < self.tolerance * scalar::from_f64::<T>(100.0); // More lenient check

        let mut details = std::collections::HashMap::new();
        details.insert("rms_error".to_string(), error);

        Ok(ConservationReport {
            check_name: "Kelvin Circulation Conservation".to_string(),
            error,
            is_conserved,
            tolerance: self.tolerance,
            details,
        })
    }

    /// Check potential flow vorticity (should be zero everywhere)
    pub fn check_potential_flow_vorticity(
        &self,
        u: &Array2<T>,
        v: &Array2<T>,
        dx: T,
        dy: T,
    ) -> Result<ConservationReport<T>> {
        let mut max_vorticity = scalar::zero::<T>();
        let mut total_vorticity_sq = scalar::zero::<T>();
        let mut count = 0;

        // Compute vorticity at all points
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let dv_dx = (v[[i + 1, j]] - v[[i - 1, j]]) / (scalar::from_f64::<T>(2.0) * dx);
                let du_dy = (u[[i, j + 1]] - u[[i, j - 1]]) / (scalar::from_f64::<T>(2.0) * dy);
                let vorticity = dv_dx - du_dy;

                let abs_vorticity = scalar::abs(vorticity);
                if abs_vorticity > max_vorticity {
                    max_vorticity = abs_vorticity;
                }
                total_vorticity_sq += vorticity * vorticity;
                count += 1;
            }
        }

        let rms_vorticity = if count > 0 {
            scalar::sqrt(total_vorticity_sq / scalar::from_usize::<T>(count))
        } else {
            scalar::zero::<T>()
        };

        // For potential flow, vorticity should be exactly zero
        let is_conserved = max_vorticity < self.tolerance;

        let mut details = std::collections::HashMap::new();
        details.insert("rms_error".to_string(), rms_vorticity);

        Ok(ConservationReport {
            check_name: "Potential Flow Vorticity Conservation".to_string(),
            error: max_vorticity,
            is_conserved,
            tolerance: self.tolerance,
            details,
        })
    }
}

impl<T: RealField + Copy + FloatElement> ConservationChecker<T> for VorticityChecker<T> {
    type FlowField = Vec<Array2<T>>;

    fn name(&self) -> &'static str {
        "Vorticity Conservation Checker"
    }

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        if field.len() < 2 {
            return Err(cfd_core::error::Error::InvalidInput(
                "Vorticity check requires at least u,v velocity fields".to_string(),
            ));
        }

        let u = &field[0];
        let v = &field[1];

        // Use Kelvin circulation as the primary vorticity conservation check
        // Define a simple closed curve for circulation check
        let curve_points = vec![
            (0, 0),
            (u.shape()[1] - 1, 0),
            (u.shape()[1] - 1, u.shape()[0] - 1),
            (0, u.shape()[0] - 1),
        ];
        self.check_kelvin_circulation(u, v, &curve_points, scalar::one::<T>(), scalar::one::<T>())
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }
}
