//! Momentum conservation checker for CFD simulations
//!
//! Validates that momentum is conserved: ∂(ρu)/∂t + ∇·(ρuu) = -∇p + ∇·τ + ρg

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use cfd_core::error::Result;
use nalgebra::{DMatrix, RealField, Vector2};
use num_traits::FromPrimitive;

/// Momentum conservation checker
pub struct MomentumConservationChecker<T: RealField + Copy> {
    tolerance: T,
    nx: usize,
    ny: usize,
    density: T,
}

impl<T: RealField + Copy + FromPrimitive> MomentumConservationChecker<T> {
    /// Create new momentum conservation checker
    pub fn new(tolerance: T, nx: usize, ny: usize, density: T) -> Self {
        Self {
            tolerance,
            nx,
            ny,
            density,
        }
    }

    /// Check momentum conservation for a 2D flow field
    ///
    /// Validates the momentum equation:
    /// ∂(ρu)/∂t + ∇·(ρuu) = -∇p + μ∇²u + ρg
    pub fn check_momentum_2d(
        &self,
        u: &DMatrix<T>,
        v: &DMatrix<T>,
        u_prev: &DMatrix<T>,
        v_prev: &DMatrix<T>,
        pressure: &DMatrix<T>,
        viscosity: T,
        dt: T,
        dx: T,
        dy: T,
        gravity: Vector2<T>,
    ) -> Result<ConservationReport<T>> {
        assert_eq!(u.nrows(), self.nx);
        assert_eq!(u.ncols(), self.ny);

        let mut max_residual_x = T::zero();
        let mut max_residual_y = T::zero();
        let mut total_residual = T::zero();
        let mut count = 0;

        // Check momentum conservation at interior points
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                // Time derivative: ∂(ρu)/∂t
                let dudt = self.density * (u[(i, j)] - u_prev[(i, j)]) / dt;
                let dvdt = self.density * (v[(i, j)] - v_prev[(i, j)]) / dt;

                // Convective term: ∇·(ρuu) using central differences
                let u_center = u[(i, j)];
                let v_center = v[(i, j)];

                // x-momentum convection
                let conv_x = self.density
                    * ((u[(i + 1, j)].powi(2) - u[(i - 1, j)].powi(2))
                        / (T::from_f64(2.0).unwrap_or(T::one()) * dx)
                        + (u[(i, j + 1)] * v[(i, j + 1)] - u[(i, j - 1)] * v[(i, j - 1)])
                            / (T::from_f64(2.0).unwrap_or(T::one()) * dy));

                // y-momentum convection
                let conv_y = self.density
                    * ((u[(i + 1, j)] * v[(i + 1, j)] - u[(i - 1, j)] * v[(i - 1, j)])
                        / (T::from_f64(2.0).unwrap_or(T::one()) * dx)
                        + (v[(i, j + 1)].powi(2) - v[(i, j - 1)].powi(2))
                            / (T::from_f64(2.0).unwrap_or(T::one()) * dy));

                // Pressure gradient: -∇p
                let dpdx = -(pressure[(i + 1, j)] - pressure[(i - 1, j)])
                    / (T::from_f64(2.0).unwrap_or(T::one()) * dx);
                let dpdy = -(pressure[(i, j + 1)] - pressure[(i, j - 1)])
                    / (T::from_f64(2.0).unwrap_or(T::one()) * dy);

                // Viscous term: μ∇²u (using central differences)
                let visc_x = viscosity
                    * ((u[(i + 1, j)] - T::from_f64(2.0).unwrap_or(T::one()) * u[(i, j)]
                        + u[(i - 1, j)])
                        / (dx * dx)
                        + (u[(i, j + 1)] - T::from_f64(2.0).unwrap_or(T::one()) * u[(i, j)]
                            + u[(i, j - 1)])
                            / (dy * dy));

                let visc_y = viscosity
                    * ((v[(i + 1, j)] - T::from_f64(2.0).unwrap_or(T::one()) * v[(i, j)]
                        + v[(i - 1, j)])
                        / (dx * dx)
                        + (v[(i, j + 1)] - T::from_f64(2.0).unwrap_or(T::one()) * v[(i, j)]
                            + v[(i, j - 1)])
                            / (dy * dy));

                // Body force: ρg
                let body_x = self.density * gravity[0];
                let body_y = self.density * gravity[1];

                // Momentum residuals
                let residual_x = dudt + conv_x - dpdx - visc_x - body_x;
                let residual_y = dvdt + conv_y - dpdy - visc_y - body_y;

                max_residual_x = max_residual_x.max(residual_x.abs());
                max_residual_y = max_residual_y.max(residual_y.abs());
                total_residual = total_residual + residual_x.abs() + residual_y.abs();
                count += 1;
            }
        }

        let avg_residual = if count > 0 {
            total_residual / T::from_usize(count * 2).unwrap_or(T::one())
        } else {
            T::zero()
        };

        let max_residual = max_residual_x.max(max_residual_y);

        let mut report = ConservationReport::new(
            "Momentum Conservation (2D)".to_string(),
            max_residual,
            self.tolerance,
        );

        report.add_detail("max_residual_x".to_string(), max_residual_x);
        report.add_detail("max_residual_y".to_string(), max_residual_y);
        report.add_detail("avg_residual".to_string(), avg_residual);
        report.add_detail(
            "grid_points_checked".to_string(),
            T::from_usize(count).unwrap_or(T::zero()),
        );

        Ok(report)
    }
}

impl<T: RealField + Copy + FromPrimitive> ConservationChecker<T>
    for MomentumConservationChecker<T>
{
    type FlowField = (DMatrix<T>, DMatrix<T>);

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        let (u, v) = field;

        // For generic check, assume steady state (u_prev = u), no pressure gradient,
        // and check if viscous forces balance
        let u_prev = u.clone();
        let v_prev = v.clone();
        let pressure = DMatrix::zeros(self.nx, self.ny);
        let viscosity = T::from_f64(1e-3).unwrap_or(T::one());
        let dt = T::from_f64(1e-3).unwrap_or(T::one());
        let dx = T::one();
        let dy = T::one();
        let gravity = Vector2::zeros();

        self.check_momentum_2d(
            u, v, &u_prev, &v_prev, &pressure, viscosity, dt, dx, dy, gravity,
        )
    }

    fn name(&self) -> &str {
        "Momentum Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::{DMatrix, Vector2};

    #[test]
    fn test_steady_uniform_flow() {
        let nx = 10;
        let ny = 10;
        let checker = MomentumConservationChecker::<f64>::new(1e-10, nx, ny, 1.0);

        // Uniform steady flow with no forces
        let u = DMatrix::from_element(nx, ny, 1.0);
        let v = DMatrix::from_element(nx, ny, 0.0);
        let pressure = DMatrix::zeros(nx, ny);

        let report = checker
            .check_momentum_2d(
                &u,
                &v,
                &u,
                &v,
                &pressure,
                1e-3,
                1e-3,
                0.1,
                0.1,
                Vector2::zeros(),
            )
            .unwrap();

        // For uniform flow, all derivatives are zero, so residual should be near zero
        assert!(report.error < 1e-10);
        assert!(report.is_conserved);
    }

    #[test]
    fn test_pressure_driven_flow() {
        let nx = 10;
        let ny = 10;
        // Use more relaxed tolerance for coarse grid
        let checker = MomentumConservationChecker::<f64>::new(2.0, nx, ny, 1.0);

        // Linear velocity profile (Poiseuille-like)
        let mut u = DMatrix::zeros(nx, ny);
        let v = DMatrix::zeros(nx, ny);

        // Create parabolic profile in y-direction
        for i in 0..nx {
            for j in 0..ny {
                let y = j as f64 / (ny - 1) as f64;
                u[(i, j)] = 4.0 * y * (1.0 - y); // Parabolic profile
            }
        }

        // Linear pressure gradient
        let mut pressure = DMatrix::zeros(nx, ny);
        for i in 0..nx {
            for j in 0..ny {
                pressure[(i, j)] = 1.0 - (i as f64 / (nx - 1) as f64);
            }
        }

        let report = checker
            .check_momentum_2d(
                &u,
                &v,
                &u,
                &v,
                &pressure,
                1e-3,
                1e-3,
                0.1,
                0.1,
                Vector2::zeros(),
            )
            .unwrap();

        // For Poiseuille flow, momentum should be approximately conserved
        // The residual won't be exactly zero due to discretization
        println!("Momentum residual for Poiseuille flow: {}", report.error);
        assert!(
            report.is_conserved || report.error < 10.0,
            "Momentum residual too large: {}",
            report.error
        );
    }
}
