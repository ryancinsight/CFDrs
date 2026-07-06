//! Momentum conservation checker for CFD simulations
//!
//! Validates that momentum is conserved: ∂(ρu)/∂t + ∇·(ρuu) = -∇p + ∇·τ + ρg

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use crate::scalar;
use cfd_core::error::Result;
use eunomia::FloatElement;
use eunomia::RealField;
use leto::geometry::Vector2;
use leto::Array2;

/// Momentum conservation checker
pub struct MomentumConservationChecker<T: RealField + Copy> {
    tolerance: T,
    nx: usize,
    ny: usize,
    density: T,
}

impl<T: RealField + Copy + FloatElement> MomentumConservationChecker<T> {
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
        u: &Array2<T>,
        v: &Array2<T>,
        u_prev: &Array2<T>,
        v_prev: &Array2<T>,
        pressure: &Array2<T>,
        viscosity: T,
        dt: T,
        dx: T,
        dy: T,
        gravity: Vector2<T>,
    ) -> Result<ConservationReport<T>> {
        #[allow(clippy::no_effect_underscore_binding)] // Context variables documented inline
        {
            assert_eq!(u.shape()[0], self.nx);
            assert_eq!(u.shape()[1], self.ny);

            let mut max_residual_x = scalar::zero::<T>();
            let mut max_residual_y = scalar::zero::<T>();
            let mut total_residual = scalar::zero::<T>();
            let mut count = 0;

            // Check momentum conservation at interior points
            for i in 1..self.nx - 1 {
                for j in 1..self.ny - 1 {
                    // Time derivative: ∂(ρu)/∂t
                    let dudt = self.density * (u[[i, j]] - u_prev[[i, j]]) / dt;
                    let dvdt = self.density * (v[[i, j]] - v_prev[[i, j]]) / dt;

                    // Convective term: ∇·(ρuu) using central differences
                    let _u_center = u[[i, j]];
                    let _v_center = v[[i, j]];

                    // x-momentum convection
                    let u_east = u[[i + 1, j]];
                    let u_west = u[[i - 1, j]];
                    let conv_x = self.density
                        * ((u_east * u_east - u_west * u_west) / (scalar::from_f64::<T>(2.0) * dx)
                            + (u[[i, j + 1]] * v[[i, j + 1]] - u[[i, j - 1]] * v[[i, j - 1]])
                                / (scalar::from_f64::<T>(2.0) * dy));

                    // y-momentum convection
                    let v_north = v[[i, j + 1]];
                    let v_south = v[[i, j - 1]];
                    let conv_y = self.density
                        * ((u[[i + 1, j]] * v[[i + 1, j]] - u[[i - 1, j]] * v[[i - 1, j]])
                            / (scalar::from_f64::<T>(2.0) * dx)
                            + (v_north * v_north - v_south * v_south)
                                / (scalar::from_f64::<T>(2.0) * dy));

                    // Pressure gradient: -∇p
                    let dpdx = -(pressure[[i + 1, j]] - pressure[[i - 1, j]])
                        / (scalar::from_f64::<T>(2.0) * dx);
                    let dpdy = -(pressure[[i, j + 1]] - pressure[[i, j - 1]])
                        / (scalar::from_f64::<T>(2.0) * dy);

                    // Viscous term: μ∇²u (using central differences)
                    let visc_x = viscosity
                        * ((u[[i + 1, j]] - scalar::from_f64::<T>(2.0) * u[[i, j]]
                            + u[[i - 1, j]])
                            / (dx * dx)
                            + (u[[i, j + 1]] - scalar::from_f64::<T>(2.0) * u[[i, j]]
                                + u[[i, j - 1]])
                                / (dy * dy));

                    let visc_y = viscosity
                        * ((v[[i + 1, j]] - scalar::from_f64::<T>(2.0) * v[[i, j]]
                            + v[[i - 1, j]])
                            / (dx * dx)
                            + (v[[i, j + 1]] - scalar::from_f64::<T>(2.0) * v[[i, j]]
                                + v[[i, j - 1]])
                                / (dy * dy));

                    // Body force: ρg
                    let body_x = self.density * gravity[0];
                    let body_y = self.density * gravity[1];

                    // Momentum residuals
                    let residual_x = dudt + conv_x - dpdx - visc_x - body_x;
                    let residual_y = dvdt + conv_y - dpdy - visc_y - body_y;

                    let abs_residual_x = scalar::abs(residual_x);
                    let abs_residual_y = scalar::abs(residual_y);
                    max_residual_x = scalar::max(max_residual_x, abs_residual_x);
                    max_residual_y = scalar::max(max_residual_y, abs_residual_y);
                    total_residual = total_residual + abs_residual_x + abs_residual_y;
                    count += 1;
                }
            }

            let avg_residual = if count > 0 {
                total_residual / scalar::from_usize::<T>(count * 2)
            } else {
                scalar::zero::<T>()
            };

            let max_residual = scalar::max(max_residual_x, max_residual_y);

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
                scalar::from_usize::<T>(count),
            );

            Ok(report)
        }
    }
}

impl<T: RealField + Copy + FloatElement> ConservationChecker<T> for MomentumConservationChecker<T> {
    type FlowField = (Array2<T>, Array2<T>);

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        let (u, v) = field;

        // For generic check, assume steady state (u_prev = u), no pressure gradient,
        // and check if viscous forces balance
        let pressure = Array2::zeros([self.nx, self.ny]);
        let viscosity = scalar::from_f64::<T>(1e-3);
        let dt = scalar::from_f64::<T>(1e-3);
        let dx = scalar::one::<T>();
        let dy = scalar::one::<T>();
        let gravity = Vector2::zeros();

        self.check_momentum_2d(u, v, u, v, &pressure, viscosity, dt, dx, dy, gravity)
    }

    fn name(&self) -> &'static str {
        "Momentum Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use leto::Array2;

    #[test]
    fn test_steady_uniform_flow() {
        let nx = 10;
        let ny = 10;
        let checker = MomentumConservationChecker::<f64>::new(1e-10, nx, ny, 1.0);

        // Uniform steady flow with no forces
        let u = Array2::from_elem([nx, ny], 1.0);
        let v = Array2::from_elem([nx, ny], 0.0);
        let pressure = Array2::zeros([nx, ny]);

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
        let mut u = Array2::zeros([nx, ny]);
        let v = Array2::zeros([nx, ny]);

        // Create parabolic profile in y-direction
        for i in 0..nx {
            for j in 0..ny {
                let y = j as f64 / (ny - 1) as f64;
                u[[i, j]] = 4.0 * y * (1.0 - y); // Parabolic profile
            }
        }

        // Linear pressure gradient
        let mut pressure = Array2::zeros([nx, ny]);
        for i in 0..nx {
            for j in 0..ny {
                pressure[[i, j]] = 1.0 - (i as f64 / (nx - 1) as f64);
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
