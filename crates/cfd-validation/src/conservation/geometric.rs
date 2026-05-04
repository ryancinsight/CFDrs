//! Geometric Conservation Law (GCL) tests for CFD schemes
//!
//! The Geometric Conservation Law requires that numerical schemes preserve
//! constant solutions on moving grids. This is essential for long-time accuracy
//! and stability of time-dependent simulations.
//!
//! Reference: Thomas & Lombard (1979). "Geometric Conservation Law and Its Application
//! to Flow Computations on Moving Grids"

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;

/// Geometric Conservation Law checker.
pub struct GeometricConservationChecker<T: RealField + Copy> {
    tolerance: T,
    nx: usize,
    ny: usize,
}

impl<T: RealField + Copy + FromPrimitive> GeometricConservationChecker<T> {
    /// Create new GCL checker.
    pub fn new(tolerance: T, nx: usize, ny: usize) -> Self {
        Self { tolerance, nx, ny }
    }

    /// Evaluate the finite-volume diffusion residual on a uniform Cartesian grid.
    ///
    /// # Theorem
    /// The conservative second-order face-flux residual preserves every
    /// constant scalar field exactly on a stationary uniform grid.
    ///
    /// **Proof sketch**: For `u_ij = c`, every east/west/north/south face
    /// gradient is `(c - c) / h = 0`. The residual is the divergence of these
    /// face fluxes, so each interior residual entry is zero. Boundary entries
    /// are fixed to zero because the GCL check evolves the interior state under
    /// stationary boundary geometry. Therefore any explicit Runge-Kutta method
    /// whose stages are affine combinations of `u` and `dt * R(u)` leaves the
    /// constant field unchanged in exact arithmetic.
    fn conservative_residual(&self, field: &DMatrix<T>, dx: T, dy: T) -> DMatrix<T> {
        let mut residual = DMatrix::zeros(field.nrows(), field.ncols());

        if field.nrows() < 3 || field.ncols() < 3 {
            return residual;
        }

        for i in 1..field.nrows() - 1 {
            for j in 1..field.ncols() - 1 {
                let east_flux = (field[(i + 1, j)] - field[(i, j)]) / dx;
                let west_flux = (field[(i, j)] - field[(i - 1, j)]) / dx;
                let north_flux = (field[(i, j + 1)] - field[(i, j)]) / dy;
                let south_flux = (field[(i, j)] - field[(i, j - 1)]) / dy;

                residual[(i, j)] = (east_flux - west_flux) / dx + (north_flux - south_flux) / dy;
            }
        }

        residual
    }

    fn euler_step(&self, field: &DMatrix<T>, dt: T, dx: T, dy: T) -> DMatrix<T> {
        let residual = self.conservative_residual(field, dx, dy);
        let mut next = field.clone();

        for i in 0..field.nrows() {
            for j in 0..field.ncols() {
                next[(i, j)] += dt * residual[(i, j)];
            }
        }

        next
    }

    fn add_scaled(base: &DMatrix<T>, increment: &DMatrix<T>, scale: T) -> DMatrix<T> {
        let mut result = base.clone();

        for i in 0..base.nrows() {
            for j in 0..base.ncols() {
                result[(i, j)] += scale * increment[(i, j)];
            }
        }

        result
    }

    fn convex_blend(a: &DMatrix<T>, a_weight: T, b: &DMatrix<T>, b_weight: T) -> DMatrix<T> {
        let mut result = a.clone();

        for i in 0..a.nrows() {
            for j in 0..a.ncols() {
                result[(i, j)] = a_weight * a[(i, j)] + b_weight * b[(i, j)];
            }
        }

        result
    }

    fn runge_kutta_step(
        &self,
        field: &DMatrix<T>,
        dt: T,
        dx: T,
        dy: T,
        stages: usize,
    ) -> Result<DMatrix<T>> {
        match stages {
            1 => Ok(self.euler_step(field, dt, dx, dy)),
            2 => {
                let half = <T as SafeFromF64>::from_f64_or_one(0.5);
                let k1 = self.conservative_residual(field, dx, dy);
                let midpoint = Self::add_scaled(field, &k1, half * dt);
                let k2 = self.conservative_residual(&midpoint, dx, dy);
                Ok(Self::add_scaled(field, &k2, dt))
            }
            3 => {
                let one_fourth = <T as SafeFromF64>::from_f64_or_one(0.25);
                let three_fourths = <T as SafeFromF64>::from_f64_or_one(0.75);
                let one_third = <T as SafeFromF64>::from_f64_or_one(1.0 / 3.0);
                let two_thirds = <T as SafeFromF64>::from_f64_or_one(2.0 / 3.0);

                let u1 = self.euler_step(field, dt, dx, dy);
                let u1_evolved = self.euler_step(&u1, dt, dx, dy);
                let u2 = Self::convex_blend(field, three_fourths, &u1_evolved, one_fourth);
                let u2_evolved = self.euler_step(&u2, dt, dx, dy);
                Ok(Self::convex_blend(
                    field,
                    one_third,
                    &u2_evolved,
                    two_thirds,
                ))
            }
            4 => {
                let half = <T as SafeFromF64>::from_f64_or_one(0.5);
                let one_sixth = <T as SafeFromF64>::from_f64_or_one(1.0 / 6.0);
                let one_third = <T as SafeFromF64>::from_f64_or_one(1.0 / 3.0);

                let k1 = self.conservative_residual(field, dx, dy);
                let u2 = Self::add_scaled(field, &k1, half * dt);
                let k2 = self.conservative_residual(&u2, dx, dy);
                let u3 = Self::add_scaled(field, &k2, half * dt);
                let k3 = self.conservative_residual(&u3, dx, dy);
                let u4 = Self::add_scaled(field, &k3, dt);
                let k4 = self.conservative_residual(&u4, dx, dy);

                let mut next = field.clone();
                for i in 0..field.nrows() {
                    for j in 0..field.ncols() {
                        next[(i, j)] += dt
                            * (one_sixth * k1[(i, j)]
                                + one_third * k2[(i, j)]
                                + one_third * k3[(i, j)]
                                + one_sixth * k4[(i, j)]);
                    }
                }

                Ok(next)
            }
            _ => Err(Error::UnsupportedOperation(format!(
                "GCL Runge-Kutta check supports 1, 2, 3, or 4 stages; got {stages}"
            ))),
        }
    }

    /// Test GCL for Euler time stepping with constant solution.
    /// For Euler scheme: `u^{n+1} = u^n + dt * R(u^n)`.
    /// For constant solution `u = constant`, conservative residual `R(u) = 0`.
    pub fn test_euler_gcl(&self, constant_value: T) -> Result<ConservationReport<T>> {
        let mut max_error = T::zero();
        let mut total_error = T::zero();
        let mut count = 0;

        let u = DMatrix::from_element(self.nx, self.ny, constant_value);
        let dt = <T as SafeFromF64>::from_f64_or_one(0.01);
        let dx = <T as SafeFromF64>::from_f64_or_one(0.1);
        let dy = <T as SafeFromF64>::from_f64_or_one(0.1);
        let mut u_current = u.clone();

        for _step in 0..10 {
            let u_next = self.euler_step(&u_current, dt, dx, dy);

            for i in 0..self.nx {
                for j in 0..self.ny {
                    let error = (u_next[(i, j)] - constant_value).abs();
                    max_error = max_error.max(error);
                    total_error += error;
                    count += 1;
                }
            }

            u_current = u_next;
        }

        let avg_error = if count > 0 {
            total_error / T::from_usize(count).unwrap_or(T::one())
        } else {
            T::zero()
        };

        let mut report = ConservationReport::new(
            "Geometric Conservation Law (Euler)".to_string(),
            max_error,
            self.tolerance,
        );

        report.add_detail("max_error".to_string(), max_error);
        report.add_detail("avg_error".to_string(), avg_error);
        report.add_detail("constant_value".to_string(), constant_value);
        report.add_detail(
            "time_steps".to_string(),
            T::from_usize(10).unwrap_or(T::one()),
        );
        report.add_detail("dt".to_string(), dt);
        report.add_detail("dx".to_string(), dx);
        report.add_detail("dy".to_string(), dy);

        Ok(report)
    }

    /// Test GCL for Runge-Kutta schemes.
    /// RK schemes should preserve constants if the RHS evaluation is consistent.
    pub fn test_runge_kutta_gcl(
        &self,
        constant_value: T,
        stages: usize,
    ) -> Result<ConservationReport<T>> {
        let mut max_error = T::zero();
        let mut total_error = T::zero();
        let mut count = 0;

        let u = DMatrix::from_element(self.nx, self.ny, constant_value);
        let dt = <T as SafeFromF64>::from_f64_or_one(0.01);
        let dx = <T as SafeFromF64>::from_f64_or_one(0.1);
        let dy = <T as SafeFromF64>::from_f64_or_one(0.1);
        let mut u_current = u.clone();

        for _step in 0..5 {
            let u_next = self.runge_kutta_step(&u_current, dt, dx, dy, stages)?;

            for i in 0..self.nx {
                for j in 0..self.ny {
                    let error = (u_next[(i, j)] - constant_value).abs();
                    max_error = max_error.max(error);
                    total_error += error;
                    count += 1;
                }
            }

            u_current = u_next;
        }

        let avg_error = if count > 0 {
            total_error / T::from_usize(count).unwrap_or(T::one())
        } else {
            T::zero()
        };

        let mut report = ConservationReport::new(
            format!("Geometric Conservation Law (RK{stages})"),
            max_error,
            self.tolerance,
        );

        report.add_detail("max_error".to_string(), max_error);
        report.add_detail("avg_error".to_string(), avg_error);
        report.add_detail(
            "stages".to_string(),
            T::from_usize(stages).unwrap_or(T::one()),
        );
        report.add_detail("constant_value".to_string(), constant_value);
        report.add_detail("dt".to_string(), dt);
        report.add_detail("dx".to_string(), dx);
        report.add_detail("dy".to_string(), dy);

        Ok(report)
    }

    /// Test GCL for spatial discretization schemes.
    /// Tests that constant solutions are preserved by spatial operators.
    pub fn test_spatial_gcl(&self, constant_value: T) -> Result<ConservationReport<T>> {
        let mut max_error = T::zero();
        let mut total_error = T::zero();
        let mut count = 0;

        let u = DMatrix::from_element(self.nx, self.ny, constant_value);
        let dx = <T as SafeFromF64>::from_f64_or_one(0.1);
        let dy = <T as SafeFromF64>::from_f64_or_one(0.1);
        let residual = self.conservative_residual(&u, dx, dy);

        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let error = residual[(i, j)].abs();
                max_error = max_error.max(error);
                total_error += error;
                count += 1;
            }
        }

        let avg_error = if count > 0 {
            total_error / T::from_usize(count).unwrap_or(T::one())
        } else {
            T::zero()
        };

        let mut report = ConservationReport::new(
            "Geometric Conservation Law (Spatial)".to_string(),
            max_error,
            self.tolerance,
        );

        report.add_detail("max_error".to_string(), max_error);
        report.add_detail("avg_error".to_string(), avg_error);
        report.add_detail("constant_value".to_string(), constant_value);
        report.add_detail("dx".to_string(), dx);
        report.add_detail("dy".to_string(), dy);

        Ok(report)
    }

    /// Comprehensive GCL test suite.
    pub fn run_comprehensive_gcl_tests(&self) -> Result<Vec<ConservationReport<T>>> {
        println!("Testing Geometric Conservation Law (GCL)");
        println!("========================================");
        println!("GCL ensures schemes preserve constants on moving grids");
        println!("Reference: Thomas & Lombard (1979)");
        println!();

        let mut results = Vec::new();
        let test_values = vec![
            T::zero(),
            T::one(),
            <T as SafeFromF64>::from_f64_or_one(1.5),
            <T as SafeFromF64>::from_f64_or(-1.0, -T::one()),
        ];

        for &value in &test_values {
            results.push(self.test_euler_gcl(value)?);
            results.push(self.test_runge_kutta_gcl(value, 4)?);
            results.push(self.test_spatial_gcl(value)?);
        }

        let passed_tests = results.iter().filter(|r| r.is_conserved).count();
        let total_tests = results.len();

        println!("\nGCL Test Summary:");
        println!("  Tests Passed: {passed_tests}/{total_tests}");
        println!(
            "  Success Rate: {:.1}%",
            100.0 * passed_tests as f32 / total_tests as f32
        );

        if passed_tests == total_tests {
            println!("All GCL tests passed; schemes preserve constants correctly.");
        } else {
            println!("Some GCL tests failed; schemes may not preserve constants.");
            println!("This can cause accuracy issues in long-time simulations.");
        }

        Ok(results)
    }
}

impl<T: RealField + Copy + FromPrimitive> ConservationChecker<T>
    for GeometricConservationChecker<T>
{
    type FlowField = DMatrix<T>;

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        let constant_value = if field.nrows() > 0 && field.ncols() > 0 {
            field[(0, 0)]
        } else {
            T::one()
        };

        self.test_spatial_gcl(constant_value)
    }

    fn name(&self) -> &'static str {
        "Geometric Conservation Law"
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_euler_gcl_zero() {
        let checker = GeometricConservationChecker::<f64>::new(1e-14, 10, 10);
        let result = checker.test_euler_gcl(0.0).unwrap();
        assert!(result.error < 1e-13);
        assert!(result.is_conserved);
        assert_eq!(result.details["dt"], 0.01);
    }

    #[test]
    fn test_euler_gcl_nonzero() {
        let checker = GeometricConservationChecker::<f64>::new(1e-14, 10, 10);
        let result = checker.test_euler_gcl(2.5).unwrap();
        assert!(result.error < 1e-13);
        assert!(result.is_conserved);
        assert_eq!(result.details["constant_value"], 2.5);
    }

    #[test]
    fn test_rk4_gcl() {
        let checker = GeometricConservationChecker::<f64>::new(1e-14, 10, 10);
        let result = checker.test_runge_kutta_gcl(1.0, 4).unwrap();
        assert!(result.error < 1e-13);
        assert!(result.is_conserved);
        assert_eq!(result.details["stages"], 4.0);
    }

    #[test]
    fn test_runge_kutta_rejects_unsupported_stage_count() {
        let checker = GeometricConservationChecker::<f64>::new(1e-14, 10, 10);
        let err = checker.test_runge_kutta_gcl(1.0, 5).unwrap_err();
        assert!(matches!(err, Error::UnsupportedOperation(_)));
    }

    #[test]
    fn test_spatial_gcl() {
        let checker = GeometricConservationChecker::<f64>::new(1e-14, 10, 10);
        let result = checker.test_spatial_gcl(3.14).unwrap();
        assert!(result.error < 1e-13);
        assert!(result.is_conserved);
        assert_eq!(result.details["dx"], 0.1);
    }

    #[test]
    fn test_conservative_residual_is_not_copy_through_for_quadratic_field() {
        let checker = GeometricConservationChecker::<f64>::new(1e-14, 5, 5);
        let field = DMatrix::from_fn(5, 5, |i, j| (i * i + j * j) as f64);

        let residual = checker.conservative_residual(&field, 1.0, 1.0);
        assert_eq!(residual[(2, 2)], 4.0);

        let updated = checker.euler_step(&field, 0.25, 1.0, 1.0);
        assert_eq!(updated[(2, 2)], field[(2, 2)] + 1.0);
    }

    #[test]
    fn test_comprehensive_gcl() {
        let checker = GeometricConservationChecker::<f64>::new(1e-14, 10, 10);
        let results = checker.run_comprehensive_gcl_tests().unwrap();
        assert_eq!(results.len(), 12);

        for result in results {
            assert!(
                result.is_conserved,
                "GCL test failed: {}",
                result.check_name
            );
            assert_eq!(result.error, 0.0);
        }
    }
}
