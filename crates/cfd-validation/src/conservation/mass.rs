//! Mass conservation checker for CFD simulations
//!
//! Validates that the continuity equation ∇·u = 0 is satisfied

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use crate::scalar;
use cfd_core::error::Result;
use eunomia::FloatElement;
use eunomia::RealField;
use leto::{Array1, Array2};

/// Mass conservation checker
pub struct MassConservationChecker<T: RealField + Copy> {
    tolerance: T,
    nx: usize,
    ny: usize,
}

impl<T: RealField + Copy + FloatElement> MassConservationChecker<T> {
    /// Create new mass conservation checker
    pub fn new(tolerance: T, nx: usize, ny: usize) -> Self {
        Self { tolerance, nx, ny }
    }

    /// Check mass conservation for a 2D velocity field
    ///
    /// Computes divergence ∇·u = ∂u/∂x + ∂v/∂y using central differences
    pub fn check_divergence_2d(
        &self,
        u: &Array2<T>,
        v: &Array2<T>,
        dx: T,
        dy: T,
    ) -> Result<ConservationReport<T>> {
        assert_eq!(u.shape()[0], self.nx);
        assert_eq!(u.shape()[1], self.ny);
        assert_eq!(v.shape()[0], self.nx);
        assert_eq!(v.shape()[1], self.ny);

        let mut max_divergence = scalar::zero::<T>();
        let mut total_divergence = scalar::zero::<T>();
        let mut count = 0;

        // Compute divergence at interior points using central differences
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                // ∂u/∂x using central difference
                let dudx = (u[[i + 1, j]] - u[[i - 1, j]]) / (scalar::from_f64::<T>(2.0) * dx);

                // ∂v/∂y using central difference
                let dvdy = (v[[i, j + 1]] - v[[i, j - 1]]) / (scalar::from_f64::<T>(2.0) * dy);

                // Divergence = ∂u/∂x + ∂v/∂y
                let div = dudx + dvdy;

                let abs_div = scalar::abs(div);
                max_divergence = scalar::max(max_divergence, abs_div);
                total_divergence += abs_div;
                count += 1;
            }
        }

        let avg_divergence = if count > 0 {
            total_divergence / scalar::from_usize::<T>(count)
        } else {
            scalar::zero::<T>()
        };

        let mut report = ConservationReport::new(
            "Mass Conservation (2D)".to_string(),
            max_divergence,
            self.tolerance,
        );

        report.add_detail("max_divergence".to_string(), max_divergence);
        report.add_detail("avg_divergence".to_string(), avg_divergence);
        report.add_detail(
            "grid_points_checked".to_string(),
            scalar::from_usize::<T>(count),
        );

        Ok(report)
    }

    /// Check mass conservation for a 1D velocity field (for pipes/channels)
    pub fn check_divergence_1d(
        &self,
        u: &Array1<T>,
        area: &Array1<T>,
        dx: T,
    ) -> Result<ConservationReport<T>> {
        let n = u.shape()[0];
        assert_eq!(area.shape()[0], n);

        let mut max_divergence = scalar::zero::<T>();
        let mut total_divergence = scalar::zero::<T>();
        let mut count = 0;

        // For 1D flow in variable area: ∂(ρAu)/∂x = 0
        // For incompressible flow: ∂(Au)/∂x = 0
        for i in 1..n - 1 {
            // Mass flux gradient using central difference
            let flux_plus = area[i + 1] * u[i + 1];
            let flux_minus = area[i - 1] * u[i - 1];
            let div = (flux_plus - flux_minus) / (scalar::from_f64::<T>(2.0) * dx);

            let abs_div = scalar::abs(div);
            max_divergence = scalar::max(max_divergence, abs_div);
            total_divergence += abs_div;
            count += 1;
        }

        let avg_divergence = if count > 0 {
            total_divergence / scalar::from_usize::<T>(count)
        } else {
            scalar::zero::<T>()
        };

        let mut report = ConservationReport::new(
            "Mass Conservation (1D)".to_string(),
            max_divergence,
            self.tolerance,
        );

        report.add_detail("max_divergence".to_string(), max_divergence);
        report.add_detail("avg_divergence".to_string(), avg_divergence);
        report.add_detail(
            "grid_points_checked".to_string(),
            scalar::from_usize::<T>(count),
        );

        Ok(report)
    }
}

impl<T: RealField + Copy + FloatElement> ConservationChecker<T> for MassConservationChecker<T> {
    type FlowField = (Array2<T>, Array2<T>);

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        let (u, v) = field;
        // Assume unit spacing for generic check
        let dx = scalar::one::<T>();
        let dy = scalar::one::<T>();

        self.check_divergence_2d(u, v, dx, dy)
    }

    fn name(&self) -> &'static str {
        "Mass Conservation"
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
    fn test_zero_divergence_uniform_flow() {
        let nx = 10;
        let ny = 10;
        let checker = MassConservationChecker::<f64>::new(1e-10, nx, ny);

        // Uniform flow has zero divergence
        let u = Array2::from_elem([nx, ny], 1.0);
        let v = Array2::from_elem([nx, ny], 0.0);

        let report = checker.check_divergence_2d(&u, &v, 0.1, 0.1).unwrap();
        assert!(report.error < 1e-10);
        assert!(report.is_conserved);
    }

    #[test]
    fn test_nonzero_divergence_source() {
        let nx = 10;
        let ny = 10;
        let checker = MassConservationChecker::<f64>::new(1e-6, nx, ny);

        // Create a velocity field with a source (non-zero divergence)
        let mut u = Array2::zeros([nx, ny]);
        let mut v = Array2::zeros([nx, ny]);

        // Radial outflow from center: u = x/r, v = y/r
        // Expected divergence: ∇·u = ∂u/∂x + ∂v/∂y = ∂(x/r)/∂x + ∂(y/r)/∂y
        // For r = sqrt(x² + y²): ∇·u = 1/r - x²/r³ + 1/r - y²/r³ = 2/r - (x² + y²)/r³ = 2/r - r²/r³ = 1/r
        let cx = nx / 2;
        let cy = ny / 2;
        for i in 0..nx {
            for j in 0..ny {
                let dx = i as f64 - cx as f64;
                let dy = j as f64 - cy as f64;
                let r = (dx * dx + dy * dy).sqrt();
                if r > 0.0 {
                    u[[i, j]] = dx / r;
                    v[[i, j]] = dy / r;
                }
            }
        }

        let report = checker.check_divergence_2d(&u, &v, 0.1, 0.1).unwrap();

        // For radial outflow u = x/r, v = y/r, the analytical divergence is 1/r
        // At the center (r ≈ 0), we expect high divergence
        // At r = 1 grid spacing, expected divergence ≈ 1/1 = 1.0
        // The numerical approximation should be close to this theoretical value

        // Verify the error is of the expected magnitude (not just > 0)
        assert!(
            report.error > 0.1,
            "Divergence should be significant for radial outflow: got {}",
            report.error
        );
        assert!(
            report.error < 50.0,
            "Divergence should be bounded for well-behaved field: got {}",
            report.error
        );

        // Should correctly detect this as non-conserved
        assert!(
            !report.is_conserved,
            "Radial outflow should not be mass-conserved"
        );
    }

    #[test]
    fn test_1d_continuity() {
        let n = 10;
        let checker = MassConservationChecker::<f64>::new(1e-10, n, n);

        // Constant area pipe with uniform flow
        let u = Array1::from_elem([n], 2.0);
        let area = Array1::from_elem([n], 1.0);

        let report = checker.check_divergence_1d(&u, &area, 0.1).unwrap();
        assert!(report.error < 1e-10);
        assert!(report.is_conserved);
    }
}
