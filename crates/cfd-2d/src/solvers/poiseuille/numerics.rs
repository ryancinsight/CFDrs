//! Numerical routines for the Poiseuille flow solver
//!
//! Thomas algorithm (tridiagonal solver), internal solver steps:
//! velocity solve, shear rate calculation, viscosity update, and
//! convergence check.

use super::{BloodModel, PoiseuilleFlow2D};
use crate::error::Error;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

impl<T: RealField + FromPrimitive + Float + Copy> PoiseuilleFlow2D<T> {
    /// Solve linear system for velocity with current viscosity
    ///
    /// Solves tridiagonal system using Thomas algorithm (direct solver)
    ///
    /// Matrix equation: A*u = b where
    /// - A is tridiagonal
    /// - u is velocity vector
    /// - b is source term (pressure gradient)
    pub(super) fn solve_velocity_linear(&mut self) -> Result<(), Error> {
        let ny = self.config.ny;
        let dy = self.dy;
        let dp_dx = self.config.pressure_gradient;

        // Allocate tridiagonal matrix coefficients
        let mut a = vec![T::zero(); ny]; // Sub-diagonal
        let mut b = vec![T::zero(); ny]; // Diagonal
        let mut c = vec![T::zero(); ny]; // Super-diagonal
        let mut d = vec![T::zero(); ny]; // RHS

        // Boundary conditions: u(0) = 0, u(ny-1) = 0
        b[0] = T::one();
        c[0] = T::zero();
        d[0] = T::zero();

        b[ny - 1] = T::one();
        a[ny - 1] = T::zero();
        d[ny - 1] = T::zero();

        // Interior points: j = 1, 2, ..., ny-2
        let dy2 = dy * dy;

        for j in 1..ny - 1 {
            // Viscosity at cell faces (harmonic mean for better stability)
            let mu_jm12 = T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero)
                / (T::one() / self.viscosity[j - 1] + T::one() / self.viscosity[j]);
            let mu_jp12 = T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero)
                / (T::one() / self.viscosity[j] + T::one() / self.viscosity[j + 1]);

            // Coefficients for interior stencil
            a[j] = -mu_jm12 / dy2;
            c[j] = -mu_jp12 / dy2;
            b[j] = mu_jm12 / dy2 + mu_jp12 / dy2;
            d[j] = dp_dx; // Source term from pressure gradient
        }

        // Solve tridiagonal system using Thomas algorithm
        self.velocity = thomas_algorithm(&a, &b, &c, &d)?;

        Ok(())
    }

    /// Calculate shear rate from velocity gradient
    ///
    /// γ̇(y) = |du/dy|
    ///
    /// Uses central differences for interior points:
    /// du/dy ≈ (u_{j+1} - u_{j-1}) / (2Δy)
    pub(super) fn calculate_shear_rate(&mut self) {
        let ny = self.config.ny;
        let dy = self.dy;
        let two = T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

        // Boundaries: use one-sided differences
        self.shear_rate[0] = Float::abs((self.velocity[1] - self.velocity[0]) / dy);
        self.shear_rate[ny - 1] = Float::abs((self.velocity[ny - 1] - self.velocity[ny - 2]) / dy);

        // Interior: central differences
        for j in 1..ny - 1 {
            let du_dy = (self.velocity[j + 1] - self.velocity[j - 1]) / (two * dy);
            self.shear_rate[j] = Float::abs(du_dy);
        }
    }

    /// Update viscosity from shear rate using blood model
    ///
    /// μ(y) = μ(γ̇(y))
    pub(super) fn update_viscosity_from_shear_rate(&mut self) {
        for j in 0..self.config.ny {
            self.viscosity[j] = match &self.blood_model {
                BloodModel::Casson(casson) => casson.apparent_viscosity(self.shear_rate[j]),
                BloodModel::CarreauYasuda(carreau) => {
                    carreau.apparent_viscosity(self.shear_rate[j])
                }
            };
        }
    }

    /// Calculate L2 norm of viscosity change (for convergence check)
    pub(super) fn calculate_viscosity_residual(&self, viscosity_old: &[T]) -> T {
        let mut sum_sq = T::zero();
        let mut sum_sq_old = T::zero();

        for j in 0..self.config.ny {
            let diff = self.viscosity[j] - viscosity_old[j];
            sum_sq += diff * diff;
            sum_sq_old += viscosity_old[j] * viscosity_old[j];
        }

        // Relative L2 norm
        Float::sqrt(
            sum_sq / (sum_sq_old + T::from_f64(1e-20).unwrap_or_else(num_traits::Zero::zero)),
        )
    }
}

/// Thomas algorithm for tridiagonal systems
///
/// Solves Ax = d where A is tridiagonal:
/// - a: sub-diagonal (a[0] unused)
/// - b: diagonal
/// - c: super-diagonal (c[n-1] unused)
/// - d: right-hand side
///
/// Returns solution vector x
pub(super) fn thomas_algorithm<T: RealField + Copy + Float>(
    a: &[T],
    b: &[T],
    c: &[T],
    d: &[T],
) -> Result<Vec<T>, Error> {
    let n = b.len();

    if a.len() != n || c.len() != n || d.len() != n {
        return Err(Error::InvalidInput("Array lengths must match".to_string()));
    }

    let mut c_prime = vec![T::zero(); n];
    let mut d_prime = vec![T::zero(); n];
    let mut x = vec![T::zero(); n];

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for i in 1..n {
        let denom = b[i] - a[i] * c_prime[i - 1];

        if Float::abs(denom) < T::from_f64(1e-14).unwrap_or_else(num_traits::Zero::zero) {
            return Err(Error::NumericalInstability(
                "Near-zero pivot in Thomas algorithm".to_string(),
            ));
        }

        c_prime[i] = c[i] / denom;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    // Back substitution
    x[n - 1] = d_prime[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    Ok(x)
}
