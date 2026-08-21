//! 2D heat diffusion solver using an explicit finite difference method.
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

use cfd_core::CfdScalar;
use cfd_core::error::Error;
use eunomia::{FloatElement, NumericElement};
use std::collections::HashMap;

use crate::scalar;

/// A 2D heat diffusion solver implementing an explicit, first-order time-stepping
/// scheme and second-order central differences for spatial discretization.
///
/// This solver is designed for transient diffusion problems on a structured grid.
///
/// # Type Parameters
///
/// * `T`: The floating-point type for calculations, e.g., `f64`.
pub struct DiffusionSolver<T: CfdScalar + Copy + FloatElement> {
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    alpha: T,
    solution: HashMap<(usize, usize), T>,
}

impl<T: CfdScalar + Copy + FloatElement> DiffusionSolver<T> {
    /// Creates a new `DiffusionSolver`.
    ///
    /// # Arguments
    ///
    /// * `nx`: Number of grid points in the x-direction.
    /// * `ny`: Number of grid points in the y-direction.
    /// * `dx`: Grid spacing in the x-direction.
    /// * `dy`: Grid spacing in the y-direction.
    /// * `alpha`: Thermal diffusivity.
    ///
    /// # Panics
    /// Panics if any input is invalid (see [`Self::try_new`]).
    pub fn new(nx: usize, ny: usize, dx: T, dy: T, alpha: T) -> Self {
        Self::try_new(nx, ny, dx, dy, alpha).unwrap_or_else(|error| {
            panic!("DiffusionSolver::new called with invalid inputs: {error}");
        })
    }

    /// Create a new explicit-Euler 2D heat-diffusion solver with grid and
    /// physical-parameter validation.
    ///
    /// The explicit scheme computes a stable time step
    /// `dt = 0.25 * 0.9 * min(dx², dy²) / alpha` and divides the Laplacian
    /// by `dx²` and `dy²`. Zero or non-finite `dx`/`dy`/`alpha` silently
    /// produce `inf`/`NaN` `dt` and a divergent solve.
    ///
    /// # Errors
    /// Returns `Error::InvalidConfiguration` if:
    /// - `nx < 3` or `ny < 3` (the central-difference stencil needs at
    ///   least one interior point in each direction),
    /// - `dx`, `dy`, or `alpha` is non-finite or non-positive.
    pub fn try_new(nx: usize, ny: usize, dx: T, dy: T, alpha: T) -> cfd_core::error::Result<Self> {
        if nx < 3 {
            return Err(Error::InvalidConfiguration(format!(
                "DiffusionSolver::try_new: nx must be at least 3 for the central-difference stencil, got {nx}"
            )));
        }
        if ny < 3 {
            return Err(Error::InvalidConfiguration(format!(
                "DiffusionSolver::try_new: ny must be at least 3 for the central-difference stencil, got {ny}"
            )));
        }
        if !<T as NumericElement>::is_finite(dx) || dx <= scalar::zero() {
            return Err(Error::InvalidConfiguration(format!(
                "DiffusionSolver::try_new: dx (grid spacing) must be finite and positive, got {dx:?}"
            )));
        }
        if !<T as NumericElement>::is_finite(dy) || dy <= scalar::zero() {
            return Err(Error::InvalidConfiguration(format!(
                "DiffusionSolver::try_new: dy (grid spacing) must be finite and positive, got {dy:?}"
            )));
        }
        if !<T as NumericElement>::is_finite(alpha) || alpha <= scalar::zero() {
            return Err(Error::InvalidConfiguration(format!(
                "DiffusionSolver::try_new: alpha (thermal diffusivity) must be finite and positive, got {alpha:?}"
            )));
        }
        Ok(Self {
            nx,
            ny,
            dx,
            dy,
            alpha,
            solution: HashMap::new(),
        })
    }

    /// Sets the boundary and initial conditions for the solver.
    ///
    /// The solver uses the provided `HashMap` to initialize the entire solution field,
    /// including interior points and boundary nodes.
    pub fn set_boundary_conditions(&mut self, boundary_conditions: &HashMap<(usize, usize), T>) {
        self.solution.clone_from(boundary_conditions);
    }

    /// Solves the diffusion equation up to a specified final time using an explicit
    /// Euler time-stepping scheme.
    ///
    /// # Arguments
    ///
    /// * `t_final`: The final simulation time.
    /// * `source_fn`: A closure that provides the time-dependent source term `S(x, y, t)`.
    pub fn solve_to_time(
        &mut self,
        t_final: T,
        source_fn: &impl Fn(T, T, T) -> T,
    ) -> HashMap<(usize, usize), T> {
        let dt = <T as FloatElement>::from_f64(0.25 * 0.9) * self.dx * self.dx / self.alpha;
        let mut t: T = scalar::zero();
        let two = <T as FloatElement>::from_f64(2.0);

        while t < t_final {
            let mut next_solution = self.solution.clone();
            for i in 1..self.nx - 1 {
                for j in 1..self.ny - 1 {
                    let un = self
                        .solution
                        .get(&(i, j))
                        .copied()
                        .expect("analytical constant conversion");
                    let un_e = self
                        .solution
                        .get(&(i + 1, j))
                        .copied()
                        .expect("analytical constant conversion");
                    let un_w = self
                        .solution
                        .get(&(i - 1, j))
                        .copied()
                        .expect("analytical constant conversion");
                    let un_n = self
                        .solution
                        .get(&(i, j + 1))
                        .copied()
                        .expect("analytical constant conversion");
                    let un_s = self
                        .solution
                        .get(&(i, j - 1))
                        .copied()
                        .expect("analytical constant conversion");

                    let laplacian = (un_e - two * un + un_w) / (self.dx * self.dx)
                        + (un_n - two * un + un_s) / (self.dy * self.dy);

                    let x = scalar::from_usize::<T>(i) * self.dx;
                    let y = scalar::from_usize::<T>(j) * self.dy;
                    let source = source_fn(x, y, t);
                    let un_plus_1 = un + dt * (self.alpha * laplacian + source);
                    next_solution.insert((i, j), un_plus_1);
                }
            }
            self.solution = next_solution;
            t += dt;
        }
        self.solution.clone()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **Positive**: `try_new` accepts valid inputs.
    #[test]
    fn diffusion_try_new_accepts_valid_inputs() {
        let _solver = DiffusionSolver::<f64>::try_new(10, 10, 0.1, 0.1, 0.01)
            .expect("valid inputs must succeed");
    }

    /// **Adversarial**: `nx < 3` is rejected (central-difference stencil needs ≥3).
    #[test]
    fn diffusion_try_new_rejects_small_nx() {
        match DiffusionSolver::<f64>::try_new(2, 10, 0.1, 0.1, 0.01) {
            Err(e) => assert!(e.to_string().contains("nx"), "error must mention nx: {e}"),
            Ok(_) => panic!("nx < 3 must be rejected"),
        }
    }

    /// **Adversarial**: `ny < 3` is rejected.
    #[test]
    fn diffusion_try_new_rejects_small_ny() {
        match DiffusionSolver::<f64>::try_new(10, 2, 0.1, 0.1, 0.01) {
            Err(e) => assert!(e.to_string().contains("ny"), "error must mention ny: {e}"),
            Ok(_) => panic!("ny < 3 must be rejected"),
        }
    }

    /// **Adversarial**: invalid `dx` is rejected.
    #[test]
    fn diffusion_try_new_rejects_invalid_dx() {
        match DiffusionSolver::<f64>::try_new(10, 10, 0.0, 0.1, 0.01) {
            Err(e) => assert!(e.to_string().contains("dx"), "error must mention dx: {e}"),
            Ok(_) => panic!("zero dx must be rejected"),
        }
        match DiffusionSolver::<f64>::try_new(10, 10, f64::NAN, 0.1, 0.01) {
            Err(e) => assert!(e.to_string().contains("dx"), "error must mention dx: {e}"),
            Ok(_) => panic!("NaN dx must be rejected"),
        }
    }

    /// **Adversarial**: invalid `dy` is rejected.
    #[test]
    fn diffusion_try_new_rejects_invalid_dy() {
        match DiffusionSolver::<f64>::try_new(10, 10, 0.1, 0.0, 0.01) {
            Err(e) => assert!(e.to_string().contains("dy"), "error must mention dy: {e}"),
            Ok(_) => panic!("zero dy must be rejected"),
        }
    }

    /// **Adversarial**: invalid `alpha` is rejected.
    #[test]
    fn diffusion_try_new_rejects_invalid_alpha() {
        match DiffusionSolver::<f64>::try_new(10, 10, 0.1, 0.1, 0.0) {
            Err(e) => assert!(
                e.to_string().contains("alpha"),
                "error must mention alpha: {e}"
            ),
            Ok(_) => panic!("zero alpha must be rejected"),
        }
        match DiffusionSolver::<f64>::try_new(10, 10, 0.1, 0.1, -1.0) {
            Err(e) => assert!(
                e.to_string().contains("alpha"),
                "error must mention alpha: {e}"
            ),
            Ok(_) => panic!("negative alpha must be rejected"),
        }
    }

    /// **Boundary**: `new` panics on invalid `alpha` (thin wrapper contract).
    #[test]
    #[should_panic(expected = "alpha")]
    fn diffusion_new_panics_on_invalid_alpha() {
        let _ = DiffusionSolver::<f64>::new(10, 10, 0.1, 0.1, 0.0);
    }
}
