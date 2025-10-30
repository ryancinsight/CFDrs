//! 2D heat diffusion solver using an explicit finite difference method.

use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// A 2D heat diffusion solver implementing an explicit, first-order time-stepping
/// scheme and second-order central differences for spatial discretization.
///
/// This solver is designed for transient diffusion problems on a structured grid.
///
/// # Type Parameters
///
/// * `T`: The floating-point type for calculations, e.g., `f64`.
pub struct DiffusionSolver<T: RealField + Copy + FromPrimitive> {
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    alpha: T,
    solution: HashMap<(usize, usize), T>,
}

impl<T: RealField + Copy + FromPrimitive> DiffusionSolver<T> {
    /// Creates a new `DiffusionSolver`.
    ///
    /// # Arguments
    ///
    /// * `nx`: Number of grid points in the x-direction.
    /// * `ny`: Number of grid points in the y-direction.
    /// * `dx`: Grid spacing in the x-direction.
    /// * `dy`: Grid spacing in the y-direction.
    /// * `alpha`: Thermal diffusivity.
    pub fn new(nx: usize, ny: usize, dx: T, dy: T, alpha: T) -> Self {
        Self {
            nx,
            ny,
            dx,
            dy,
            alpha,
            solution: HashMap::new(),
        }
    }

    /// Sets the boundary and initial conditions for the solver.
    ///
    /// The solver uses the provided `HashMap` to initialize the entire solution field,
    /// including interior points and boundary nodes.
    pub fn set_boundary_conditions(&mut self, boundary_conditions: &HashMap<(usize, usize), T>) {
        self.solution = boundary_conditions.clone();
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
        let dt = T::from_f64(0.25 * 0.9).unwrap() * self.dx * self.dx / self.alpha;
        let mut t = T::zero();

        while t < t_final {
            let mut next_solution = self.solution.clone();
            for i in 1..self.nx - 1 {
                for j in 1..self.ny - 1 {
                    let un = self.solution.get(&(i, j)).copied().unwrap_or_else(T::zero);
                    let un_e = self.solution.get(&(i + 1, j)).copied().unwrap_or_else(T::zero);
                    let un_w = self.solution.get(&(i - 1, j)).copied().unwrap_or_else(T::zero);
                    let un_n = self.solution.get(&(i, j + 1)).copied().unwrap_or_else(T::zero);
                    let un_s = self.solution.get(&(i, j - 1)).copied().unwrap_or_else(T::zero);

                    let laplacian = (un_e - T::from_f64(2.0).unwrap() * un + un_w) / (self.dx * self.dx)
                        + (un_n - T::from_f64(2.0).unwrap() * un + un_s) / (self.dy * self.dy);

                    let x = T::from_usize(i).unwrap() * self.dx;
                    let y = T::from_usize(j).unwrap() * self.dy;
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
