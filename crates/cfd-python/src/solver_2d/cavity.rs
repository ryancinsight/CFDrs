//! 2D Lid-Driven Cavity solver `PyO3` wrapper — Ghia et al. (1982) benchmark.

use super::pyarray2_from_leto;
use cfd_validation::benchmarks::cavity::LidDrivenCavity;
use leto::Array2;
use numpy::PyArray2;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

// ── Solver ───────────────────────────────────────────────────────────────────

/// 2D Lid-Driven Cavity solver for Ghia et al. (1982) benchmark validation
///
/// Solves incompressible Navier-Stokes equations in a square cavity:
///
/// ```text
/// ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u
/// ∇·u = 0
/// ```
///
/// Boundary conditions:
/// - Top wall: u = `U_lid`, v = 0 (moving lid)
/// - Other walls: u = v = 0 (no-slip)
///
/// Reference:
///   Ghia, U.K.N.G., Ghia, K.N., Shin, C.T. (1982). "High-Re solutions for
///   incompressible flow using the Navier-Stokes equations and a multigrid
///   method". Journal of Computational Physics, 48(3):387-411.
#[pyclass(name = "CavitySolver2D", skip_from_py_object)]
pub struct PyCavitySolver2D {
    #[pyo3(get)]
    reynolds: f64,
    #[pyo3(get)]
    nx: usize,
    #[pyo3(get)]
    ny: usize,
    #[pyo3(get)]
    lid_velocity: f64,
    #[pyo3(get)]
    cavity_size: f64,
}

#[pymethods]
impl PyCavitySolver2D {
    /// Create new lid-driven cavity solver
    #[new]
    #[pyo3(signature = (reynolds=100.0, nx=129, ny=129, lid_velocity=1.0, cavity_size=1.0))]
    fn new(reynolds: f64, nx: usize, ny: usize, lid_velocity: f64, cavity_size: f64) -> Self {
        PyCavitySolver2D {
            reynolds,
            nx,
            ny,
            lid_velocity,
            cavity_size,
        }
    }

    /// Get kinematic viscosity from Reynolds number
    fn viscosity(&self) -> f64 {
        self.lid_velocity * self.cavity_size / self.reynolds
    }

    /// Get Ghia benchmark data for U-velocity along vertical centerline
    fn ghia_u_centerline<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let reference = self.ghia_u_reference();
        if reference.is_empty() {
            return Err(PyValueError::new_err(format!(
                "Ghia Table I has no u-centerline column for Re={}",
                self.reynolds
            )));
        }
        let array = Array2::from_shape_fn([reference.len(), 2], |[i, j]| {
            let &(coordinate, velocity) = reference
                .get(i)
                .expect("invariant: array shape matches canonical Ghia data");
            if j == 0 {
                coordinate
            } else {
                velocity
            }
        });
        pyarray2_from_leto(py, array)
    }

    /// Get Ghia benchmark data for V-velocity along horizontal centerline
    fn ghia_v_centerline<'py>(&self, py: Python<'py>) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let reference = self.ghia_v_reference();
        if reference.is_empty() {
            return Err(PyValueError::new_err(format!(
                "Ghia Table II has no v-centerline column for Re={}",
                self.reynolds
            )));
        }
        let array = Array2::from_shape_fn([reference.len(), 2], |[i, j]| {
            let &(coordinate, velocity) = reference
                .get(i)
                .expect("invariant: array shape matches canonical Ghia data");
            if j == 0 {
                coordinate
            } else {
                velocity
            }
        });
        pyarray2_from_leto(py, array)
    }

    /// Solve lid-driven cavity flow
    fn solve(&self) -> PyResult<PyCavityResult2D> {
        use cfd_2d::solvers::cavity_solver;

        let max_iterations = if self.nx * self.ny > 4000 { 5000 } else { 3000 };
        let result = cavity_solver::solve_lid_driven_cavity(
            self.nx,
            self.ny,
            self.reynolds,
            self.lid_velocity,
            self.cavity_size,
            max_iterations,
            1e-6,
            0.5,
            0.3,
        );

        let l2_error = self.calculate_ghia_error(&result.y_coords, &result.u_centerline)?;

        Ok(PyCavityResult2D {
            l2_error,
            u_centerline: result.u_centerline,
            v_centerline: result.v_centerline,
            y_coords: result.y_coords,
            x_coords: result.x_coords,
            converged: result.converged,
        })
    }

    fn __str__(&self) -> String {
        format!(
            "CavitySolver2D(Re={:.0}, grid={}×{}, U_lid={:.2} m/s)",
            self.reynolds, self.nx, self.ny, self.lid_velocity
        )
    }
}

/// Internal implementation for `PyCavitySolver2D` (non-PyO3 methods)
impl PyCavitySolver2D {
    fn ghia_u_reference(&self) -> Vec<(f64, f64)> {
        LidDrivenCavity::new(1.0, 1.0, self.reynolds).ghia_u_centerline(self.reynolds)
    }

    fn ghia_v_reference(&self) -> Vec<(f64, f64)> {
        LidDrivenCavity::new(1.0, 1.0, self.reynolds).ghia_v_centerline(self.reynolds)
    }

    /// Calculate L2 error between computed and Ghia benchmark data
    fn calculate_ghia_error(&self, y_computed: &[f64], u_computed: &[f64]) -> PyResult<f64> {
        let reference = self.ghia_u_reference();
        if reference.is_empty() {
            return Err(PyValueError::new_err(format!(
                "Ghia Table I has no u-centerline column for Re={}",
                self.reynolds
            )));
        }
        let mut sum_sq_error = 0.0;
        let mut count = 0;

        for &(y_ghia, u_ghia) in &reference {
            let mut best_idx = 0;
            let mut best_dist = f64::MAX;
            for (i, &y) in y_computed.iter().enumerate() {
                let dist = (y - y_ghia).abs();
                if dist < best_dist {
                    best_dist = dist;
                    best_idx = i;
                }
            }

            let u_interp = *u_computed.get(best_idx).ok_or_else(|| {
                PyValueError::new_err("cavity result has no u-centerline samples")
            })?;
            sum_sq_error += (u_interp - u_ghia).powi(2);
            count += 1;
        }

        Ok((sum_sq_error / f64::from(count)).sqrt())
    }
}

// ── Result ───────────────────────────────────────────────────────────────────

/// Result from lid-driven cavity simulation
#[pyclass(name = "CavityResult2D", skip_from_py_object)]
#[derive(Debug, Clone)]
pub struct PyCavityResult2D {
    #[pyo3(get)]
    pub l2_error: f64,
    #[pyo3(get)]
    pub u_centerline: Vec<f64>,
    #[pyo3(get)]
    pub v_centerline: Vec<f64>,
    #[pyo3(get)]
    pub y_coords: Vec<f64>,
    #[pyo3(get)]
    pub x_coords: Vec<f64>,
    #[pyo3(get)]
    pub converged: bool,
}

#[pymethods]
impl PyCavityResult2D {
    fn __str__(&self) -> String {
        format!(
            "CavityResult2D(L2_error={:.4}, converged={})",
            self.l2_error, self.converged
        )
    }
}
