//! 2D CFD solver PyO3 wrappers
//!
//! This module exposes 2D finite volume solvers for validation against
//! Python CFD packages like FEniCS, SfePy, etc.

use nalgebra::DMatrix;
use numpy::PyArray2;
use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;


use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::simplec_pimple::solver::SimplecPimpleSolver;
use cfd_2d::simplec_pimple::config::SimplecPimpleConfig;
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};
use cfd_core::physics::fluid::Fluid;
use cfd_validation::benchmarks::{Benchmark, BenchmarkConfig, BifurcationFlow, TrifurcationFlow, VenturiFlow};
use cfd_validation::geometry::{Venturi2D, Point2D, Geometry2D};

// ============================================================================
// 2D Poiseuille Flow Solver (Analytical Validation)
// ============================================================================

/// 2D Poiseuille flow solver for channel flow validation
///
/// Solves laminar flow between parallel plates with analytical solution:
///
/// ```text
/// u(y) = (1/(2μ)) * (dp/dx) * y(H - y)
/// u_max = (H²/8μ) * (dp/dx)
/// Q = (H³W/12μ) * (dp/dx)
/// ```
///
/// Where:
/// - H: channel height
/// - W: channel width
/// - μ: dynamic viscosity
/// - dp/dx: pressure gradient
#[pyclass]
pub struct PyPoiseuille2DSolver {
    #[pyo3(get)]
    height: f64,
    #[pyo3(get)]
    width: f64,
    #[pyo3(get)]
    length: f64,
    #[pyo3(get)]
    nx: usize,
    #[pyo3(get)]
    ny: usize,
}

#[pymethods]
impl PyPoiseuille2DSolver {
    /// Create new 2D Poiseuille solver
    #[new]
    fn new(height: f64, width: f64, length: f64, nx: usize, ny: usize) -> Self {
        PyPoiseuille2DSolver {
            height,
            width,
            length,
            nx,
            ny,
        }
    }

    /// Compute analytical velocity profile
    fn analytical_velocity_profile<'py>(
        &self,
        py: Python<'py>,
        pressure_gradient: f64,
        viscosity: f64,
    ) -> PyResult<Bound<'py, PyArray2<f64>>> {
        let mut u = DMatrix::<f64>::zeros(self.ny, self.nx);
        let dy = self.height / (self.ny - 1) as f64;

        for j in 0..self.ny {
            let y = j as f64 * dy;
            let u_val = (-pressure_gradient / (2.0 * viscosity)) * y * (self.height - y);
            for i in 0..self.nx {
                u[(j, i)] = u_val;
            }
        }

        // Convert DMatrix to ndarray Array2
        let array = ndarray::Array2::from_shape_fn((self.ny, self.nx), |(j, i)| u[(j, i)]);
        Ok(PyArray2::from_owned_array_bound(py, array))
    }

    /// Compute analytical maximum velocity
    fn analytical_max_velocity(&self, pressure_gradient: f64, viscosity: f64) -> f64 {
        (-pressure_gradient / (8.0 * viscosity)) * self.height.powi(2)
    }

    /// Compute analytical flow rate
    fn analytical_flow_rate(&self, pressure_gradient: f64, viscosity: f64) -> f64 {
        (-pressure_gradient / (12.0 * viscosity)) * self.height.powi(3) * self.width
    }

    /// Get grid spacing
    fn grid_spacing(&self) -> (f64, f64) {
        let dx = self.length / (self.nx - 1) as f64;
        let dy = self.height / (self.ny - 1) as f64;
        (dx, dy)
    }

    fn __str__(&self) -> String {
        format!(
            "Poiseuille2DSolver(H={:.1} μm, L={:.2} mm, grid={}×{})",
            self.height * 1e6,
            self.length * 1e3,
            self.nx,
            self.ny
        )
    }

    /// Solve Pooiseuille flow simulation
    fn solve(&self, pressure_drop: f64, blood_type: &str) -> PyResult<PyPoiseuille2DResult> {
        let grid = StructuredGrid2D::new(
            self.nx,
            self.ny,
            0.0,
            self.length,
            0.0,
            self.height,
        ).map_err(|e| PyRuntimeError::new_err(format!("Grid error: {}", e)))?;

        let mut fields = SimulationFields::new(self.nx, self.ny);
        let rho = 1060.0;
        let pressure_gradient = pressure_drop / self.length;

        // Set up boundary conditions and initial viscosity
        let viscosity = match blood_type {
            "casson" => CassonBlood::<f64>::normal_blood().apparent_viscosity(100.0),
            "carreau_yasuda" => CarreauYasudaBlood::<f64>::normal_blood().apparent_viscosity(100.0),
            _ => 0.0035, // Default water-like
        };

        for i in 0..self.nx {
            for j in 0..self.ny {
                fields.viscosity.set(i, j, viscosity);
                fields.mask.set(i, j, true);
            }
        }

        let mut solver = SimplecPimpleSolver::new(grid, SimplecPimpleConfig::default())
            .map_err(|e| PyRuntimeError::new_err(format!("Solver error: {}", e)))?;

        // Simplified steady-state solve for Poiseuille
        for _ in 0..100 {
            solver.solve_time_step(&mut fields, 0.01, 0.0, rho)
                .map_err(|e| PyRuntimeError::new_err(format!("Step error: {}", e)))?;
        }

        let max_u = fields.u.data().iter().copied().fold(0.0, f64::max);
        let flow_rate = (fields.u.at(self.nx / 2, self.ny / 2) * self.height * self.width) / 1.5; // Rough estimate

        Ok(PyPoiseuille2DResult {
            max_velocity: max_u,
            flow_rate,
            reynolds_number: (rho * max_u * self.height) / viscosity,
            pressure_drop,
            wall_shear_stress: (pressure_gradient * self.height) / 2.0,
        })
    }
}

// ============================================================================
// 2D Venturi Flow Solver
// ============================================================================

/// 2D Venturi throat flow solver with Bernoulli validation
#[pyclass]
pub struct PyVenturiSolver2D {
    #[pyo3(get)]
    w_inlet: f64,
    #[pyo3(get)]
    w_throat: f64,
    #[pyo3(get)]
    l_inlet: f64,
    #[pyo3(get)]
    l_converge: f64,
    #[pyo3(get)]
    l_throat: f64,
    #[pyo3(get)]
    l_diverge: f64,
    #[pyo3(get)]
    nx: usize,
    #[pyo3(get)]
    ny: usize,
}

#[pymethods]
impl PyVenturiSolver2D {
    #[new]
    #[pyo3(signature = (w_inlet, w_throat, l_inlet, l_converge, l_throat, l_diverge, nx=200, ny=100))]
    fn new(
        w_inlet: f64,
        w_throat: f64,
        l_inlet: f64,
        l_converge: f64,
        l_throat: f64,
        l_diverge: f64,
        nx: usize,
        ny: usize,
    ) -> Self {
        PyVenturiSolver2D {
            w_inlet,
            w_throat,
            l_inlet,
            l_converge,
            l_throat,
            l_diverge,
            nx,
            ny,
        }
    }

    /// Create ISO 5167 standard Venturi (area ratio 0.5)
    #[staticmethod]
    fn iso_5167_standard(nx: usize, ny: usize) -> Self {
        Self {
            w_inlet: 10e-3,
            w_throat: 7.07e-3, // √0.5 ratio
            l_inlet: 1e-3,
            l_converge: 1e-3,
            l_throat: 2e-3,
            l_diverge: 5e-3,
            nx,
            ny,
        }
    }

    /// Compute analytical pressure coefficient in throat (Bernoulli)
    fn pressure_coefficient_analytical(&self) -> f64 {
        let area_ratio = self.w_throat / self.w_inlet;
        1.0 - area_ratio.powi(2)
    }

    /// Get area ratio (throat/inlet)
    fn area_ratio(&self) -> f64 {
        self.w_throat / self.w_inlet
    }

    /// Get total length
    fn total_length(&self) -> f64 {
        self.l_inlet + self.l_converge + self.l_throat + self.l_diverge
    }

    fn __str__(&self) -> String {
        format!(
            "VenturiSolver2D(β={:.3}, L={:.1} mm, grid={}×{})",
            self.area_ratio(),
            self.total_length() * 1e3,
            self.nx,
            self.ny
        )
    }

    /// Solve Venturi flow simulation
    fn solve(&self, inlet_velocity: f64, _blood_type: &str) -> PyResult<PyVenturiResult2D> {
        let bench = VenturiFlow::new(self.w_inlet, self.w_throat);
        let config = BenchmarkConfig {
            resolution: self.nx,
            max_iterations: 100,
            ..Default::default()
        };

        let result = bench.run(&config)
            .map_err(|e| PyRuntimeError::new_err(format!("Benchmark error: {}", e)))?;

        let cp = result.metrics.get("pressure_coefficient").cloned().unwrap_or(0.0);
        let recovery = result.metrics.get("pressure_recovery").cloned().unwrap_or(0.0);

        Ok(PyVenturiResult2D {
            cp_throat: cp,
            pressure_recovery: recovery,
            velocity_ratio: self.w_inlet / self.w_throat,
            mass_conservation_error: result.convergence.last().copied().unwrap_or(0.0),
        })
    }
}

// ============================================================================
// 2D Trifurcation Solver
// ============================================================================

/// 2D Trifurcation flow solver
#[pyclass]
pub struct PyTrifurcationSolver2D {
    #[pyo3(get)]
    pub width: f64,
    #[pyo3(get)]
    pub length: f64,
    #[pyo3(get)]
    pub angle: f64,
    #[pyo3(get)]
    pub nx: usize,
}

#[pymethods]
impl PyTrifurcationSolver2D {
    #[new]
    #[pyo3(signature = (width, length, angle, nx=128))]
    fn new(width: f64, length: f64, angle: f64, nx: usize) -> Self {
        PyTrifurcationSolver2D { width, length, angle, nx }
    }

    /// Solve 2D trifurcation simulation
    fn solve(&self, flow_rate: f64, _blood_type: &str) -> PyResult<PyTrifurcationResult2D> {
        let bench = TrifurcationFlow::new(self.width, self.length, self.angle);
        let config = BenchmarkConfig {
            resolution: self.nx,
            max_iterations: 100,
            ..Default::default()
        };

        let result = bench.run(&config)
            .map_err(|e| PyRuntimeError::new_err(format!("Benchmark error: {}", e)))?;

        Ok(PyTrifurcationResult2D {
            execution_time: result.execution_time,
            mass_conservation_error: result.convergence.last().copied().unwrap_or(0.0),
        })
    }
}

/// Result from 2D trifurcation simulation
#[pyclass]
#[derive(Debug, Clone, Copy)]
pub struct PyTrifurcationResult2D {
    #[pyo3(get)]
    pub execution_time: f64,
    #[pyo3(get)]
    pub mass_conservation_error: f64,
}

#[pymethods]
impl PyTrifurcationResult2D {
    fn __str__(&self) -> String {
        format!(
            "TrifurcationResult2D(time={:.2}s, mass_err={:.2e})",
            self.execution_time, self.mass_conservation_error
        )
    }
}

// ============================================================================
// 2D Result Types
// ============================================================================

/// Result from 2D Poiseuille flow simulation
#[pyclass]
pub struct PyPoiseuille2DResult {
    #[pyo3(get)]
    pub max_velocity: f64,
    #[pyo3(get)]
    pub flow_rate: f64,
    #[pyo3(get)]
    pub reynolds_number: f64,
    #[pyo3(get)]
    pub pressure_drop: f64,
    #[pyo3(get)]
    pub wall_shear_stress: f64,
}

#[pymethods]
impl PyPoiseuille2DResult {
    fn __str__(&self) -> String {
        format!(
            "Poiseuille2DResult(u_max={:.3e} m/s, Q={:.3e} m³/s, Re={:.1})",
            self.max_velocity, self.flow_rate, self.reynolds_number
        )
    }
}

/// Result from 2D Venturi flow simulation
#[pyclass]
pub struct PyVenturiResult2D {
    #[pyo3(get)]
    pub cp_throat: f64,
    #[pyo3(get)]
    pub pressure_recovery: f64,
    #[pyo3(get)]
    pub velocity_ratio: f64,
    #[pyo3(get)]
    pub mass_conservation_error: f64,
}

#[pymethods]
impl PyVenturiResult2D {
    fn __str__(&self) -> String {
        format!(
            "VenturiResult2D(Cp_throat={:.3}, recovery={:.1}%, mass_error={:.2e})",
            self.cp_throat,
            self.pressure_recovery * 100.0,
            self.mass_conservation_error
        )
    }
}
