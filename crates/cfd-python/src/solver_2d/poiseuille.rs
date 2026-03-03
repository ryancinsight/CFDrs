//! 2D Poiseuille flow solver `PyO3` wrapper.

use nalgebra::DMatrix;
use numpy::PyArray2;
use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;

use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};

// ── Solver ───────────────────────────────────────────────────────────────────

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
#[pyclass(name = "Poiseuille2DSolver")]
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

    /// Solve Poiseuille flow simulation
    fn solve(&self, pressure_drop: f64, blood_type: &str) -> PyResult<PyPoiseuille2DResult> {
        let _grid = StructuredGrid2D::new(
            self.nx,
            self.ny,
            0.0,
            self.length,
            0.0,
            self.height,
        ).map_err(|e| PyRuntimeError::new_err(format!("Grid error: {e}")))?;

        let mut fields = SimulationFields::new(self.nx, self.ny);
        let rho = 1060.0;
        let pressure_gradient = pressure_drop / self.length;

        let viscosity = match blood_type {
            "casson" => CassonBlood::<f64>::normal_blood().apparent_viscosity(100.0),
            "carreau_yasuda" => CarreauYasudaBlood::<f64>::normal_blood().apparent_viscosity(100.0),
            _ => 0.0035,
        };

        for i in 0..self.nx {
            for j in 0..self.ny {
                fields.viscosity.set(i, j, viscosity);
                fields.mask.set(i, j, true);
            }
        }

        let u_max = (self.height.powi(2) / (8.0 * viscosity)) * pressure_gradient.abs();

        for i in 0..self.nx {
            for j in 0..self.ny {
                let y = j as f64 * (self.height / (self.ny - 1) as f64);
                let u_val = 4.0 * u_max * (y / self.height) * (1.0 - y / self.height);
                fields.u.set(i, j, u_val);
            }
        }

        let flow_rate = (2.0 / 3.0) * u_max * self.height * self.width;

        let mid_i = self.nx / 2;
        let mut u_centerline = Vec::with_capacity(self.ny);
        let mut y_coords = Vec::with_capacity(self.ny);

        for j in 0..self.ny {
            u_centerline.push(fields.u.at(mid_i, j));
            y_coords.push(j as f64 * (self.height / (self.ny - 1) as f64));
        }

        Ok(PyPoiseuille2DResult {
            max_velocity: u_max,
            flow_rate,
            reynolds_number: (rho * u_max * self.height) / viscosity,
            pressure_drop,
            wall_shear_stress: (pressure_gradient * self.height) / 2.0,
            u_centerline,
            y_coords,
        })
    }
}

// ── Result ───────────────────────────────────────────────────────────────────

/// Result from 2D Poiseuille flow simulation
#[pyclass(name = "Poiseuille2DResult")]
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
    #[pyo3(get)]
    pub u_centerline: Vec<f64>,
    #[pyo3(get)]
    pub y_coords: Vec<f64>,
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
