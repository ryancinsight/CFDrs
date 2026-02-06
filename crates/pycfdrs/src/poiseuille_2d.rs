//! Python bindings for 2D Poiseuille flow solver with non-Newtonian blood

use cfd_2d::solvers::{BloodModel as RustBloodModel, PoiseuilleConfig, PoiseuilleFlow2D};
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};
use pyo3::prelude::*;

use crate::blood::{PyCarreauYasudaBlood, PyCassonBlood};

/// Configuration for 2D Poiseuille flow solver
#[pyclass(name = "PoiseuilleConfig2D")]
#[derive(Clone)]
pub struct PyPoiseuilleConfig {
    /// Channel height [m]
    pub height: f64,
    /// Channel width [m]
    pub width: f64,
    /// Channel length [m]
    pub length: f64,
    /// Number of grid points in y-direction
    pub ny: usize,
    /// Pressure gradient dP/dx [Pa/m]
    pub pressure_gradient: f64,
    /// Convergence tolerance
    pub tolerance: f64,
    /// Maximum iterations
    pub max_iterations: usize,
    /// Relaxation factor (0 < alpha ≤ 1)
    pub relaxation_factor: f64,
}

#[pymethods]
impl PyPoiseuilleConfig {
    /// Create new Poiseuille configuration
    ///
    /// # Arguments
    ///
    /// * `height` - Channel height [m]
    /// * `width` - Channel width [m] (for flow rate calculation)
    /// * `length` - Channel length [m] (for reference)
    /// * `ny` - Number of grid points in y-direction
    /// * `pressure_gradient` - Pressure gradient dP/dx [Pa/m]
    /// * `tolerance` - Convergence tolerance (default: 1e-6)
    /// * `max_iterations` - Maximum iterations (default: 1000)
    /// * `relaxation_factor` - Relaxation factor (default: 0.5)
    ///
    /// # Returns
    ///
    /// Configuration object
    ///
    /// # Example
    ///
    /// ```python
    /// config = PoiseuilleConfig2D(
    ///     height=0.001,      # 1 mm
    ///     width=0.01,        # 10 mm
    ///     length=0.05,       # 50 mm
    ///     ny=101,
    ///     pressure_gradient=1000.0,  # Pa/m
    ///     tolerance=1e-8,
    ///     max_iterations=1000,
    ///     relaxation_factor=0.5
    /// )
    /// ```
    #[new]
    #[pyo3(signature = (height, width, length, ny, pressure_gradient, tolerance=1e-6, max_iterations=1000, relaxation_factor=0.5))]
    fn new(
        height: f64,
        width: f64,
        length: f64,
        ny: usize,
        pressure_gradient: f64,
        tolerance: f64,
        max_iterations: usize,
        relaxation_factor: f64,
    ) -> Self {
        PyPoiseuilleConfig {
            height,
            width,
            length,
            ny,
            pressure_gradient,
            tolerance,
            max_iterations,
            relaxation_factor,
        }
    }

    fn __repr__(&self) -> String {
        format!(
            "PoiseuilleConfig2D(height={:.6e}, width={:.6e}, ny={}, dP/dx={:.3e} Pa/m)",
            self.height, self.width, self.ny, self.pressure_gradient
        )
    }
}

/// Result from 2D Poiseuille flow solver
#[pyclass(name = "PoiseuilleResult2D")]
#[derive(Clone)]
pub struct PyPoiseuilleResult {
    /// Y coordinates [m]
    pub y_coords: Vec<f64>,
    /// Velocity profile u(y) [m/s]
    pub velocity: Vec<f64>,
    /// Shear rate profile γ(y) [1/s]
    pub shear_rate: Vec<f64>,
    /// Viscosity profile μ(y) [Pa·s]
    pub viscosity: Vec<f64>,
    /// Flow rate Q [m³/s]
    pub flow_rate: f64,
    /// Wall shear stress τ_w [Pa]
    pub wall_shear_stress: f64,
    /// Number of iterations to converge
    pub iterations: usize,
}

#[pymethods]
impl PyPoiseuilleResult {
    fn __repr__(&self) -> String {
        format!(
            "PoiseuilleResult2D(Q={:.6e} m³/s, τ_w={:.3e} Pa, iterations={})",
            self.flow_rate, self.wall_shear_stress, self.iterations
        )
    }

    /// Get y coordinates as Python list
    #[getter]
    fn get_y_coords(&self) -> Vec<f64> {
        self.y_coords.clone()
    }

    /// Get velocity profile as Python list
    #[getter]
    fn get_velocity(&self) -> Vec<f64> {
        self.velocity.clone()
    }

    /// Get shear rate profile as Python list
    #[getter]
    fn get_shear_rate(&self) -> Vec<f64> {
        self.shear_rate.clone()
    }

    /// Get viscosity profile as Python list
    #[getter]
    fn get_viscosity(&self) -> Vec<f64> {
        self.viscosity.clone()
    }

    /// Get flow rate [m³/s]
    #[getter]
    fn get_flow_rate(&self) -> f64 {
        self.flow_rate
    }

    /// Get wall shear stress [Pa]
    #[getter]
    fn get_wall_shear_stress(&self) -> f64 {
        self.wall_shear_stress
    }

    /// Get number of iterations
    #[getter]
    fn get_iterations(&self) -> usize {
        self.iterations
    }
}

/// 2D Poiseuille flow solver with non-Newtonian blood rheology
#[pyclass(name = "PoiseuilleSolver2D")]
pub struct PyPoiseuilleSolver {
    config: PyPoiseuilleConfig,
}

#[pymethods]
impl PyPoiseuilleSolver {
    /// Create new Poiseuille solver
    ///
    /// # Arguments
    ///
    /// * `config` - Solver configuration
    ///
    /// # Example
    ///
    /// ```python
    /// config = PoiseuilleConfig2D(height=0.001, width=0.01, length=0.05, ny=101, pressure_gradient=1000.0)
    /// solver = PoiseuilleSolver2D(config)
    /// ```
    #[new]
    fn new(config: PyPoiseuilleConfig) -> Self {
        PyPoiseuilleSolver { config }
    }

    /// Solve Poiseuille flow with specified blood model
    ///
    /// # Arguments
    ///
    /// * `blood` - Blood rheology model (Casson or Carreau-Yasuda)
    ///
    /// # Returns
    ///
    /// PoiseuilleResult2D with velocity, shear rate, viscosity profiles and flow rate
    ///
    /// # Example
    ///
    /// ```python
    /// # Using Casson blood model
    /// blood = CassonBlood.normal_blood()
    /// result = solver.solve(blood)
    /// print(f"Flow rate: {result.flow_rate:.6e} m³/s")
    /// print(f"Wall shear stress: {result.wall_shear_stress:.3e} Pa")
    /// print(f"Converged in {result.iterations} iterations")
    ///
    /// # Access profiles
    /// import matplotlib.pyplot as plt
    /// plt.plot(result.y_coords, result.velocity)
    /// plt.xlabel('y [m]')
    /// plt.ylabel('u [m/s]')
    /// plt.show()
    /// ```
    fn solve(&self, blood: &Bound<'_, PyAny>) -> PyResult<PyPoiseuilleResult> {
        // Convert Python blood model to Rust blood model
        // Try to downcast to each blood model type
        let rust_blood_model = if blood.downcast::<PyCassonBlood>().is_ok() {
            // Use the normal_blood() constructor which has all fields set
            RustBloodModel::Casson(CassonBlood::<f64>::normal_blood())
        } else if blood.downcast::<PyCarreauYasudaBlood>().is_ok() {
            // Use the normal_blood() constructor which has all fields set
            RustBloodModel::CarreauYasuda(CarreauYasudaBlood::<f64>::normal_blood())
        } else {
            return Err(pyo3::exceptions::PyTypeError::new_err(
                "Blood model must be CassonBlood or CarreauYasudaBlood",
            ));
        };

        // Convert Python config to Rust config
        let mut rust_config = PoiseuilleConfig::<f64>::default();
        rust_config.height = self.config.height;
        rust_config.width = self.config.width;
        rust_config.length = self.config.length;
        rust_config.ny = self.config.ny;
        rust_config.pressure_gradient = self.config.pressure_gradient;
        rust_config.tolerance = self.config.tolerance;
        rust_config.max_iterations = self.config.max_iterations;
        rust_config.relaxation_factor = self.config.relaxation_factor;

        // Create and run solver
        let mut solver = PoiseuilleFlow2D::new(rust_config, rust_blood_model);

        let iterations = solver.solve().map_err(|e| {
            pyo3::exceptions::PyRuntimeError::new_err(format!("Solver failed: {}", e))
        })?;

        // Extract results
        let y_coords = solver.y_coordinates().to_vec();
        let velocity = solver.velocity_profile().to_vec();
        let shear_rate = solver.shear_rate_profile().to_vec();
        let viscosity = solver.viscosity_profile().to_vec();
        let flow_rate = solver.flow_rate();
        let wall_shear_stress = solver.wall_shear_stress();

        Ok(PyPoiseuilleResult {
            y_coords,
            velocity,
            shear_rate,
            viscosity,
            flow_rate,
            wall_shear_stress,
            iterations,
        })
    }

    fn __repr__(&self) -> String {
        format!("PoiseuilleSolver2D({})", self.config.__repr__())
    }
}
