//! 2D Bifurcation and Trifurcation solver `PyO3` wrappers.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};
use cfd_validation::benchmarks::{Benchmark, BenchmarkConfig, TrifurcationFlow};

// ── Trifurcation Solver ──────────────────────────────────────────────────────

/// 2D Trifurcation flow solver
#[pyclass(name = "TrifurcationSolver2D", skip_from_py_object)]
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
        PyTrifurcationSolver2D {
            width,
            length,
            angle,
            nx,
        }
    }

    /// Solve 2D trifurcation simulation
    fn solve(
        &self,
        py: Python<'_>,
        _flow_rate: f64,
        _blood_type: &str,
    ) -> PyResult<PyTrifurcationResult2D> {
        let width = self.width;
        let length = self.length;
        let angle = self.angle;
        let resolution = self.nx;

        py.detach(move || {
            let bench = TrifurcationFlow::new(width, length, angle);
            let config = BenchmarkConfig {
                resolution,
                max_iterations: 100,
                ..Default::default()
            };

            let result = bench
                .run(&config)
                .map_err(|e| format!("Benchmark error: {e}"))?;

            Ok::<_, String>(PyTrifurcationResult2D {
                execution_time: result.execution_time,
                mass_conservation_error: result.convergence.last().copied().unwrap_or(0.0),
            })
        })
        .map_err(PyRuntimeError::new_err)
    }
}

/// Result from 2D trifurcation simulation
#[pyclass(name = "TrifurcationResult2D", skip_from_py_object)]
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

// ── Bifurcation Solver ───────────────────────────────────────────────────────

/// 2D bifurcation flow solver for branching microfluidic channels.
///
/// Solves Navier-Stokes equations on a staggered grid with geometry mask
/// defining the Y-shaped bifurcation domain. Validates mass conservation
/// at the junction and computes flow split between daughter branches.
///
/// # Physics
///
/// Mass conservation: `Q_parent` = `Q_daughter1` + `Q_daughter2`
/// Murray's law (optimal branching): `r_p^3` = `r_d1^3` + `r_d2^3`
#[pyclass(name = "BifurcationSolver2D", skip_from_py_object)]
pub struct PyBifurcationSolver2D {
    #[pyo3(get)]
    parent_width: f64,
    #[pyo3(get)]
    parent_length: f64,
    #[pyo3(get)]
    daughter_width: f64,
    #[pyo3(get)]
    daughter_length: f64,
    #[pyo3(get)]
    angle: f64,
    #[pyo3(get)]
    nx: usize,
    #[pyo3(get)]
    ny: usize,
}

#[pymethods]
impl PyBifurcationSolver2D {
    /// Create symmetric 2D bifurcation solver.
    #[new]
    #[pyo3(signature = (parent_width, parent_length, daughter_width, daughter_length, angle=0.5, nx=50, ny=30))]
    fn new(
        parent_width: f64,
        parent_length: f64,
        daughter_width: f64,
        daughter_length: f64,
        angle: f64,
        nx: usize,
        ny: usize,
    ) -> Self {
        PyBifurcationSolver2D {
            parent_width,
            parent_length,
            daughter_width,
            daughter_length,
            angle,
            nx,
            ny,
        }
    }

    /// Murray's law deviation for symmetric bifurcation.
    ///
    /// Returns relative deviation from Murray's law: |`w_p^3` - 2*`w_d^3`| / `w_p^3`
    fn murray_law_deviation(&self) -> f64 {
        let wp3 = self.parent_width.powi(3);
        let wd3 = 2.0 * self.daughter_width.powi(3);
        (wp3 - wd3).abs() / wp3
    }

    /// Solve 2D bifurcation flow.
    fn solve(
        &self,
        py: Python<'_>,
        inlet_velocity: f64,
        blood_type: &str,
    ) -> PyResult<PyBifurcationResult2D> {
        use cfd_2d::solvers::bifurcation_flow::{BifurcationGeometry, BifurcationSolver2D};
        use cfd_2d::solvers::ns_fvm::{BloodModel, SIMPLEConfig};

        let parent_width = self.parent_width;
        let parent_length = self.parent_length;
        let daughter_width = self.daughter_width;
        let daughter_length = self.daughter_length;
        let angle = self.angle;
        let nx = self.nx;
        let ny = self.ny;
        let blood_type = blood_type.to_owned();

        py.detach(move || {
            let geom = BifurcationGeometry::new_symmetric(
                parent_width,
                parent_length,
                daughter_width,
                daughter_length,
                angle,
            );

            let blood = match blood_type.as_str() {
                "casson" => BloodModel::Casson(CassonBlood::normal_blood()),
                "carreau_yasuda" => BloodModel::CarreauYasuda(CarreauYasudaBlood::normal_blood()),
                _ => BloodModel::Newtonian(0.0035),
            };

            let density = 1060.0;
            let mut config = SIMPLEConfig::default();
            config.max_iterations = 5000;
            config.tolerance = 1e-5;
            config.alpha_u = 0.5;
            config.alpha_p = 0.2;

            let mut solver = BifurcationSolver2D::new(geom, blood, density, nx, ny, config);

            let sol = solver
                .solve(inlet_velocity)
                .map_err(|e| format!("Bifurcation solver error: {e}"))?;

            let flow_split = if sol.q_parent.abs() > 1e-30 {
                sol.q_daughter1 / sol.q_parent
            } else {
                0.5
            };

            Ok::<_, String>(PyBifurcationResult2D {
                q_parent: sol.q_parent,
                q_daughter1: sol.q_daughter1,
                q_daughter2: sol.q_daughter2,
                mass_balance_error: sol.mass_balance_error,
                flow_split_ratio: flow_split,
            })
        })
        .map_err(PyRuntimeError::new_err)
    }

    fn __str__(&self) -> String {
        format!(
            "BifurcationSolver2D(w_p={:.1} μm, w_d={:.1} μm, θ={:.1}°, grid={}×{})",
            self.parent_width * 1e6,
            self.daughter_width * 1e6,
            self.angle.to_degrees(),
            self.nx,
            self.ny
        )
    }
}

// ── Bifurcation Result ───────────────────────────────────────────────────────

/// Result from 2D bifurcation simulation
#[pyclass(name = "BifurcationResult2D", skip_from_py_object)]
#[derive(Debug, Clone, Copy)]
pub struct PyBifurcationResult2D {
    #[pyo3(get)]
    pub q_parent: f64,
    #[pyo3(get)]
    pub q_daughter1: f64,
    #[pyo3(get)]
    pub q_daughter2: f64,
    #[pyo3(get)]
    pub mass_balance_error: f64,
    #[pyo3(get)]
    pub flow_split_ratio: f64,
}

#[pymethods]
impl PyBifurcationResult2D {
    fn __str__(&self) -> String {
        format!(
            "BifurcationResult2D(Q_p={:.2e}, split={:.3}, mass_err={:.2e})",
            self.q_parent, self.flow_split_ratio, self.mass_balance_error
        )
    }
}
