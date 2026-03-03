//! 3D Poiseuille flow (pipe flow) `PyO3` wrapper.

use pyo3::prelude::*;
use pyo3::exceptions::PyRuntimeError;

use cfd_core::physics::fluid::blood::CassonBlood;

/// 3D Poiseuille flow in circular pipe
///
/// **IMPORTANT - Sign Convention:**
/// This solver uses pressure **gradient** convention:
/// - `pressure_drop < 0` → forward flow (Q > 0)
/// - `pressure_drop > 0` → backward flow (Q < 0)
/// - dp/dx = `pressure_drop` / length = (`P_out` - `P_in`) / L
///
/// Analytical solution:
/// ```text
/// u(r) = u_max(1 - (r/R)²)
/// u_max = (R²/4μ)(dp/dx)
/// Q = (πR⁴/8μ)(dp/dx)
/// ```
///
/// **Note:** This is an analytical calculator using Hagen-Poiseuille formula
/// with constant viscosity evaluated at γ̇=100 s⁻¹. Grid parameters
/// (nr, ntheta, nz) are stored but not used in computation.
#[pyclass(name = "Poiseuille3DSolver")]
pub struct PyPoiseuille3DSolver {
    #[pyo3(get)]
    diameter: f64,
    #[pyo3(get)]
    length: f64,
    #[pyo3(get)]
    nr: usize,
    #[pyo3(get)]
    ntheta: usize,
    #[pyo3(get)]
    nz: usize,
}

#[pymethods]
impl PyPoiseuille3DSolver {
    #[new]
    fn new(diameter: f64, length: f64, nr: usize, ntheta: usize, nz: usize) -> Self {
        PyPoiseuille3DSolver {
            diameter,
            length,
            nr,
            ntheta,
            nz,
        }
    }

    /// Analytical maximum velocity
    fn analytical_max_velocity(&self, pressure_gradient: f64, viscosity: f64) -> f64 {
        let r = self.diameter / 2.0;
        (-pressure_gradient / (4.0 * viscosity)) * r.powi(2)
    }

    /// Analytical flow rate
    fn analytical_flow_rate(&self, pressure_gradient: f64, viscosity: f64) -> f64 {
        let r = self.diameter / 2.0;
        let pi = std::f64::consts::PI;
        (-pressure_gradient / (8.0 * viscosity)) * pi * r.powi(4)
    }

    fn __str__(&self) -> String {
        format!(
            "Poiseuille3DSolver(D={:.1} μm, L={:.2} mm, grid={}×{}×{})",
            self.diameter * 1e6,
            self.length * 1e3,
            self.nr,
            self.ntheta,
            self.nz
        )
    }

    /// Solve 3D pipe flow simulation
    ///
    /// # Arguments
    /// * `pressure_drop` - Pressure drop (Pa). Use NEGATIVE for forward flow!
    /// * `blood_type` - Blood model: "newtonian", "casson", or "`carreau_yasuda`"
    fn solve(&self, pressure_drop: f64, blood_type: &str) -> PyResult<PyPoiseuille3DResult> {
        use cfd_core::physics::fluid::blood::CarreauYasudaBlood;

        let rho = 1060.0;
        let gamma_dot = 100.0;

        let mu = match blood_type {
            "newtonian" => 0.0035,
            "casson" => {
                let fluid = CassonBlood::<f64>::normal_blood();
                fluid.apparent_viscosity(gamma_dot)
            },
            "carreau_yasuda" | "carreau_yasuda_blood" => {
                let fluid = CarreauYasudaBlood::<f64>::normal_blood();
                fluid.apparent_viscosity(gamma_dot)
            },
            _ => {
                return Err(PyRuntimeError::new_err(
                    format!("Unknown blood type '{blood_type}'. Use: newtonian, casson, or carreau_yasuda")
                ));
            }
        };

        let dp_dx = pressure_drop / self.length;
        let u_max = self.analytical_max_velocity(dp_dx, mu);
        let q = self.analytical_flow_rate(dp_dx, mu);

        Ok(PyPoiseuille3DResult {
            max_velocity: u_max,
            flow_rate: q,
            reynolds_number: (rho * u_max * self.diameter) / mu,
            wall_shear_stress: (dp_dx * self.diameter) / 4.0,
        })
    }
}

#[pyclass(name = "Poiseuille3DResult")]
pub struct PyPoiseuille3DResult {
    #[pyo3(get)]
    pub max_velocity: f64,
    #[pyo3(get)]
    pub flow_rate: f64,
    #[pyo3(get)]
    pub reynolds_number: f64,
    #[pyo3(get)]
    pub wall_shear_stress: f64,
}

#[pymethods]
impl PyPoiseuille3DResult {
    fn __str__(&self) -> String {
        format!(
            "Poiseuille3DResult(u_max={:.3e} m/s, Q={:.3e} m³/s, Re={:.1})",
            self.max_velocity, self.flow_rate, self.reynolds_number
        )
    }
}
