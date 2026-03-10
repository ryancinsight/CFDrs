//! 2D and 1D Venturi solver `PyO3` wrappers.

use pyo3::exceptions::PyRuntimeError;
use pyo3::prelude::*;

use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};

// ── 2D Venturi Solver ────────────────────────────────────────────────────────

/// 2D Venturi throat flow solver with Bernoulli validation
#[pyclass(name = "VenturiSolver2D")]
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
            w_throat: 7.07e-3,
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

    /// Solve Venturi flow simulation using FVM Navier-Stokes solver.
    fn solve(&self, inlet_velocity: f64, blood_type: &str) -> PyResult<PyVenturiResult2D> {
        use cfd_2d::solvers::ns_fvm::BloodModel;
        use cfd_2d::solvers::venturi_flow::{VenturiGeometry as VGeom, VenturiSolver2D as VSolver};

        let geom = VGeom::new(
            self.w_inlet,
            self.w_throat,
            self.l_inlet,
            self.l_converge,
            self.l_throat,
            self.l_diverge,
            1.0e-3,
        );

        let blood = match blood_type {
            "casson" => BloodModel::Casson(CassonBlood::normal_blood()),
            "carreau_yasuda" => BloodModel::CarreauYasuda(CarreauYasudaBlood::normal_blood()),
            _ => BloodModel::Newtonian(0.0035),
        };

        let density = 1060.0;
        let mut solver = VSolver::new(geom, blood, density, self.nx, self.ny);

        let sol = solver
            .solve(inlet_velocity)
            .map_err(|e| PyRuntimeError::new_err(format!("Venturi solver error: {e}")))?;

        Ok(PyVenturiResult2D {
            cp_throat: sol.cp_throat,
            pressure_recovery: sol.cp_recovery,
            velocity_ratio: sol.u_throat / sol.u_inlet.max(1e-30),
            mass_conservation_error: ((sol.u_inlet - sol.u_outlet) / sol.u_inlet.max(1e-30)).abs(),
        })
    }
}

// ── 2D Venturi Result ────────────────────────────────────────────────────────

/// Result from 2D Venturi flow simulation
#[pyclass(name = "VenturiResult2D")]
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

// ── 1D Venturi Resistance Solver ─────────────────────────────────────────────

/// 1D Venturi tube resistance model.
///
/// Computes total pressure drop through a converging-diverging channel:
/// - Contraction loss via discharge coefficient (ISO 5167)
/// - Throat friction (Darcy)
/// - Expansion recovery loss (Borda-Carnot)
///
/// # Physics
///
/// `ΔP_total` = `ΔP_contraction` + `ΔP_friction` + `ΔP_expansion`
/// `ΔP_contraction` = (`ρV_t²/2)(1` - `β⁴)/C_d²` where β = `D_t/D_1`
/// `ΔP_friction` = f*(`L_t/D_t`)*(ρV_t²/2)
/// `ΔP_expansion` = `K_exp`*(ρ/2)*(`V_t` - `V_3)²`
#[pyclass(name = "VenturiSolver1D")]
pub struct PyVenturiSolver1D {
    #[pyo3(get)]
    inlet_diameter: f64,
    #[pyo3(get)]
    throat_diameter: f64,
    #[pyo3(get)]
    throat_length: f64,
    #[pyo3(get)]
    total_length: f64,
}

#[pymethods]
impl PyVenturiSolver1D {
    #[new]
    #[pyo3(signature = (inlet_diameter, throat_diameter, throat_length, total_length=0.01))]
    fn new(
        inlet_diameter: f64,
        throat_diameter: f64,
        throat_length: f64,
        total_length: f64,
    ) -> Self {
        PyVenturiSolver1D {
            inlet_diameter,
            throat_diameter,
            throat_length,
            total_length,
        }
    }

    /// Area ratio β = `D_throat` / `D_inlet`
    fn beta(&self) -> f64 {
        self.throat_diameter / self.inlet_diameter
    }

    /// Solve Venturi resistance for given flow conditions.
    fn solve(&self, velocity: f64, blood_type: &str) -> PyResult<PyVenturiResult1D> {
        use cfd_1d::{FlowConditions, ResistanceModel, VenturiModel};

        let model = VenturiModel::symmetric(
            self.inlet_diameter,
            self.throat_diameter,
            self.throat_length,
            self.total_length,
        );

        let density = 1060.0;
        let dh = self.inlet_diameter;

        let mu = match blood_type {
            "casson" => {
                let blood = CassonBlood::<f64>::normal_blood();
                blood.apparent_viscosity(velocity / dh * 8.0)
            }
            "carreau_yasuda" => {
                let blood = CarreauYasudaBlood::<f64>::normal_blood();
                blood.apparent_viscosity(velocity / dh * 8.0)
            }
            _ => 0.0035,
        };

        let re = density * velocity * dh / mu;
        let mut conditions = FlowConditions::new(velocity);
        conditions.reynolds_number = Some(re);

        let fluid = cfd_core::physics::fluid::ConstantPropertyFluid::new(
            "blood".to_string(),
            density,
            mu,
            3617.0,
            0.52,
            1570.0,
        );

        let resistance = model
            .calculate_resistance(&fluid, &conditions)
            .map_err(|e| PyRuntimeError::new_err(format!("Venturi solver error: {e}")))?;

        let area = std::f64::consts::PI / 4.0 * self.inlet_diameter.powi(2);
        let flow_rate = velocity * area;
        let pressure_drop = resistance * flow_rate;

        let beta = self.beta();
        let dp_bernoulli = 0.5 * density * velocity.powi(2) * (1.0 / beta.powi(4) - 1.0);

        Ok(PyVenturiResult1D {
            pressure_drop,
            resistance,
            dp_bernoulli,
            reynolds_number: re,
            beta,
            apparent_viscosity: mu,
        })
    }

    fn __str__(&self) -> String {
        format!(
            "VenturiSolver1D(D_in={:.0} μm, D_t={:.0} μm, β={:.3})",
            self.inlet_diameter * 1e6,
            self.throat_diameter * 1e6,
            self.beta()
        )
    }
}

// ── 1D Venturi Result ────────────────────────────────────────────────────────

/// Result from 1D Venturi resistance calculation
#[pyclass(name = "VenturiResult1D")]
#[derive(Debug, Clone, Copy)]
pub struct PyVenturiResult1D {
    #[pyo3(get)]
    pub pressure_drop: f64,
    #[pyo3(get)]
    pub resistance: f64,
    #[pyo3(get)]
    pub dp_bernoulli: f64,
    #[pyo3(get)]
    pub reynolds_number: f64,
    #[pyo3(get)]
    pub beta: f64,
    #[pyo3(get)]
    pub apparent_viscosity: f64,
}

#[pymethods]
impl PyVenturiResult1D {
    fn __str__(&self) -> String {
        format!(
            "VenturiResult1D(ΔP={:.2} Pa, ΔP_bernoulli={:.2} Pa, Re={:.1})",
            self.pressure_drop, self.dp_bernoulli, self.reynolds_number
        )
    }
}
