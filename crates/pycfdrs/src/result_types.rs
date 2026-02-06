//! Result data structures for PyO3 bindings
//!
//! These structures wrap Rust solver results into Python-friendly classes.

use pyo3::prelude::*;

/// Bifurcation solver result
///
/// Contains all outputs from a bifurcation flow calculation including
/// flow rates, pressures, shear rates, and validation errors.
#[pyclass]
#[derive(Debug, Clone)]
pub struct PyBifurcationResult {
    /// Parent inlet flow rate [m³/s]
    #[pyo3(get)]
    pub q_parent: f64,
    /// Daughter 1 flow rate [m³/s]
    #[pyo3(get)]
    pub q_1: f64,
    /// Daughter 2 flow rate [m³/s]
    #[pyo3(get)]
    pub q_2: f64,

    /// Parent inlet pressure [Pa]
    #[pyo3(get)]
    pub p_parent: f64,
    /// Daughter 1 outlet pressure [Pa]
    #[pyo3(get)]
    pub p_1: f64,
    /// Daughter 2 outlet pressure [Pa]
    #[pyo3(get)]
    pub p_2: f64,

    /// Pressure drop in daughter 1 [Pa]
    #[pyo3(get)]
    pub dp_1: f64,
    /// Pressure drop in daughter 2 [Pa]
    #[pyo3(get)]
    pub dp_2: f64,

    /// Wall shear rate in daughter 1 [s⁻¹]
    #[pyo3(get)]
    pub gamma_1: f64,
    /// Wall shear rate in daughter 2 [s⁻¹]
    #[pyo3(get)]
    pub gamma_2: f64,

    /// Apparent viscosity in daughter 1 [Pa·s]
    #[pyo3(get)]
    pub mu_1: f64,
    /// Apparent viscosity in daughter 2 [Pa·s]
    #[pyo3(get)]
    pub mu_2: f64,

    /// Wall shear stress in daughter 1 [Pa]
    #[pyo3(get)]
    pub wss_1: f64,
    /// Wall shear stress in daughter 2 [Pa]
    #[pyo3(get)]
    pub wss_2: f64,

    /// Mass conservation error |Q_1 + Q_2 - Q_parent| / Q_parent
    #[pyo3(get)]
    pub mass_conservation_error: f64,
    /// Pressure continuity error |P_1 - P_2| / P_parent
    #[pyo3(get)]
    pub pressure_continuity_error: f64,
}

#[pymethods]
impl PyBifurcationResult {
    /// Create new bifurcation result
    #[new]
    pub fn new(
        q_parent: f64,
        q_1: f64,
        q_2: f64,
        p_parent: f64,
        p_1: f64,
        p_2: f64,
        dp_1: f64,
        dp_2: f64,
        gamma_1: f64,
        gamma_2: f64,
        mu_1: f64,
        mu_2: f64,
        wss_1: f64,
        wss_2: f64,
        mass_conservation_error: f64,
        pressure_continuity_error: f64,
    ) -> Self {
        PyBifurcationResult {
            q_parent,
            q_1,
            q_2,
            p_parent,
            p_1,
            p_2,
            dp_1,
            dp_2,
            gamma_1,
            gamma_2,
            mu_1,
            mu_2,
            wss_1,
            wss_2,
            mass_conservation_error,
            pressure_continuity_error,
        }
    }

    /// Get flow split ratio Q_1 / Q_parent
    fn flow_split_ratio(&self) -> f64 {
        self.q_1 / (self.q_parent + 1e-15)
    }

    /// Check if solution is valid (conservation errors < 1e-6)
    fn is_valid(&self, tolerance: f64) -> bool {
        self.mass_conservation_error < tolerance && self.pressure_continuity_error < tolerance
    }

    /// Format result as string
    fn __str__(&self) -> String {
        format!(
            "BifurcationResult {{\n  Q: parent={:.3e}, d1={:.3e}, d2={:.3e}\n  \
             P: parent={:.1}, d1={:.1}, d2={:.1} Pa\n  \
             Errors: mass={:.2e}, pressure={:.2e}\n}}",
            self.q_parent,
            self.q_1,
            self.q_2,
            self.p_parent,
            self.p_1,
            self.p_2,
            self.mass_conservation_error,
            self.pressure_continuity_error
        )
    }

    /// Format result nicely for display
    fn __repr__(&self) -> String {
        self.__str__()
    }
}
