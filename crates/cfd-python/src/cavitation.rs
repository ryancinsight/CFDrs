//! Cavitation model wrappers for `PyO3`

use cfd_core::physics::cavitation::RayleighPlesset;
use pyo3::prelude::*;

/// Calculate critical Blake radius for unstable growth
#[pyfunction]
#[pyo3(name = "blake_critical_radius")]
pub fn py_blake_critical_radius(p_inf: f64, p_v: f64, sigma: f64) -> f64 {
    // R_c = 0.85 * (2σ/(p_∞ - p_v))
    let model = RayleighPlesset {
        initial_radius: 10e-6, // Dummy value
        liquid_density: 997.0, // Dummy value
        liquid_viscosity: 0.001, // Dummy value
        surface_tension: sigma,
        vapor_pressure: p_v,
        polytropic_index: 1.4, // Dummy value
    };
    model.blake_critical_radius(p_inf)
}

/// Calculate Blake threshold pressure
#[pyfunction]
#[pyo3(name = "blake_threshold")]
pub fn py_blake_threshold(p_inf: f64, p_v: f64, sigma: f64) -> f64 {
    let r_critical = py_blake_critical_radius(p_inf, p_v, sigma);
    p_v + (4.0 / 3.0) * sigma / r_critical
}
