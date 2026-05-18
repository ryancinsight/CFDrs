//! Cavitation models exposed to Python
use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset;
use pyo3::prelude::*;

#[pyclass(name = "RayleighPlesset")]
pub struct PyRayleighPlesset {
    inner: RayleighPlesset<f64>,
}

impl PyRayleighPlesset {
    pub fn inner(&self) -> RayleighPlesset<f64> {
        self.inner.clone()
    }
}

#[pymethods]
impl PyRayleighPlesset {
    #[new]
    fn new(
        initial_radius: f64,
        liquid_density: f64,
        liquid_viscosity: f64,
        surface_tension: f64,
        vapor_pressure: f64,
        polytropic_index: f64,
    ) -> Self {
        PyRayleighPlesset {
            inner: RayleighPlesset {
                initial_radius,
                liquid_density,
                liquid_viscosity,
                surface_tension,
                vapor_pressure,
                polytropic_index,
            },
        }
    }

    fn blake_critical_radius(&self, ambient_pressure: f64) -> f64 {
        self.inner.blake_critical_radius(ambient_pressure)
    }
}
