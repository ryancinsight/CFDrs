use cfd_core::physics::cavitation::{CavitationRegimeClassifier, RayleighPlesset};
use pyo3::prelude::*;

/// Rayleigh-Plesset bubble model
#[pyclass(name = "RayleighPlesset")]
#[derive(Clone)]
pub struct PyRayleighPlesset {
    pub(crate) inner: RayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    /// Create a new Rayleigh-Plesset model
    #[new]
    #[pyo3(signature = (initial_radius, liquid_density, liquid_viscosity, surface_tension, vapor_pressure, polytropic_index))]
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

    /// Calculate critical Blake radius for unstable growth
    fn blake_critical_radius(&self, ambient_pressure: f64) -> f64 {
        self.inner.blake_critical_radius(ambient_pressure)
    }
}

/// Cavitation regime classifier
#[pyclass(name = "CavitationRegimeClassifier")]
pub struct PyCavitationRegimeClassifier {
    inner: CavitationRegimeClassifier<f64>,
}

#[pymethods]
impl PyCavitationRegimeClassifier {
    /// Create new cavitation regime classifier
    #[new]
    #[pyo3(signature = (bubble_model, ambient_pressure, acoustic_pressure=None, acoustic_frequency=None))]
    fn new(
        bubble_model: PyRayleighPlesset,
        ambient_pressure: f64,
        acoustic_pressure: Option<f64>,
        acoustic_frequency: Option<f64>,
    ) -> Self {
        PyCavitationRegimeClassifier {
            inner: CavitationRegimeClassifier::new(
                bubble_model.inner,
                ambient_pressure,
                acoustic_pressure,
                acoustic_frequency,
            ),
        }
    }

    /// Calculate Blake threshold pressure
    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }
}
