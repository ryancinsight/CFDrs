use pyo3::prelude::*;
use cfd_core::physics::cavitation::{
    RayleighPlesset as RustRayleighPlesset,
    CavitationRegimeClassifier as RustCavitationRegimeClassifier
};

/// Rayleigh-Plesset bubble dynamics model
#[pyclass(name = "RayleighPlesset")]
pub struct PyRayleighPlesset {
    inner: RustRayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    /// Create new Rayleigh-Plesset model
    #[new]
    #[pyo3(signature = (initial_radius, liquid_density=997.0, liquid_viscosity=0.001, surface_tension=0.0728, vapor_pressure=2339.0, polytropic_index=1.4))]
    fn new(
        initial_radius: f64,
        liquid_density: f64,
        liquid_viscosity: f64,
        surface_tension: f64,
        vapor_pressure: f64,
        polytropic_index: f64,
    ) -> Self {
        PyRayleighPlesset {
            inner: RustRayleighPlesset {
                initial_radius,
                liquid_density,
                liquid_viscosity,
                surface_tension,
                vapor_pressure,
                polytropic_index,
            },
        }
    }

    /// Calculate critical Blake radius
    fn blake_critical_radius(&self, ambient_pressure: f64) -> f64 {
        self.inner.blake_critical_radius(ambient_pressure)
    }
}

/// Cavitation regime classifier
#[pyclass(name = "CavitationRegimeClassifier")]
pub struct PyCavitationRegimeClassifier {
    inner: RustCavitationRegimeClassifier<f64>,
}

#[pymethods]
impl PyCavitationRegimeClassifier {
    /// Create new regime classifier
    #[new]
    #[pyo3(signature = (bubble_model, ambient_pressure, acoustic_pressure=None, acoustic_frequency=None))]
    fn new(
        bubble_model: &PyRayleighPlesset,
        ambient_pressure: f64,
        acoustic_pressure: Option<f64>,
        acoustic_frequency: Option<f64>,
    ) -> Self {
        PyCavitationRegimeClassifier {
            inner: RustCavitationRegimeClassifier::new(
                bubble_model.inner.clone(),
                ambient_pressure,
                acoustic_pressure,
                acoustic_frequency,
            ),
        }
    }

    /// Calculate Blake threshold
    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }

    /// Calculate inertial threshold
    fn inertial_threshold(&self) -> f64 {
        self.inner.inertial_threshold()
    }
}
