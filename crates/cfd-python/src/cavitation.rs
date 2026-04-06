use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset;
use cfd_core::physics::cavitation::regimes::CavitationRegimeClassifier;
use pyo3::prelude::*;

/// Rayleigh-Plesset wrapper for PyO3
#[pyclass(name = "RayleighPlesset")]
pub struct PyRayleighPlesset {
    inner: RayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    #[new]
    pub fn new(initial_radius: f64) -> Self {
        Self {
            inner: RayleighPlesset {
                initial_radius,
                liquid_density: 998.0,
                liquid_viscosity: 1.002e-3,
                surface_tension: 0.0728,
                vapor_pressure: 2339.0,
                polytropic_index: 1.4,
            },
        }
    }

    pub fn blake_critical_radius(&self, ambient_pressure: f64) -> f64 {
        self.inner.blake_critical_radius(ambient_pressure)
    }
}

/// Cavitation Regime Classifier wrapper for PyO3
#[pyclass(name = "CavitationRegimeClassifier")]
pub struct PyCavitationRegimeClassifier {
    inner: CavitationRegimeClassifier<f64>,
}

#[pymethods]
impl PyCavitationRegimeClassifier {
    #[new]
    #[pyo3(signature = (bubble_model, ambient_pressure, acoustic_pressure=None, acoustic_frequency=None))]
    pub fn new(
        bubble_model: &PyRayleighPlesset,
        ambient_pressure: f64,
        acoustic_pressure: Option<f64>,
        acoustic_frequency: Option<f64>,
    ) -> Self {
        Self {
            inner: CavitationRegimeClassifier::new(
                bubble_model.inner.clone(),
                ambient_pressure,
                acoustic_pressure,
                acoustic_frequency,
            ),
        }
    }

    pub fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }
}
