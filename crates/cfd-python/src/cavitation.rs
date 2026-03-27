//! Cavitation model wrappers for `PyO3`

use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset as RustRayleighPlesset;
use cfd_core::physics::cavitation::regimes::CavitationRegimeClassifier as RustCavitationRegimeClassifier;
use pyo3::prelude::*;

/// Rayleigh-Plesset bubble model
#[pyclass(name = "RayleighPlesset")]
#[derive(Clone)]
pub struct PyRayleighPlesset {
    inner: RustRayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    /// Create new Rayleigh-Plesset model for water
    #[new]
    fn new() -> Self {
        PyRayleighPlesset {
            inner: RustRayleighPlesset::<f64> {
                initial_radius: 10e-6,
                liquid_density: 997.0,
                liquid_viscosity: 0.001,
                surface_tension: 0.0728,
                vapor_pressure: 2339.0,
                polytropic_index: 1.4,
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
    inner: RustCavitationRegimeClassifier<f64>,
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
            inner: RustCavitationRegimeClassifier::new(
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
