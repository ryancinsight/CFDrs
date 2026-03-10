use cfd_core::physics::cavitation::{RayleighPlesset, CavitationRegimeClassifier};
use pyo3::prelude::*;

#[pyclass(name = "RayleighPlesset")]
pub struct PyRayleighPlesset {
    pub(crate) inner: RayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    fn blake_critical_radius(&self, ambient_pressure: f64) -> f64 {
        self.inner.blake_critical_radius(ambient_pressure)
    }

    #[new]
    #[pyo3(signature = (initial_radius, liquid_density, liquid_viscosity, surface_tension, vapor_pressure, polytropic_index=1.33))]
    fn new(
        initial_radius: f64,
        liquid_density: f64,
        liquid_viscosity: f64,
        surface_tension: f64,
        vapor_pressure: f64,
        polytropic_index: f64,
    ) -> Self {
        Self {
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
}

#[pyclass(name = "CavitationRegimeClassifier")]
pub struct PyCavitationRegimeClassifier {
    inner: CavitationRegimeClassifier<f64>,
}

#[pymethods]
impl PyCavitationRegimeClassifier {
    #[new]
    #[pyo3(signature = (bubble_model, ambient_pressure, acoustic_pressure=None, acoustic_frequency=None))]
    fn new(
        bubble_model: &Bound<'_, PyRayleighPlesset>,
        ambient_pressure: f64,
        acoustic_pressure: Option<f64>,
        acoustic_frequency: Option<f64>,
    ) -> Self {
        let model = bubble_model.borrow().inner.clone();
        Self {
            inner: CavitationRegimeClassifier::new(
                model,
                ambient_pressure,
                acoustic_pressure,
                acoustic_frequency,
            ),
        }
    }

    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }
}
