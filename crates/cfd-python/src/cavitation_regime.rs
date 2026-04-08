//! Cavitation regime classifier exposed to Python
use cfd_core::physics::cavitation::CavitationRegimeClassifier;
use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset;
use pyo3::prelude::*;
use crate::cavitation::PyRayleighPlesset;

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
        let borrowed = bubble_model.borrow();
        let rust_bubble_model = borrowed.inner();
        PyCavitationRegimeClassifier {
            inner: CavitationRegimeClassifier::new(
                rust_bubble_model,
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
