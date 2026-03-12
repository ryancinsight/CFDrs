//! Cavitation regime models wrapper for `PyO3`

use cfd_core::physics::cavitation::RayleighPlesset;
use cfd_core::physics::cavitation::CavitationRegimeClassifier;
use pyo3::prelude::*;

/// Cavitation Regime Classifier for predicting cavitation behavior
#[pyclass(name = "CavitationRegimeClassifier")]
pub struct PyCavitationRegimeClassifier {
    inner: CavitationRegimeClassifier<f64>,
}

#[pymethods]
impl PyCavitationRegimeClassifier {
    /// Create a new Cavitation Regime Classifier
    #[new]
    #[pyo3(signature = (initial_radius, liquid_density, liquid_viscosity, surface_tension, vapor_pressure, polytropic_index, ambient_pressure))]
    fn new(
        initial_radius: f64,
        liquid_density: f64,
        liquid_viscosity: f64,
        surface_tension: f64,
        vapor_pressure: f64,
        polytropic_index: f64,
        ambient_pressure: f64,
    ) -> Self {
        let bubble_model = RayleighPlesset {
            initial_radius,
            liquid_density,
            liquid_viscosity,
            surface_tension,
            vapor_pressure,
            polytropic_index,
        };

        Self {
            inner: CavitationRegimeClassifier::new(
                bubble_model,
                ambient_pressure,
                None, // no acoustic pressure by default
                None, // no acoustic frequency by default
            ),
        }
    }

    /// Calculate Blake threshold pressure [Pa]
    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }

    /// Calculate critical Blake radius for unstable growth [m]
    fn blake_critical_radius(&self) -> f64 {
        self.inner.bubble_model.blake_critical_radius(self.inner.ambient_pressure)
    }

    fn __str__(&self) -> String {
        format!(
            "CavitationRegimeClassifier(P_ambient={:.1} Pa)",
            self.inner.ambient_pressure
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}
