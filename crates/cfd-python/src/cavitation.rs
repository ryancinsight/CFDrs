//! Cavitation model wrappers for `PyO3`

use cfd_core::physics::cavitation::regimes::CavitationRegime;
use cfd_core::physics::cavitation::{CavitationRegimeClassifier, RayleighPlesset};
use pyo3::prelude::*;

/// Rayleigh-Plesset bubble dynamics model
///
/// This model predicts the evolution of a spherical bubble in an infinite liquid
/// domain under varying far-field pressure.
#[pyclass(name = "RayleighPlesset")]
#[derive(Clone, Copy)]
pub struct PyRayleighPlesset {
    inner: RayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    /// Create Rayleigh-Plesset model
    ///
    /// # Arguments
    /// - `initial_radius`: Initial bubble radius [m]
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

    /// Calculate the critical radius R_c for the Blake threshold
    ///
    /// The critical radius is the size above which a bubble will grow explosively
    /// when the ambient pressure drops below the vapor pressure.
    ///
    /// Formula: R_c = 0.85 * 2 * sigma / (P_inf - P_v)
    ///
    /// # Arguments
    /// - `p_inf`: Far-field ambient pressure [Pa]
    fn blake_critical_radius(&self, p_inf: f64) -> f64 {
        self.inner.blake_critical_radius(p_inf)
    }

}

/// Classifier for cavitation flow regimes
#[pyclass(name = "CavitationRegimeClassifier")]
pub struct PyCavitationRegimeClassifier {
    inner: CavitationRegimeClassifier<f64>,
}

#[pymethods]
impl PyCavitationRegimeClassifier {
    /// Create classifier
    ///
    /// # Arguments
    /// - `bubble_model`: RayleighPlesset instance
    /// - `ambient_pressure`: Ambient pressure [Pa]
    /// - `acoustic_pressure`: Acoustic pressure amplitude [Pa]
    /// - `acoustic_frequency`: Acoustic frequency [Hz]
    #[new]
    #[pyo3(signature = (bubble_model, ambient_pressure, acoustic_pressure=None, acoustic_frequency=None))]
    fn new(
        bubble_model: &PyRayleighPlesset,
        ambient_pressure: f64,
        acoustic_pressure: Option<f64>,
        acoustic_frequency: Option<f64>,
    ) -> Self {
        PyCavitationRegimeClassifier {
            inner: CavitationRegimeClassifier::<f64>::new(
                bubble_model.inner,
                ambient_pressure,
                acoustic_pressure,
                acoustic_frequency,
            ),
        }
    }

    /// Calculate the Blake threshold pressure
    ///
    /// The Blake threshold is the critical tension at which a bubble of a given
    /// initial size becomes unstable and grows rapidly.
    ///
    /// Formula: P_Blake = P_v + 4/3 * sigma / R_c
    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }

    /// Classify flow regime into stable, inertial, or mixed cavitation
    fn classify_regime(&self) -> String {
        match self.inner.classify_regime() {
            CavitationRegime::None => "None".to_string(),
            CavitationRegime::Stable => "Stable".to_string(),
            CavitationRegime::Inertial => "Inertial".to_string(),
            CavitationRegime::Mixed => "Mixed".to_string(),
        }
    }
}
