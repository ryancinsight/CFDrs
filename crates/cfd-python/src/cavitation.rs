//! Cavitation physics bindings for PyO3

use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset as RustRayleighPlesset;
use cfd_core::physics::cavitation::regimes::{
    CavitationRegime as RustCavitationRegime,
    CavitationRegimeAnalysis as RustCavitationRegimeAnalysis,
    CavitationRegimeClassifier as RustCavitationRegimeClassifier,
};
use pyo3::prelude::*;

/// Rayleigh-Plesset bubble dynamics model
#[pyclass(name = "RayleighPlesset")]
#[derive(Clone)]
pub struct PyRayleighPlesset {
    pub inner: RustRayleighPlesset<f64>,
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

    /// Calculate Blake critical radius [m]
    fn blake_critical_radius(&self, ambient_pressure: f64) -> f64 {
        self.inner.blake_critical_radius(ambient_pressure)
    }
}

/// Cavitation regime types
#[pyclass(name = "CavitationRegime", eq)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PyCavitationRegime {
    None,
    Stable,
    Inertial,
    Mixed,
}

impl From<RustCavitationRegime> for PyCavitationRegime {
    fn from(regime: RustCavitationRegime) -> Self {
        match regime {
            RustCavitationRegime::None => PyCavitationRegime::None,
            RustCavitationRegime::Stable => PyCavitationRegime::Stable,
            RustCavitationRegime::Inertial => PyCavitationRegime::Inertial,
            RustCavitationRegime::Mixed => PyCavitationRegime::Mixed,
        }
    }
}

/// Cavitation regime classifier
#[pyclass(name = "CavitationRegimeClassifier")]
pub struct PyCavitationRegimeClassifier {
    inner: RustCavitationRegimeClassifier<f64>,
}

#[pymethods]
impl PyCavitationRegimeClassifier {
    #[new]
    fn new(
        bubble_model: &PyRayleighPlesset,
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

    /// Classify cavitation regime
    fn classify_regime(&self) -> PyCavitationRegime {
        self.inner.classify_regime().into()
    }

    /// Calculate Blake threshold pressure [Pa]
    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }

    /// Calculate inertial threshold pressure amplitude [Pa]
    fn inertial_threshold(&self) -> f64 {
        self.inner.inertial_threshold()
    }

    /// Analyze cavitation regime and properties
    fn analyze(&self) -> PyResult<PyCavitationRegimeAnalysis> {
        let analysis = self.inner.analyze().map_err(|e| {
            pyo3::exceptions::PyValueError::new_err(format!("Analysis failed: {}", e))
        })?;

        Ok(PyCavitationRegimeAnalysis { inner: analysis })
    }
}

/// Detailed cavitation regime analysis
#[pyclass(name = "CavitationRegimeAnalysis")]
pub struct PyCavitationRegimeAnalysis {
    inner: RustCavitationRegimeAnalysis<f64>,
}

#[pymethods]
impl PyCavitationRegimeAnalysis {
    #[getter]
    fn regime(&self) -> PyCavitationRegime {
        self.inner.regime.into()
    }

    #[getter]
    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold
    }

    #[getter]
    fn inertial_threshold(&self) -> f64 {
        self.inner.inertial_threshold
    }

    #[getter]
    fn cavitation_number(&self) -> f64 {
        self.inner.cavitation_number
    }

    #[getter]
    fn mechanical_index(&self) -> f64 {
        self.inner.mechanical_index
    }

    #[getter]
    fn max_bubble_radius(&self) -> f64 {
        self.inner.max_bubble_radius
    }

    #[getter]
    fn sonoluminescence_probability(&self) -> f64 {
        self.inner.sonoluminescence_probability
    }

    #[getter]
    fn damage_potential(&self) -> f64 {
        self.inner.damage_potential
    }

    #[getter]
    fn hemolysis_risk(&self) -> f64 {
        self.inner.hemolysis_risk
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }
}
