//! Cavitation models and regime classification wrappers for `PyO3`

use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset;
use cfd_core::physics::cavitation::regimes::CavitationRegimeClassifier;
use pyo3::prelude::*;

/// Rayleigh-Plesset bubble model
#[pyclass(name = "RayleighPlesset")]
#[derive(Clone)]
pub struct PyRayleighPlesset {
    pub(crate) inner: RayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    /// Create new Rayleigh-Plesset bubble model
    #[new]
    #[pyo3(signature = (initial_radius=1e-6, liquid_density=997.0, liquid_viscosity=0.001, surface_tension=0.0728, vapor_pressure=2339.0, polytropic_index=1.4))]
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

    /// Calculate bubble growth rate for inviscid case
    fn growth_rate_inviscid(&self, radius: f64, ambient_pressure: f64) -> f64 {
        self.inner.growth_rate_inviscid(radius, ambient_pressure)
    }

    /// Calculate collapse time from Rayleigh collapse formula
    fn collapse_time(&self, initial_radius: f64, pressure_difference: f64) -> f64 {
        self.inner.collapse_time(initial_radius, pressure_difference)
    }

    /// Calculate maximum bubble radius during growth (Rayleigh-Plesset)
    fn maximum_radius(&self, pressure_ratio: f64) -> f64 {
        self.inner.maximum_radius(pressure_ratio)
    }

    /// Calculate bubble natural frequency
    fn natural_frequency(&self, radius: f64, ambient_pressure: f64) -> f64 {
        self.inner.natural_frequency(radius, ambient_pressure)
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

    /// Calculate inertial cavitation threshold (Apfel & Holland 1991)
    fn inertial_threshold(&self) -> f64 {
        self.inner.inertial_threshold()
    }

    /// Estimate damage potential (0-1 scale)
    fn damage_potential(&self) -> f64 {
        self.inner.damage_potential()
    }

    /// Estimate hemolysis risk for blood flow
    fn hemolysis_risk(&self) -> f64 {
        self.inner.hemolysis_risk()
    }

    /// Estimate sonoluminescence probability based on regime
    fn sonoluminescence_probability(&self) -> f64 {
        self.inner.sonoluminescence_probability()
    }
}
