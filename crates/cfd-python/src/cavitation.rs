use cfd_core::physics::api::{CavitationRegime, CavitationRegimeClassifier, RayleighPlesset};
use pyo3::prelude::*;

#[pyclass(name = "CavitationRegime", eq, eq_int)]
#[derive(Clone, PartialEq, Eq)]
pub enum PyCavitationRegime {
    None,
    Stable,
    Inertial,
    Mixed,
}

impl From<CavitationRegime> for PyCavitationRegime {
    fn from(r: CavitationRegime) -> Self {
        match r {
            CavitationRegime::None => PyCavitationRegime::None,
            CavitationRegime::Stable => PyCavitationRegime::Stable,
            CavitationRegime::Inertial => PyCavitationRegime::Inertial,
            CavitationRegime::Mixed => PyCavitationRegime::Mixed,
        }
    }
}

#[pyclass(name = "RayleighPlesset")]
#[derive(Clone)]
pub struct PyRayleighPlesset {
    inner: RayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    #[new]
    #[pyo3(signature = (initial_radius=10e-6, liquid_density=997.0, liquid_viscosity=0.001, surface_tension=0.0728, vapor_pressure=2339.0, polytropic_index=1.4))]
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

    #[staticmethod]
    fn water_standard() -> Self {
        Self::new(10e-6, 997.0, 0.001, 0.0728, 2339.0, 1.4)
    }

    fn blake_critical_radius(&self, p_inf: f64) -> f64 {
        self.inner.blake_critical_radius(p_inf)
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
        bubble_model: PyRayleighPlesset,
        ambient_pressure: f64,
        acoustic_pressure: Option<f64>,
        acoustic_frequency: Option<f64>,
    ) -> Self {
        Self {
            inner: CavitationRegimeClassifier::new(
                bubble_model.inner,
                ambient_pressure,
                acoustic_pressure,
                acoustic_frequency,
            ),
        }
    }

    #[staticmethod]
    fn water_standard() -> Self {
        let rp = PyRayleighPlesset::water_standard();
        Self::new(rp, 101325.0, None, None)
    }

    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }
}
