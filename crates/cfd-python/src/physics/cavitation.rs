use pyo3::prelude::*;
use cfd_core::physics::cavitation::regimes::CavitationRegimeClassifier as RustCavitationRegimeClassifier;
use cfd_core::physics::cavitation::rayleigh_plesset::RayleighPlesset as RustRayleighPlesset;

#[pyclass(name = "CavitationRegimeClassifier")]
pub struct PyCavitationRegimeClassifier {
    inner: RustCavitationRegimeClassifier<f64>,
}

#[pymethods]
impl PyCavitationRegimeClassifier {
    #[new]
    fn new(ambient_pressure: f64, vapor_pressure: f64, surface_tension: f64) -> Self {
        let bubble_model = RustRayleighPlesset {
            initial_radius: 10e-6,
            vapor_pressure,
            surface_tension,
            liquid_density: 1000.0,
            liquid_viscosity: 0.001,
            polytropic_index: 1.4,
        };
        Self {
            inner: RustCavitationRegimeClassifier::new(bubble_model, ambient_pressure, None, None),
        }
    }

    // Removing blake_critical_radius from CavitationRegimeClassifier.
    // It's already exposed on RayleighPlesset, and we'll use that instead.

    fn blake_threshold(&self) -> f64 {
        self.inner.blake_threshold()
    }
}

#[pyclass(name = "RayleighPlesset")]
pub struct PyRayleighPlesset {
    inner: RustRayleighPlesset<f64>,
}

#[pymethods]
impl PyRayleighPlesset {
    #[new]
    fn new(vapor_pressure: f64, surface_tension: f64, fluid_density: f64, fluid_viscosity: f64) -> Self {
        Self {
            inner: RustRayleighPlesset {
                initial_radius: 10e-6, // Some default, it isn't used for blake_critical_radius
                vapor_pressure,
                surface_tension,
                liquid_density: fluid_density,
                liquid_viscosity: fluid_viscosity,
                polytropic_index: 1.4,
            },
        }
    }

    fn blake_critical_radius(&self, ambient_pressure: f64) -> f64 {
        self.inner.blake_critical_radius(ambient_pressure)
    }
}
