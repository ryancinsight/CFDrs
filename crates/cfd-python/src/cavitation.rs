use cfd_core::physics::cavitation::regimes::CavitationRegimeClassifier;
use cfd_core::physics::cavitation::RayleighPlesset;
use pyo3::prelude::*;

/// Calculator for Blake Threshold for cavitation
#[pyclass(name = "BlakeThreshold")]
pub struct PyBlakeThreshold {
    classifier: CavitationRegimeClassifier<f64>,
    ambient_pressure: f64,
}

#[pymethods]
impl PyBlakeThreshold {
    /// Create new Blake Threshold calculator
    #[new]
    fn new(ambient_pressure: f64) -> Self {
        let model = RayleighPlesset {
            initial_radius: 10e-6,
            liquid_density: 997.0,
            liquid_viscosity: 0.001,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };

        PyBlakeThreshold {
            classifier: CavitationRegimeClassifier::new(model, ambient_pressure, None, None),
            ambient_pressure,
        }
    }

    /// Calculate critical radius (Rc)
    fn critical_radius(&self) -> f64 {
        let model = RayleighPlesset {
            initial_radius: 10e-6,
            liquid_density: 997.0,
            liquid_viscosity: 0.001,
            surface_tension: 0.0728,
            vapor_pressure: 2339.0,
            polytropic_index: 1.4,
        };
        model.blake_critical_radius(self.ambient_pressure)
    }

    /// Calculate Blake threshold pressure
    fn threshold(&self) -> f64 {
        self.classifier.blake_threshold()
    }
}
