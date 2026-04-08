//! Hemolysis models wrapper for `PyO3`

use cfd_core::physics::api::HemolysisModel;
use pyo3::prelude::*;

/// Hemolysis model (Giersiepen)
#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: HemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    /// Create Giersiepen model with standard constants
    #[staticmethod]
    fn giersiepen_standard() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::giersiepen_standard(),
        }
    }

    /// Create Giersiepen model for turbulent flow
    #[staticmethod]
    fn giersiepen_turbulent() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::giersiepen_turbulent(),
        }
    }

    /// Create Giersiepen model for laminar flow
    #[staticmethod]
    fn giersiepen_laminar() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::giersiepen_laminar(),
        }
    }

    /// Create Zhang model for Couette flow
    #[staticmethod]
    fn zhang() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::zhang(),
        }
    }

    /// Create Heuser-Opitz threshold model
    #[staticmethod]
    fn heuser_opitz() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::heuser_opitz(),
        }
    }

    /// Giersiepen (1990) model validated for millifluidic and blood-processing devices
    #[staticmethod]
    fn giersiepen_millifluidic() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::giersiepen_millifluidic(),
        }
    }

    /// Amplify a baseline Haemolysis Index by local cavitation potential
    #[staticmethod]
    fn cavitation_amplified(base_hi: f64, cav_potential: f64) -> f64 {
        HemolysisModel::cavitation_amplified(base_hi, cav_potential)
    }

    /// Calculate blood damage index from shear stress and exposure time
    fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e: cfd_core::error::Error| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
