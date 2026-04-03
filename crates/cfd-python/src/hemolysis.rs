//! Hemolysis model wrappers for `PyO3`

use cfd_core::physics::api::HemolysisModel;
use pyo3::prelude::*;

/// Hemolysis model for predicting blood damage
#[pyclass(name = "HemolysisModel")]
#[derive(Clone)]
pub struct PyHemolysisModel {
    inner: HemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    /// Create the standard Giersiepen hemolysis model
    #[staticmethod]
    fn giersiepen_standard() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::giersiepen_standard(),
        }
    }

    /// Create the Heuser-Opitz hemolysis model
    #[staticmethod]
    fn heuser_opitz() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::heuser_opitz(),
        }
    }

    /// Calculate the damage index for a given shear stress and exposure time
    ///
    /// # Arguments
    /// - `shear_stress`: Shear stress [Pa]
    /// - `exposure_time`: Exposure time [s]
    ///
    /// # Returns
    /// - Damage index (Delta Hb / Hb)
    fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> f64 {
        self.inner.damage_index(shear_stress, exposure_time).expect("HemolysisModel parameters were invalid")
    }
}
