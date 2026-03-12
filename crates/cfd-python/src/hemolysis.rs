//! Hemolysis model wrappers for `PyO3`

use cfd_core::physics::api::HemolysisModel;
use pyo3::prelude::*;

/// Giersiepen hemolysis power law model
#[pyclass(name = "GiersiepenModel")]
pub struct PyGiersiepenModel {
    inner: HemolysisModel,
}

#[pymethods]
impl PyGiersiepenModel {
    /// Create standard Giersiepen model
    #[new]
    fn new() -> Self {
        PyGiersiepenModel {
            inner: HemolysisModel::giersiepen_standard(),
        }
    }

    /// Calculate hemolysis damage index from shear stress and exposure time
    ///
    /// # Arguments
    /// - `shear_stress`: Shear stress [Pa]
    /// - `exposure_time`: Exposure time [s]
    ///
    /// # Returns
    /// - Blood damage index [-]
    fn calculate_damage(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e: cfd_core::error::Error| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }

    fn __str__(&self) -> String {
        "GiersiepenModel()".to_string()
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}
