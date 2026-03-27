//! Hemolysis model wrappers for `PyO3`

use cfd_core::physics::api::HemolysisModel as RustHemolysisModel;
use pyo3::prelude::*;

/// Hemolysis damage model (Giersiepen-Wurzinger)
///
/// Power law model for estimating shear-induced red blood cell damage:
///
/// ```text
/// HI = C · τ^α · t^β
/// ```
#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: RustHemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    /// Create Giersiepen model with standard constants
    #[new]
    fn new() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_standard(),
        }
    }

    /// Calculate blood damage index from shear stress and exposure time
    ///
    /// # Arguments
    /// * `shear_stress` - Shear stress in Pascals
    /// * `exposure_time` - Exposure time in seconds
    fn calculate_damage(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
