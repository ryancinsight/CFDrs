//! Hemolysis model wrappers for `PyO3`

use cfd_core::physics::api::HemolysisModel as RustHemolysisModel;
use pyo3::exceptions::PyValueError;
use pyo3::prelude::*;

/// Hemolysis damage model
///
/// Wraps the underlying `cfd_core` HemolysisModel (Giersiepen power law by default).
#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: RustHemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    /// Create standard Giersiepen hemolysis model
    #[new]
    fn new() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_standard(),
        }
    }

    /// Calculate blood damage index from shear stress and exposure time
    ///
    /// # Arguments
    /// - `shear_stress`: Shear stress magnitude [Pa]
    /// - `exposure_time`: Duration of exposure [s]
    ///
    /// # Returns
    /// - Fractional damage index (ΔHb / Hb)
    fn calculate_damage(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e: cfd_core::error::Error| PyValueError::new_err(e.to_string()))
    }

    fn __str__(&self) -> String {
        "HemolysisModel(Giersiepen)".to_string()
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}
