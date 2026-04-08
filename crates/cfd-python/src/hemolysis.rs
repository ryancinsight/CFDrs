//! Hemolysis physics model wrappers for `PyO3`

use cfd_core::physics::api::HemolysisModel as RustHemolysisModel;
use pyo3::prelude::*;

/// Hemolysis model wrapper for blood damage estimation
#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: RustHemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    /// Create standard Giersiepen power-law model
    #[staticmethod]
    fn giersiepen_standard() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_standard(),
        }
    }

    /// Calculate blood damage index from shear stress and exposure time
    fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e: cfd_core::error::Error| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
