use cfd_core::physics::api::HemolysisModel as RustHemolysisModel;
use pyo3::prelude::*;
use pyo3::exceptions::PyValueError;

/// Hemolysis model wrapper for PyO3
#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: RustHemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    /// Create Giersiepen model with standard constants
    #[staticmethod]
    fn giersiepen_standard() -> Self {
        Self {
            inner: RustHemolysisModel::giersiepen_standard(),
        }
    }

    /// Calculate blood damage index from shear stress and exposure time
    fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e| PyValueError::new_err(format!("Hemolysis calculation failed: {}", e)))
    }
}
