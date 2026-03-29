use cfd_core::physics::api::HemolysisModel as RustHemolysisModel;
use pyo3::prelude::*;

/// Giersiepen Hemolysis Model
#[pyclass(name = "GiersiepenModel")]
pub struct PyGiersiepenModel {
    inner: RustHemolysisModel,
}

#[pymethods]
impl PyGiersiepenModel {
    #[new]
    fn new() -> Self {
        PyGiersiepenModel {
            inner: RustHemolysisModel::giersiepen_standard(),
        }
    }

    /// Calculate blood damage index from shear stress and exposure time
    fn calculate_damage(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e: cfd_core::error::Error| {
                pyo3::exceptions::PyValueError::new_err(e.to_string())
            })
    }
}
