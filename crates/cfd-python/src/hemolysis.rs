use pyo3::prelude::*;
use cfd_core::physics::api::HemolysisModel as RustHemolysisModel;

/// Giersiepen Hemolysis Model
///
/// Estimates shear stress-related blood damage using the empirical power law:
/// `HI = C * tau^alpha * t^beta`
#[pyclass(name = "GiersiepenModel")]
pub struct PyGiersiepenModel {
    inner: RustHemolysisModel,
}

#[pymethods]
impl PyGiersiepenModel {
    /// Create standard Giersiepen model
    #[new]
    fn new() -> Self {
        PyGiersiepenModel {
            inner: RustHemolysisModel::giersiepen_standard(),
        }
    }

    /// Calculate hemolysis damage index
    fn calculate_damage(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner.damage_index(shear_stress, exposure_time)
            .map_err(|e: cfd_core::error::Error| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
