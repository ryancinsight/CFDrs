use cfd_core::physics::hemolysis::HemolysisModel;
use pyo3::prelude::*;

/// Hemolysis model types
#[pyclass(name = "HemolysisModel")]
#[derive(Clone)]
pub struct PyHemolysisModel {
    pub(crate) inner: HemolysisModel,
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

    /// Calculate blood damage index from shear stress and exposure time
    fn calculate_damage(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
