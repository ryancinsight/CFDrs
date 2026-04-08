//! Hemolysis models for `PyO3`

use cfd_core::physics::api::HemolysisModel as RustHemolysisModel;
use pyo3::prelude::*;

#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: RustHemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    #[staticmethod]
    fn giersiepen_standard() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_standard(),
        }
    }

    #[staticmethod]
    fn giersiepen_millifluidic() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_millifluidic(),
        }
    }

    fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner.damage_index(shear_stress, exposure_time).map_err(|e: cfd_core::error::Error| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
