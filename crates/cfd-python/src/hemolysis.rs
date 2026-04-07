use cfd_core::physics::hemolysis::models::HemolysisModel;
use pyo3::prelude::*;

#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: HemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    #[staticmethod]
    fn giersiepen_standard() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::giersiepen_standard(),
        }
    }

    #[staticmethod]
    fn giersiepen_millifluidic() -> Self {
        PyHemolysisModel {
            inner: HemolysisModel::giersiepen_millifluidic(),
        }
    }

    fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
