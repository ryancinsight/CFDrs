use cfd_core::physics::api::HemolysisModel;
use pyo3::prelude::*;

/// Giersiepen Hemolysis Model (Power Law)
#[pyclass(name = "GiersiepenModel")]
pub struct PyGiersiepenModel {
    model: HemolysisModel,
}

#[pymethods]
impl PyGiersiepenModel {
    /// Create new Giersiepen model
    #[new]
    fn new() -> Self {
        PyGiersiepenModel {
            model: HemolysisModel::giersiepen_standard(),
        }
    }

    /// Calculate damage based on shear stress and exposure time
    fn calculate_damage(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.model.damage_index(shear_stress, exposure_time)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
