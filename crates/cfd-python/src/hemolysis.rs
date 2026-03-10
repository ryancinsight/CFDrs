use cfd_core::physics::api::HemolysisModel;
use pyo3::prelude::*;

#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: HemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    #[staticmethod]
    fn giersiepen_standard() -> Self {
        Self {
            inner: HemolysisModel::giersiepen_standard(),
        }
    }

    #[staticmethod]
    fn giersiepen_millifluidic() -> Self {
        Self {
            inner: HemolysisModel::giersiepen_millifluidic(),
        }
    }

    #[staticmethod]
    fn giersiepen_turbulent() -> Self {
        Self {
            inner: HemolysisModel::giersiepen_turbulent(),
        }
    }

    #[staticmethod]
    fn giersiepen_laminar() -> Self {
        Self {
            inner: HemolysisModel::giersiepen_laminar(),
        }
    }

    #[staticmethod]
    fn zhang() -> Self {
        Self {
            inner: HemolysisModel::zhang(),
        }
    }

    #[staticmethod]
    fn heuser_opitz() -> Self {
        Self {
            inner: HemolysisModel::heuser_opitz(),
        }
    }

    fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
