//! Hemolysis model wrappers for `PyO3`

use cfd_core::physics::hemolysis::HemolysisModel as RustHemolysisModel;
use pyo3::prelude::*;

#[pyclass(name = "HemolysisModel")]
pub struct PyHemolysisModel {
    inner: RustHemolysisModel,
}

#[pymethods]
impl PyHemolysisModel {
    /// Create Giersiepen model with standard constants
    #[staticmethod]
    pub fn giersiepen_standard() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_standard(),
        }
    }

    /// Create Giersiepen model for turbulent flow
    #[staticmethod]
    pub fn giersiepen_turbulent() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_turbulent(),
        }
    }

    /// Create Giersiepen model for laminar flow
    #[staticmethod]
    pub fn giersiepen_laminar() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_laminar(),
        }
    }

    /// Create Zhang model for Couette flow
    #[staticmethod]
    pub fn zhang() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::zhang(),
        }
    }

    /// Create Heuser-Opitz threshold model
    #[staticmethod]
    pub fn heuser_opitz() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::heuser_opitz(),
        }
    }

    /// Calculate blood damage index from shear stress and exposure time
    pub fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner.damage_index(shear_stress, exposure_time).map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}
