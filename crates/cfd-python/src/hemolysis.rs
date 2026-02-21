//! Hemolysis physics bindings for PyO3

use cfd_core::physics::hemolysis::{
    BloodTrauma as RustBloodTrauma, BloodTraumaSeverity as RustBloodTraumaSeverity,
    HemolysisCalculator as RustHemolysisCalculator, HemolysisModel as RustHemolysisModel,
};
use pyo3::prelude::*;

/// Hemolysis model for blood damage prediction
#[pyclass(name = "HemolysisModel")]
#[derive(Clone)]
pub struct PyHemolysisModel {
    pub inner: RustHemolysisModel,
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
    fn giersiepen_turbulent() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_turbulent(),
        }
    }

    #[staticmethod]
    fn giersiepen_laminar() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::giersiepen_laminar(),
        }
    }

    #[staticmethod]
    fn zhang() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::zhang(),
        }
    }

    #[staticmethod]
    fn heuser_opitz() -> Self {
        PyHemolysisModel {
            inner: RustHemolysisModel::heuser_opitz(),
        }
    }

    fn damage_index(&self, shear_stress: f64, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .damage_index(shear_stress, exposure_time)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}

/// Hemolysis calculator
#[pyclass(name = "HemolysisCalculator")]
pub struct PyHemolysisCalculator {
    inner: RustHemolysisCalculator<f64>,
}

#[pymethods]
impl PyHemolysisCalculator {
    #[new]
    fn new(
        model: &PyHemolysisModel,
        hematocrit: f64,
        hemoglobin_initial: f64,
        blood_volume: f64,
        flow_rate: f64,
    ) -> PyResult<Self> {
        let inner = RustHemolysisCalculator::new(
            model.inner.clone(),
            hematocrit,
            hemoglobin_initial,
            blood_volume,
            flow_rate,
        ).map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))?;

        Ok(PyHemolysisCalculator { inner })
    }

    fn normalized_index(&self, delta_hemoglobin: f64) -> f64 {
        self.inner.normalized_index(delta_hemoglobin)
    }

    fn modified_index(&self, delta_hemoglobin: f64, exposure_time: f64) -> f64 {
        self.inner.modified_index(delta_hemoglobin, exposure_time)
    }

    fn hemoglobin_release(&self, damage_index: f64) -> f64 {
        self.inner.hemoglobin_release(damage_index)
    }

    fn critical_shear_stress(&self, exposure_time: f64) -> PyResult<f64> {
        self.inner
            .critical_shear_stress(exposure_time)
            .map_err(|e| pyo3::exceptions::PyValueError::new_err(e.to_string()))
    }
}

/// Blood trauma severity classification
#[pyclass(name = "BloodTraumaSeverity", eq, eq_int)]
#[derive(Clone, Copy, PartialEq, Eq)]
pub enum PyBloodTraumaSeverity {
    Minimal,
    Moderate,
    Severe,
    Critical,
}

impl From<RustBloodTraumaSeverity> for PyBloodTraumaSeverity {
    fn from(severity: RustBloodTraumaSeverity) -> Self {
        match severity {
            RustBloodTraumaSeverity::Minimal => PyBloodTraumaSeverity::Minimal,
            RustBloodTraumaSeverity::Moderate => PyBloodTraumaSeverity::Moderate,
            RustBloodTraumaSeverity::Severe => PyBloodTraumaSeverity::Severe,
            RustBloodTraumaSeverity::Critical => PyBloodTraumaSeverity::Critical,
        }
    }
}

/// Blood trauma assessment
#[pyclass(name = "BloodTrauma")]
pub struct PyBloodTrauma {
    inner: RustBloodTrauma,
}

#[pymethods]
impl PyBloodTrauma {
    #[getter]
    fn hemolysis_level(&self) -> f64 {
        self.inner.hemolysis_level
    }

    #[getter]
    fn platelet_activation(&self) -> f64 {
        self.inner.platelet_activation
    }

    #[getter]
    fn thrombosis_risk(&self) -> f64 {
        self.inner.thrombosis_risk
    }

    #[getter]
    fn max_shear_stress(&self) -> f64 {
        self.inner.max_shear_stress
    }

    #[getter]
    fn avg_exposure_time(&self) -> f64 {
        self.inner.avg_exposure_time
    }

    fn severity(&self) -> PyBloodTraumaSeverity {
        self.inner.severity().into()
    }

    fn meets_fda_guidance(&self) -> bool {
        self.inner.meets_fda_guidance()
    }

    fn __str__(&self) -> String {
        self.inner.to_string()
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}
