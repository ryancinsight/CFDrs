//! Blood rheology model wrappers for PyO3

use cfd_core::physics::fluid::blood::{
    CarreauYasudaBlood as RustCarreauYasudaBlood, CassonBlood as RustCassonBlood,
};
use cfd_core::physics::fluid::traits::{Fluid, NonNewtonianFluid};
use num_traits::FromPrimitive;
use pyo3::prelude::*;

/// Casson blood model with yield stress
///
/// The Casson model represents blood as a non-Newtonian fluid with yield stress:
///
/// ```text
/// μ(γ̇) = (√τ_y + √(μ_∞ · γ̇))²
/// ```
///
/// Where:
/// - τ_y: Yield stress (~0.5 Pa for normal blood)
/// - μ_∞: Viscosity at high shear (~0.003 Pa·s)
/// - γ̇: Shear rate [s⁻¹]
///
/// # Parameters
/// - Normal blood at 37°C
/// - Yield stress: 0.5 Pa (typical range 0.4-0.6 Pa)
/// - Viscosity: 3-10 cP depending on shear rate
#[pyclass(name = "CassonBlood")]
pub struct PyCassonBlood {
    inner: RustCassonBlood<f64>,
}

#[pymethods]
impl PyCassonBlood {
    /// Create Casson blood model with normal parameters
    #[new]
    fn new() -> Self {
        PyCassonBlood {
            inner: RustCassonBlood::<f64>::normal_blood(),
        }
    }

    /// Compute apparent viscosity for given shear rate
    ///
    /// # Arguments
    /// - `gamma`: Wall shear rate [s⁻¹]
    ///
    /// # Returns
    /// - Apparent viscosity [Pa·s]
    fn viscosity(&self, gamma: f64) -> f64 {
        let gamma_val = <f64 as FromPrimitive>::from_f64(gamma).unwrap();
        self.inner.apparent_viscosity(gamma_val)
    }

    /// Get yield stress
    fn yield_stress(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.yield_stress.to_f64().unwrap_or(0.0)
    }

    /// Get density
    fn density(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.density.to_f64().unwrap_or(1060.0)
    }

    /// Get viscosity at high shear (asymptotic value)
    fn viscosity_high_shear(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner
            .infinite_shear_viscosity
            .to_f64()
            .unwrap_or(0.003)
    }

    fn __str__(&self) -> String {
        format!(
            "CassonBlood(τ_y={:.4} Pa, μ_∞={:.6} Pa·s, ρ={:.0} kg/m³)",
            self.yield_stress(),
            self.viscosity_high_shear(),
            self.density()
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

/// Carreau-Yasuda blood model with temperature dependence
///
/// More general non-Newtonian model valid over wider shear rate range:
///
/// ```text
/// μ(γ̇) = μ_∞ + (μ_0 - μ_∞)[1 + (λ·γ̇)ᵃ]^((n-1)/a)
/// ```
///
/// Where:
/// - μ_0: Viscosity at zero shear
/// - μ_∞: Viscosity at infinite shear
/// - λ: Relaxation time
/// - n: Power law index
/// - a: Carreau-Yasuda parameter
#[pyclass(name = "CarreauYasudaBlood")]
pub struct PyCarreauYasudaBlood {
    inner: RustCarreauYasudaBlood<f64>,
}

#[pymethods]
impl PyCarreauYasudaBlood {
    /// Create Carreau-Yasuda blood model with normal parameters
    #[new]
    fn new() -> Self {
        PyCarreauYasudaBlood {
            inner: RustCarreauYasudaBlood::<f64>::normal_blood(),
        }
    }

    /// Compute apparent viscosity for given shear rate
    ///
    /// # Arguments
    /// - `gamma`: Wall shear rate [s⁻¹]
    ///
    /// # Returns
    /// - Apparent viscosity [Pa·s]
    fn viscosity(&self, gamma: f64) -> f64 {
        let gamma_val = <f64 as FromPrimitive>::from_f64(gamma).unwrap();
        self.inner.apparent_viscosity(gamma_val)
    }

    /// Get density
    fn density(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.density.to_f64().unwrap_or(1060.0)
    }

    /// Get viscosity at zero shear (limit)
    fn viscosity_zero_shear(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.zero_shear_viscosity.to_f64().unwrap_or(0.08)
    }

    /// Get viscosity at infinite shear (limit)
    fn viscosity_high_shear(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner
            .infinite_shear_viscosity
            .to_f64()
            .unwrap_or(0.003)
    }

    fn __str__(&self) -> String {
        format!(
            "CarreauYasudaBlood(μ_0={:.4} Pa·s, μ_∞={:.6} Pa·s, ρ={:.0} kg/m³)",
            self.viscosity_zero_shear(),
            self.viscosity_high_shear(),
            self.density()
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}
