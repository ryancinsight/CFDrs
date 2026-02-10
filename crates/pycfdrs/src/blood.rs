//! Blood rheology model wrappers for PyO3

use cfd_core::physics::fluid::blood::{
    CarreauYasudaBlood as RustCarreauYasudaBlood, CassonBlood as RustCassonBlood,
    CrossBlood as RustCrossBlood, FahraeuasLindqvist as RustFahraeuasLindqvist,
};
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
    /// Compute apparent viscosity for given shear rate
    fn apparent_viscosity(&self, gamma: f64) -> f64 {
        self.viscosity(gamma)
    }

    /// Compute apparent viscosity for given shear rate
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
    /// Compute apparent viscosity for given shear rate
    fn apparent_viscosity(&self, gamma: f64) -> f64 {
        self.viscosity(gamma)
    }

    /// Compute apparent viscosity for given shear rate
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

/// Cross blood rheology model
///
/// Simpler alternative to Carreau-Yasuda:
///
/// ```text
/// mu(gamma) = mu_inf + (mu_0 - mu_inf) / (1 + (K * gamma)^n)
/// ```
///
/// Parameters fitted for normal human blood.
#[pyclass(name = "CrossBlood")]
pub struct PyCrossBlood {
    inner: RustCrossBlood<f64>,
}

#[pymethods]
impl PyCrossBlood {
    /// Create Cross blood model with normal parameters
    #[new]
    fn new() -> Self {
        PyCrossBlood {
            inner: RustCrossBlood::<f64>::normal_blood(),
        }
    }

    /// Compute apparent viscosity for given shear rate [Pa.s]
    fn apparent_viscosity(&self, gamma: f64) -> f64 {
        self.inner.apparent_viscosity(gamma)
    }

    /// Alias for apparent_viscosity
    fn viscosity(&self, gamma: f64) -> f64 {
        self.inner.apparent_viscosity(gamma)
    }

    /// Get density [kg/m^3]
    fn density(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.density.to_f64().unwrap_or(1060.0)
    }

    /// Get viscosity at zero shear (limit) [Pa.s]
    fn viscosity_zero_shear(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.zero_shear_viscosity.to_f64().unwrap_or(0.056)
    }

    /// Get viscosity at infinite shear (limit) [Pa.s]
    fn viscosity_high_shear(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner
            .infinite_shear_viscosity
            .to_f64()
            .unwrap_or(0.00345)
    }

    /// Get time constant K [s]
    fn time_constant(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.time_constant.to_f64().unwrap_or(1.007)
    }

    /// Get rate index n [-]
    fn rate_index(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.rate_index.to_f64().unwrap_or(1.028)
    }

    fn __str__(&self) -> String {
        format!(
            "CrossBlood(mu_0={:.4} Pa.s, mu_inf={:.6} Pa.s, K={:.4} s, n={:.4})",
            self.viscosity_zero_shear(),
            self.viscosity_high_shear(),
            self.time_constant(),
            self.rate_index()
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}

/// Fahraeus-Lindqvist effect calculator for microvascular blood flow
///
/// In microvessels (D < 300 um), blood exhibits apparent viscosity reduction
/// due to axial migration of RBCs (cell-free layer near wall).
///
/// Uses simplified Pries et al. (1992) correlation.
#[pyclass(name = "FahraeuasLindqvist")]
pub struct PyFahraeuasLindqvist {
    inner: RustFahraeuasLindqvist<f64>,
}

#[pymethods]
impl PyFahraeuasLindqvist {
    /// Create Fahraeus-Lindqvist calculator
    ///
    /// # Arguments
    /// - `diameter`: Vessel diameter [m]
    /// - `hematocrit`: Volume fraction of RBCs [-] (default 0.45)
    #[new]
    #[pyo3(signature = (diameter, hematocrit=0.45))]
    fn new(diameter: f64, hematocrit: f64) -> Self {
        PyFahraeuasLindqvist {
            inner: RustFahraeuasLindqvist::<f64>::new(diameter, hematocrit),
        }
    }

    /// Check if Fahraeus-Lindqvist effect is significant (D < 300 um)
    fn is_significant(&self) -> bool {
        self.inner.is_significant()
    }

    /// Calculate relative apparent viscosity mu_app / mu_plasma [-]
    fn relative_viscosity(&self) -> f64 {
        self.inner.relative_viscosity()
    }

    /// Calculate apparent viscosity in microvessel [Pa.s]
    fn apparent_viscosity(&self) -> f64 {
        self.inner.apparent_viscosity()
    }

    /// Calculate tube hematocrit (Fahraeus effect) [-]
    fn tube_hematocrit(&self) -> f64 {
        self.inner.tube_hematocrit()
    }

    /// Get plasma viscosity [Pa.s]
    fn plasma_viscosity(&self) -> f64 {
        use num_traits::ToPrimitive;
        self.inner.plasma_viscosity.to_f64().unwrap_or(0.00122)
    }

    fn __str__(&self) -> String {
        use num_traits::ToPrimitive;
        format!(
            "FahraeuasLindqvist(D={:.1} um, Ht={:.2}, significant={})",
            self.inner.diameter.to_f64().unwrap_or(0.0) * 1e6,
            self.inner.hematocrit.to_f64().unwrap_or(0.45),
            self.is_significant()
        )
    }

    fn __repr__(&self) -> String {
        self.__str__()
    }
}
