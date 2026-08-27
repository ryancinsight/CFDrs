//! FEM solver configuration

use cfd_core::error::{Error, Result};
use eunomia::RealField;
use serde::{Deserialize, Serialize};

use eunomia::{FloatElement, NumericElement};

use super::constants;

// Use ElementType from cfd-core as the single source of truth
pub use cfd_core::geometry::ElementType;

/// FEM configuration for fluid dynamics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FemConfig<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Base solver configuration
    pub base: cfd_core::compute::solver::SolverConfig<T>,
    /// Use SUPG/PSPG stabilization
    pub use_stabilization: bool,
    /// Stabilization parameter
    pub tau: T,
    /// Time step (for transient problems)
    pub dt: Option<T>,
    /// Reynolds number (for scaling)
    pub reynolds: Option<T>,
    /// Element type to use
    pub element_type: ElementType,
    /// Quadrature order
    pub quadrature_order: usize,
    /// Grad-div stabilization penalty (0 to disable)
    pub grad_div_penalty: T,
    /// Grad-div stabilization gamma parameter (Olshanskii & Reusken 2004).
    ///
    /// When non-zero, the element-level grad-div parameter is computed as
    /// `tau_div = gamma * h_e^2` (see [`super::stabilization::grad_div_parameter`]).
    /// This provides a mesh-size-aware alternative to the direct `grad_div_penalty`.
    ///
    /// Set to 0.0 to disable. Default: 0.0 (disabled for backward compatibility).
    ///
    /// **Note**: This parameter is available for configuration but is not yet
    /// applied in the assembly loop. The existing `grad_div_penalty` field
    /// provides a direct (non-h-scaled) penalty that is already wired into
    /// the element assembly. A future enhancement will compute
    /// `grad_div_parameter(h_e, grad_div_gamma)` per element and add the
    /// corresponding `tau_div * (div(v), div(w))` contribution.
    pub grad_div_gamma: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy> Default for FemConfig<T> {
    fn default() -> Self {
        Self {
            base: cfd_core::compute::solver::SolverConfig::default(),
            use_stabilization: true,
            tau: <T as FloatElement>::from_f64(constants::DEFAULT_STABILIZATION),
            dt: Some(<T as FloatElement>::from_f64(constants::DEFAULT_TIME_STEP)),
            reynolds: Some(<T as FloatElement>::from_f64(constants::DEFAULT_REYNOLDS)),
            element_type: ElementType::Tetrahedron,
            quadrature_order: constants::DEFAULT_QUADRATURE_ORDER,
            grad_div_penalty: <T as NumericElement>::ZERO,
            grad_div_gamma: <T as NumericElement>::ZERO,
        }
    }
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy> FemConfig<T> {
    /// Construct a `FemConfig` with invariant validation.
    ///
    /// # Errors
    /// Returns `Error::InvalidConfiguration` if any of the following hold:
    /// - `tau` is non-finite or negative,
    /// - `dt = Some(v)` where `v` is non-finite or non-positive,
    /// - `reynolds = Some(v)` where `v` is non-finite or non-positive,
    /// - `quadrature_order` is zero,
    /// - `grad_div_penalty` or `grad_div_gamma` is non-finite or negative.
    pub fn try_new(
        base: cfd_core::compute::solver::SolverConfig<T>,
        tau: T,
        dt: Option<T>,
        reynolds: Option<T>,
        element_type: ElementType,
        quadrature_order: usize,
        grad_div_penalty: T,
        grad_div_gamma: T,
    ) -> Result<Self> {
        if !<T as NumericElement>::is_finite(tau) || tau < <T as NumericElement>::ZERO {
            return Err(Error::InvalidConfiguration(format!(
                "FemConfig::try_new: tau (stabilization parameter) must be finite and non-negative, got {tau:?}"
            )));
        }
        if let Some(dt_value) = dt {
            if !<T as NumericElement>::is_finite(dt_value)
                || dt_value <= <T as NumericElement>::ZERO
            {
                return Err(Error::InvalidConfiguration(format!(
                    "FemConfig::try_new: dt must be finite and positive when provided, got {dt_value:?}"
                )));
            }
        }
        if let Some(re) = reynolds {
            if !<T as NumericElement>::is_finite(re) || re <= <T as NumericElement>::ZERO {
                return Err(Error::InvalidConfiguration(format!(
                    "FemConfig::try_new: reynolds must be finite and positive when provided, got {re:?}"
                )));
            }
        }
        if quadrature_order == 0 {
            return Err(Error::InvalidConfiguration(
                "FemConfig::try_new: quadrature_order must be at least 1".to_string(),
            ));
        }
        if !<T as NumericElement>::is_finite(grad_div_penalty)
            || grad_div_penalty < <T as NumericElement>::ZERO
        {
            return Err(Error::InvalidConfiguration(format!(
                "FemConfig::try_new: grad_div_penalty must be finite and non-negative, got {grad_div_penalty:?}"
            )));
        }
        if !<T as NumericElement>::is_finite(grad_div_gamma)
            || grad_div_gamma < <T as NumericElement>::ZERO
        {
            return Err(Error::InvalidConfiguration(format!(
                "FemConfig::try_new: grad_div_gamma must be finite and non-negative, got {grad_div_gamma:?}"
            )));
        }
        Ok(Self {
            base,
            use_stabilization: true,
            tau,
            dt,
            reynolds,
            element_type,
            quadrature_order,
            grad_div_penalty,
            grad_div_gamma,
        })
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_fem_config_grad_div_defaults_zero() {
        let config = FemConfig::<f64>::default();
        assert_eq!(
            config.grad_div_gamma, 0.0,
            "grad_div_gamma must default to 0.0 (disabled)"
        );
        assert_eq!(
            config.grad_div_penalty, 0.0,
            "grad_div_penalty must default to 0.0 (disabled)"
        );
    }

    #[test]
    fn test_fem_config_grad_div_configurable() {
        let mut config = FemConfig::<f64>::default();
        config.grad_div_gamma = 1.0;
        assert_eq!(config.grad_div_gamma, 1.0);

        config.grad_div_gamma = 5.0;
        assert_eq!(config.grad_div_gamma, 5.0);

        // Setting back to zero disables it
        config.grad_div_gamma = 0.0;
        assert_eq!(config.grad_div_gamma, 0.0);
    }

    /// **Positive**: `try_new` accepts valid configurations.
    #[test]
    fn fem_config_try_new_accepts_defaults() {
        let base = cfd_core::compute::solver::SolverConfig::<f64>::default();
        let config = FemConfig::<f64>::try_new(
            base,
            <f64 as NumericElement>::ZERO,
            None,
            None,
            ElementType::Tetrahedron,
            constants::DEFAULT_QUADRATURE_ORDER,
            <f64 as NumericElement>::ZERO,
            <f64 as NumericElement>::ZERO,
        )
        .expect("zero-stabilization, no-dt, no-reynolds, zero grad_div is valid");
        assert_eq!(config.tau, 0.0);
        assert_eq!(config.dt, None);
        assert_eq!(config.reynolds, None);
        assert_eq!(config.quadrature_order, constants::DEFAULT_QUADRATURE_ORDER);
    }

    /// **Adversarial**: invalid `tau` is rejected.
    #[test]
    fn fem_config_try_new_rejects_invalid_tau() {
        let base = cfd_core::compute::solver::SolverConfig::<f64>::default();
        match FemConfig::<f64>::try_new(
            base.clone(),
            -1.0,
            None,
            None,
            ElementType::Tetrahedron,
            1,
            0.0,
            0.0,
        ) {
            Err(e) => assert!(e.to_string().contains("tau"), "error must mention tau: {e}"),
            Ok(_) => panic!("negative tau must be rejected"),
        }
        match FemConfig::<f64>::try_new(
            base,
            f64::NAN,
            None,
            None,
            ElementType::Tetrahedron,
            1,
            0.0,
            0.0,
        ) {
            Err(e) => assert!(e.to_string().contains("tau"), "error must mention tau: {e}"),
            Ok(_) => panic!("NaN tau must be rejected"),
        }
    }

    /// **Adversarial**: invalid `dt` (when Some) is rejected.
    #[test]
    fn fem_config_try_new_rejects_invalid_dt() {
        let base = cfd_core::compute::solver::SolverConfig::<f64>::default();
        match FemConfig::<f64>::try_new(
            base.clone(),
            0.0,
            Some(0.0),
            None,
            ElementType::Tetrahedron,
            1,
            0.0,
            0.0,
        ) {
            Err(e) => assert!(e.to_string().contains("dt"), "error must mention dt: {e}"),
            Ok(_) => panic!("zero dt must be rejected"),
        }
        match FemConfig::<f64>::try_new(
            base,
            0.0,
            Some(-1.0),
            None,
            ElementType::Tetrahedron,
            1,
            0.0,
            0.0,
        ) {
            Err(e) => assert!(e.to_string().contains("dt"), "error must mention dt: {e}"),
            Ok(_) => panic!("negative dt must be rejected"),
        }
    }

    /// **Adversarial**: invalid `reynolds` (when Some) is rejected.
    #[test]
    fn fem_config_try_new_rejects_invalid_reynolds() {
        let base = cfd_core::compute::solver::SolverConfig::<f64>::default();
        match FemConfig::<f64>::try_new(
            base,
            0.0,
            None,
            Some(-1.0),
            ElementType::Tetrahedron,
            1,
            0.0,
            0.0,
        ) {
            Err(e) => assert!(
                e.to_string().contains("reynolds"),
                "error must mention reynolds: {e}"
            ),
            Ok(_) => panic!("negative reynolds must be rejected"),
        }
    }

    /// **Adversarial**: zero `quadrature_order` is rejected.
    #[test]
    fn fem_config_try_new_rejects_zero_quadrature_order() {
        let base = cfd_core::compute::solver::SolverConfig::<f64>::default();
        match FemConfig::<f64>::try_new(
            base,
            0.0,
            None,
            None,
            ElementType::Tetrahedron,
            0,
            0.0,
            0.0,
        ) {
            Err(e) => assert!(
                e.to_string().contains("quadrature_order"),
                "error must mention quadrature_order: {e}"
            ),
            Ok(_) => panic!("zero quadrature_order must be rejected"),
        }
    }

    /// **Adversarial**: invalid `grad_div_penalty` or `grad_div_gamma` is rejected.
    #[test]
    fn fem_config_try_new_rejects_invalid_grad_div() {
        let base = cfd_core::compute::solver::SolverConfig::<f64>::default();
        match FemConfig::<f64>::try_new(
            base.clone(),
            0.0,
            None,
            None,
            ElementType::Tetrahedron,
            1,
            -1.0,
            0.0,
        ) {
            Err(e) => assert!(
                e.to_string().contains("grad_div_penalty"),
                "error must mention grad_div_penalty: {e}"
            ),
            Ok(_) => panic!("negative grad_div_penalty must be rejected"),
        }
        match FemConfig::<f64>::try_new(
            base,
            0.0,
            None,
            None,
            ElementType::Tetrahedron,
            1,
            0.0,
            f64::NAN,
        ) {
            Err(e) => assert!(
                e.to_string().contains("grad_div_gamma"),
                "error must mention grad_div_gamma: {e}"
            ),
            Ok(_) => panic!("NaN grad_div_gamma must be rejected"),
        }
    }
}
