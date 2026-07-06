//! FEM solver configuration

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

use eunomia::{FloatElement, NumericElement, RealField as EunomiaRealField};

use super::{constants, scalar};

// Use ElementType from cfd-core as the single source of truth
pub use cfd_core::geometry::ElementType;

/// FEM configuration for fluid dynamics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FemConfig<T: cfd_mesh::domain::core::Scalar + RealField + EunomiaRealField + Copy> {
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

impl<T: cfd_mesh::domain::core::Scalar + RealField + EunomiaRealField + FloatElement + Copy> Default
    for FemConfig<T>
{
    fn default() -> Self {
        Self {
            base: cfd_core::compute::solver::SolverConfig::default(),
            use_stabilization: true,
            tau: scalar::constant(constants::DEFAULT_STABILIZATION),
            dt: Some(scalar::constant(constants::DEFAULT_TIME_STEP)),
            reynolds: Some(scalar::constant(constants::DEFAULT_REYNOLDS)),
            element_type: ElementType::Tetrahedron,
            quadrature_order: constants::DEFAULT_QUADRATURE_ORDER,
            grad_div_penalty: <T as NumericElement>::ZERO,
            grad_div_gamma: <T as NumericElement>::ZERO,
        }
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
}
