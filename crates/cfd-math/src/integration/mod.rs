//! Numerical integration methods for CFD applications.
use cfd_core::error::{Error, NumericalErrorKind};
//!
//! This module provides various quadrature rules and integration schemes
//! optimized for CFD simulations with support for adaptive integration.

pub mod composite;
pub mod quadrature;
pub mod tensor;
pub mod traits;
pub mod utils;
pub mod variable;
// Re-export main types for convenience
pub use composite::CompositeQuadrature;
pub use quadrature::{GaussQuadrature, SimpsonsRule, TrapezoidalRule};
pub use tensor::TensorProductQuadrature;
pub use traits::Quadrature;
pub use utils::IntegrationUtils;
pub use variable::VariableQuadrature;
// Type alias for backward compatibility while eliminating adjective-based naming
/// Type alias for variable quadrature (backwards compatibility)
pub type AdaptiveQuadrature<Q> = VariableQuadrature<Q>;
