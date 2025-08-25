//! Numerical integration methods for CFD applications.
//!
//! This module provides various quadrature rules and integration schemes
//! optimized for CFD simulations with support for adaptive integration.

pub mod traits;
pub mod quadrature;
pub mod composite;
pub mod variable;
pub mod tensor;
pub mod utils;

// Re-export main types for convenience
pub use traits::Quadrature;
pub use quadrature::{TrapezoidalRule, SimpsonsRule, GaussQuadrature};
pub use composite::CompositeQuadrature;
pub use variable::VariableQuadrature;
pub use tensor::TensorProductQuadrature;
pub use utils::IntegrationUtils;

// Type alias for backward compatibility while eliminating adjective-based naming
/// Type alias for variable quadrature (backwards compatibility)
pub type AdaptiveQuadrature<Q> = VariableQuadrature<Q>;