//! Domain-specific error types for boundary condition operations
//!
//! Provides type-safe error handling with clear failure modes

use thiserror::Error;

/// Boundary condition application errors
#[derive(Debug, Error)]
pub enum BoundaryError {
    /// Insufficient interior values for requested stencil order
    #[error(
        "Insufficient interior values: need {required} for order {order} scheme, got {actual}"
    )]
    InsufficientStencil {
        /// Required number of interior values
        required: usize,
        /// Stencil order
        order: usize,
        /// Actual number provided
        actual: usize,
    },

    /// Unsupported stencil order
    #[error("Stencil order {0} not implemented (supported: 1-4)")]
    UnsupportedOrder(usize),

    /// Robin condition singularity
    #[error("Robin condition denominator near zero: β + α*dx = {value}")]
    RobinSingularity {
        /// Denominator value
        value: f64,
    },

    /// Invalid boundary region
    #[error("Invalid boundary region: {0}")]
    InvalidRegion(String),

    /// Dimension mismatch
    #[error("Dimension mismatch: expected {expected}, got {actual}")]
    DimensionMismatch {
        /// Expected dimensions
        expected: String,
        /// Actual dimensions
        actual: String,
    },
}

impl BoundaryError {
    /// Create insufficient stencil error
    pub fn insufficient_stencil(required: usize, order: usize, actual: usize) -> Self {
        Self::InsufficientStencil {
            required,
            order,
            actual,
        }
    }

    /// Create unsupported order error
    pub fn unsupported_order(order: usize) -> Self {
        Self::UnsupportedOrder(order)
    }

    /// Create Robin singularity error
    pub fn robin_singularity(value: f64) -> Self {
        Self::RobinSingularity { value }
    }
}
