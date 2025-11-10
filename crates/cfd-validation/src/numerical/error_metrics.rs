//! Error metrics computation for validation

use nalgebra::{ComplexField, DVector, RealField};
use cfd_core::conversion::SafeFromF64;

/// Error metrics for validation
#[derive(Debug, Clone)]
pub struct ErrorMetrics<T: RealField + Copy> {
    /// L2 norm of error
    pub l2_error: T,
    /// L∞ norm of error
    pub linf_error: T,
    /// Relative L2 error
    pub relative_l2_error: T,
    /// Root mean square error
    pub rmse: T,
}

/// Compute error metrics between computed and analytical solutions
pub fn compute_error_metrics<T: RealField + Copy>(
    computed: &DVector<T>,
    analytical: &DVector<T>,
) -> ErrorMetrics<T> {
    assert_eq!(
        computed.len(),
        analytical.len(),
        "Solution vectors must have the same length"
    );

    let error = computed - analytical;
    let n = T::from_usize(computed.len()).unwrap_or(T::one());

    let l2_error = error.norm();
    let linf_error = error.amax();

    let analytical_norm = analytical.norm();
    let relative_l2_error = if analytical_norm > <T as SafeFromF64>::from_f64_or_zero(1e-12) {
        l2_error / analytical_norm
    } else {
        l2_error
    };

    let rmse = ComplexField::sqrt(l2_error * l2_error / n);

    ErrorMetrics {
        l2_error,
        linf_error,
        relative_l2_error,
        rmse,
    }
}
