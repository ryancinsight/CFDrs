//! Error metrics computation for validation

use crate::scalar;
use eunomia::{FloatElement, RealField};
use leto::Array1;

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
pub fn compute_error_metrics<T: RealField + Copy + FloatElement>(
    computed: &Array1<T>,
    analytical: &Array1<T>,
) -> ErrorMetrics<T> {
    assert_eq!(
        computed.shape()[0],
        analytical.shape()[0],
        "Solution vectors must have the same length"
    );

    let mut l2_sum = scalar::zero::<T>();
    let mut linf_error = scalar::zero::<T>();
    let mut analytical_l2_sum = scalar::zero::<T>();

    for idx in 0..computed.shape()[0] {
        let error = scalar::abs(computed[idx] - analytical[idx]);
        l2_sum += error * error;
        analytical_l2_sum += analytical[idx] * analytical[idx];
        if error > linf_error {
            linf_error = error;
        }
    }

    let n = scalar::from_usize::<T>(computed.shape()[0]);
    let l2_error = scalar::sqrt(l2_sum);
    let analytical_norm = scalar::sqrt(analytical_l2_sum);
    let relative_l2_error = if analytical_norm > scalar::from_f64::<T>(1e-12) {
        l2_error / analytical_norm
    } else {
        l2_error
    };

    let rmse = scalar::sqrt(l2_error * l2_error / n);

    ErrorMetrics {
        l2_error,
        linf_error,
        relative_l2_error,
        rmse,
    }
}
