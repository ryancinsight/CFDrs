//! Tests for k-ω SST model components.

use super::limiter::limit_production;
use crate::physics::turbulence::constants::SST_BETA_STAR;

/// When P_k is smaller than the limiter threshold (10 * beta_star * k * omega),
/// the production term should pass through unchanged.
#[test]
fn test_sst_limiter_passes_small_production() {
    let beta_star = SST_BETA_STAR; // 0.09
    let k = 1.0;
    let omega = 100.0;
    // Limit = 10 * 0.09 * 1.0 * 100.0 = 90.0
    let p_k = 50.0; // well below limit
    let result = limit_production(p_k, k, omega, beta_star);
    assert!(
        (result - p_k).abs() < 1e-12,
        "Small production should pass through unchanged: got {result}, expected {p_k}"
    );
}

/// When P_k exceeds the limiter threshold, it should be clipped to
/// C_lim * beta_star * k * omega.
#[test]
fn test_sst_limiter_clips_large_production() {
    let beta_star = SST_BETA_STAR; // 0.09
    let k = 1.0;
    let omega = 100.0;
    let limit = 10.0 * beta_star * k * omega; // 90.0
    let p_k = 500.0; // well above limit
    let result = limit_production(p_k, k, omega, beta_star);
    assert!(
        (result - limit).abs() < 1e-12,
        "Large production should be clipped to {limit}: got {result}"
    );
}

/// When k = 0, the limiter should clip production to zero.
#[test]
fn test_sst_limiter_zero_k() {
    let beta_star = SST_BETA_STAR;
    let k = 0.0;
    let omega = 100.0;
    let p_k = 50.0;
    let result = limit_production(p_k, k, omega, beta_star);
    assert!(
        result.abs() < 1e-12,
        "With k=0, limiter should clip to zero: got {result}"
    );
}
