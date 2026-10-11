//! SOR conformance tests for the provider preconditioner.
//!
//! The saddle-point and multigrid families live in `athena-leto` (Atlas ADR
//! 0062); the consumer-facing re-exports retired with them. This module keeps
//! only the SOR contract tests that pin the provider's relaxation behavior
//! from the consumer side.

#[cfg(test)]
mod ssor_tests;
