//! SIMD optimizations for numerical operations
//!
//! Delegates to hermes-simd for generic runtime-dispatched SIMD operations.
//! hermes handles architecture detection (AVX2/SSE4.2/NEON/Scalar) internally,
//! eliminating the prior hand-rolled x86/arm/SWAR dispatch duplication.

pub mod cfd;
mod ops;
pub mod vector;
pub mod vectorization;

#[cfg(test)]
mod tests;

pub use ops::SimdOps;
