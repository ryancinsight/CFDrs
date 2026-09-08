//! SIMD compute backend with architecture-specific dispatch

#[cfg(target_arch = "aarch64")]
pub mod aarch64;
#[cfg(any(target_arch = "x86", target_arch = "x86_64"))]
pub mod x86;

// SIMD implementations moved to architecture-specific modules.
// Scalar fallback is used where SIMD performance is lower than scalar.
