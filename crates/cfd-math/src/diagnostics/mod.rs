//! Adaptive performance diagnostics for CFD numerical operations.
//!
//! ## Architectural Role
//!
//! This module is the **Single Source of Truth** for runtime performance measurement,
//! adaptive threshold tuning, and algorithm selection diagnostics. It enforces SRP:
//! all monitoring concerns are isolated here, away from numerical algorithms.
//!
//! ## Theorem — Adaptive Threshold Optimality
//!
//! **Theorem (Piecewise Linear Crossover)**: Given throughput measurements T(n) for
//! operation sizes n, the optimal SIMD crossover threshold n* minimises total
//! Residual Sum of Squares across a bilinear fit:
//!
//! ```text
//! n* = argmin_{k} [ RSS(T[1..k]) + RSS(T[k+1..N]) ]
//! ```
//!
//! where RSS(S) = Σ(T_i - (α·nᵢ + β))². The algorithm runs in O(N²) with N
//! data points; for the typical calibration grid of 8 sizes this is negligible.
//!
//! ## Invariants
//!
//! - `PerformanceMetrics::avg_time_ns` is a running mean — never decreases below 0.
//! - `cache_efficiency` ∈ [0, 1] — enforced via `.clamp(0.0, 1.0)`.
//! - `AdaptivePerformanceMonitor` is `Arc<Mutex<...>>`-guarded — safe to share across threads.

pub mod performance_monitor;

pub use performance_monitor::{AdaptivePerformanceMonitor, PerformanceAwareSelector, PerformanceMetrics};
