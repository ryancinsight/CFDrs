//! Analytical benchmarks for CFD validation
//!
//! Exact reference data used to validate numerical implementations.
//!
//! Analytical solution implementations live in [`crate::analytical`]. This
//! module retains benchmark data whose source contract is a normalized sample
//! table rather than a physical-unit solution model.

/// Lid-driven cavity benchmark points
///
/// Reference: Ghia et al. (1982). J. Comput. Phys., 48(3), 387-411
pub mod lid_driven_cavity {
    /// Benchmark data for Re=100
    pub const RE100_U_CENTERLINE: &[(f64, f64)] = &[
        (0.0000, 0.00000),
        (0.0625, -0.03717),
        (0.1250, -0.04192),
        (0.1875, -0.04775),
        (0.2500, -0.05641),
        (0.3125, -0.04299),
        (0.3750, -0.02388),
        (0.4375, -0.00737),
        (0.5000, 0.00669),
        (0.5625, 0.01911),
        (0.6250, 0.03055),
        (0.6875, 0.04108),
        (0.7500, 0.05052),
        (0.8125, 0.05864),
        (0.8750, 0.06500),
        (0.9375, 0.06898),
        (1.0000, 1.00000),
    ];

    /// Benchmark data for Re=1000
    pub const RE1000_U_CENTERLINE: &[(f64, f64)] = &[
        (0.0000, 0.00000),
        (0.0625, -0.18109),
        (0.1250, -0.20196),
        (0.1875, -0.22220),
        (0.2500, -0.24533),
        (0.3125, -0.14612),
        (0.3750, -0.10150),
        (0.4375, -0.06547),
        (0.5000, -0.03827),
        (0.5625, -0.01860),
        (0.6250, -0.00570),
        (0.6875, 0.00190),
        (0.7500, 0.00570),
        (0.8125, 0.00760),
        (0.8750, 0.00950),
        (0.9375, 0.01040),
        (1.0000, 1.00000),
    ];
}
