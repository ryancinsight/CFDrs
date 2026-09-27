//! Optimization algorithm constants for solver tuning

use super::primitives;

/// Optimization algorithm constants previously hardcoded
///
/// The three factor tables are borrowed slices rather than owned `Vec`s: they
/// are fixed for the lifetime of the program, so owning them only forced an
/// allocation per registry construction and a clone per accessor call. The
/// values themselves live in [`primitives`].
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OptimizationConstants {
    /// Branch factor scaling exponent (was hardcoded as 0.75)
    pub branch_factor_exponent: f64,

    /// Maximum optimization iterations
    pub max_optimization_iterations: usize,

    /// Convergence tolerance for optimization
    pub convergence_tolerance: f64,

    /// Fast optimization wavelength factors
    pub fast_wavelength_factors: &'static [f64],

    /// Fast optimization wave density factors
    pub fast_wave_density_factors: &'static [f64],

    /// Fast optimization fill factors
    pub fast_fill_factors: &'static [f64],
}

impl OptimizationConstants {
    /// Canonical default optimization constants.
    pub const DEFAULT: Self = Self {
        branch_factor_exponent: primitives::BRANCH_FACTOR_EXPONENT,
        max_optimization_iterations: primitives::MAX_OPTIMIZATION_ITERATIONS,
        convergence_tolerance: primitives::CONVERGENCE_TOLERANCE,
        fast_wavelength_factors: primitives::FAST_WAVELENGTH_FACTORS,
        fast_wave_density_factors: primitives::FAST_WAVE_DENSITY_FACTORS,
        fast_fill_factors: primitives::FAST_FILL_FACTORS,
    };
}

impl Default for OptimizationConstants {
    fn default() -> Self {
        Self::DEFAULT
    }
}
