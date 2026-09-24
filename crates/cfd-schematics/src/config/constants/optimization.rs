//! Optimization algorithm constants for solver tuning

/// Optimization algorithm constants previously hardcoded
///
/// The three factor tables are borrowed slices rather than owned `Vec`s: they
/// are fixed for the lifetime of the program, so owning them only forced an
/// allocation per registry construction and a clone per accessor call.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct OptimizationConstants {
    /// Branch factor scaling exponent (was hardcoded as 0.75)
    pub branch_factor_exponent: f64,

    /// Fill factor enhancement multiplier (was hardcoded as 1.5)
    pub fill_factor_enhancement: f64,

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
        branch_factor_exponent: 0.75,
        fill_factor_enhancement: 1.5,
        max_optimization_iterations: 100,
        convergence_tolerance: 1e-6,
        fast_wavelength_factors: &[1.0, 2.0, 3.0, 4.0],
        fast_wave_density_factors: &[1.0, 2.0, 3.0],
        fast_fill_factors: &[0.7, 0.8, 0.9],
    };
}

impl Default for OptimizationConstants {
    fn default() -> Self {
        Self::DEFAULT
    }
}
