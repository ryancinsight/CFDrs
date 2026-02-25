//! Optimization algorithm constants for solver tuning

use crate::state_management::{ConfigurableParameter, ParameterConstraints, ParameterMetadata};

/// Optimization algorithm constants previously hardcoded
pub struct OptimizationConstants {
    /// Branch factor scaling exponent (was hardcoded as 0.75)
    pub branch_factor_exponent: ConfigurableParameter<f64>,

    /// Fill factor enhancement multiplier (was hardcoded as 1.5)
    pub fill_factor_enhancement: ConfigurableParameter<f64>,

    /// Maximum optimization iterations
    pub max_optimization_iterations: ConfigurableParameter<usize>,

    /// Convergence tolerance for optimization
    pub convergence_tolerance: ConfigurableParameter<f64>,

    /// Fast optimization wavelength factors
    pub fast_wavelength_factors: ConfigurableParameter<Vec<f64>>,

    /// Fast optimization wave density factors
    pub fast_wave_density_factors: ConfigurableParameter<Vec<f64>>,

    /// Fast optimization fill factors
    pub fast_fill_factors: ConfigurableParameter<Vec<f64>>,
}

impl OptimizationConstants {
    /// Create default optimization constants
    #[must_use]
    fn make_default() -> Self {
        Self {
            branch_factor_exponent: ConfigurableParameter::new(
                0.75,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(0.1, 2.0),
                ]),
                ParameterMetadata::new(
                    "branch_factor_exponent",
                    "Exponent for branch factor scaling (was hardcoded as 0.75)",
                    "optimization",
                )
                .with_units("exponent"),
            ),

            fill_factor_enhancement: ConfigurableParameter::new(
                1.5,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1.0, 3.0),
                ]),
                ParameterMetadata::new(
                    "fill_factor_enhancement",
                    "Enhancement multiplier for fill factor (was hardcoded as 1.5)",
                    "optimization",
                )
                .with_units("multiplier"),
            ),

            max_optimization_iterations: ConfigurableParameter::new(
                100usize,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<usize>::positive(),
                    ParameterConstraints::range(10, 1000),
                ]),
                ParameterMetadata::new(
                    "max_optimization_iterations",
                    "Maximum number of iterations for optimization algorithms",
                    "optimization",
                )
                .with_units("iterations"),
            ),

            convergence_tolerance: ConfigurableParameter::new(
                1e-6,
                ParameterConstraints::all(vec![
                    ParameterConstraints::<f64>::positive(),
                    ParameterConstraints::range(1e-10, 1e-2),
                ]),
                ParameterMetadata::new(
                    "convergence_tolerance",
                    "Tolerance for optimization convergence detection",
                    "optimization",
                )
                .with_units("tolerance"),
            ),

            fast_wavelength_factors: ConfigurableParameter::new(
                vec![1.0, 2.0, 3.0, 4.0],
                ParameterConstraints::all(vec![]),
                ParameterMetadata::new(
                    "fast_wavelength_factors",
                    "Wavelength factors for fast optimization (was hardcoded array)",
                    "optimization",
                )
                .with_units("factors"),
            ),

            fast_wave_density_factors: ConfigurableParameter::new(
                vec![1.0, 2.0, 3.0],
                ParameterConstraints::all(vec![]),
                ParameterMetadata::new(
                    "fast_wave_density_factors",
                    "Wave density factors for fast optimization (was hardcoded array)",
                    "optimization",
                )
                .with_units("factors"),
            ),

            fast_fill_factors: ConfigurableParameter::new(
                vec![0.7, 0.8, 0.9],
                ParameterConstraints::all(vec![]),
                ParameterMetadata::new(
                    "fast_fill_factors",
                    "Fill factors for fast optimization (was hardcoded array)",
                    "optimization",
                )
                .with_units("factors"),
            ),
        }
    }
}

impl Default for OptimizationConstants {
    fn default() -> Self {
        Self::make_default()
    }
}
