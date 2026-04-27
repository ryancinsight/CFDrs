//! Performance benchmark validation methods for [`TurbulenceValidator`].

use super::ValidationResult;
use crate::physics::turbulence::traits::LESTurbulenceModel;
use crate::physics::turbulence::{DetachedEddySimulation, SmagorinskyLES};
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

use super::TurbulenceValidator;

/// Supported turbulence models for validation benchmarks.
///
/// # Theorem - Exhaustive Benchmark Dispatch
///
/// The performance benchmark contains no placeholder execution path when model
/// selection is represented by this closed enum. Every constructed
/// [`BenchmarkModel`] maps to one concrete production update loop, and parsing
/// rejects unsupported names before dispatch.
///
/// **Proof sketch**: `run_performance_benchmark` matches all variants of
/// `BenchmarkModel`. Rust checks the match exhaustively. Since the string parser
/// returns `None` for unsupported names, the benchmark either runs a concrete
/// supported model or returns a typed rejection; it cannot enter an
/// unimplemented model branch.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BenchmarkModel {
    SmagorinskyLes,
    DetachedEddySimulation,
}

impl BenchmarkModel {
    fn parse(model_name: &str) -> Option<Self> {
        match model_name {
            "smagorinsky-les" => Some(Self::SmagorinskyLes),
            "des" => Some(Self::DetachedEddySimulation),
            _ => None,
        }
    }

    fn name(self) -> &'static str {
        match self {
            Self::SmagorinskyLes => "smagorinsky-les",
            Self::DetachedEddySimulation => "des",
        }
    }
}

impl<T: RealField + FromPrimitive + ToPrimitive + Copy> TurbulenceValidator<T> {
    /// Validate turbulence model performance
    pub fn validate_model_performance(&self, model_name: &str) -> ValidationResult {
        let Some(model) = BenchmarkModel::parse(model_name) else {
            return ValidationResult {
                test_name: format!("{model_name} Performance"),
                passed: false,
                metric: "Unsupported turbulence benchmark model".to_string(),
                details: "Supported models: smagorinsky-les, des".to_string(),
            };
        };
        self.run_performance_benchmark(model)
    }

    fn run_performance_benchmark(&self, model: BenchmarkModel) -> ValidationResult {
        let start_time = std::time::Instant::now();

        let iterations = 100;
        let nx = 32;
        let ny = 32;

        match model {
            BenchmarkModel::SmagorinskyLes => {
                let config = crate::physics::turbulence::les_smagorinsky::SmagorinskyConfig {
                    smagorinsky_constant: 0.1,
                    dynamic_procedure: false,
                    wall_damping: false,
                    van_driest_constant: 0.0,
                    min_sgs_viscosity: 1e-10,
                    use_gpu: false,
                };
                let mut model = SmagorinskyLES::new(nx, ny, 0.1, 0.1, config);

                let velocity_u = nalgebra::DMatrix::from_element(nx, ny, 1.0);
                let velocity_v = nalgebra::DMatrix::from_element(nx, ny, 0.0);
                let pressure = nalgebra::DMatrix::zeros(nx, ny);

                for _ in 0..iterations {
                    model
                        .update(
                            &velocity_u,
                            &velocity_v,
                            &pressure,
                            1.0,
                            1e-5,
                            0.001,
                            0.1,
                            0.1,
                        )
                        .unwrap();
                }
            }
            BenchmarkModel::DetachedEddySimulation => {
                let config = crate::physics::turbulence::des::DESConfig {
                    des_constant: 0.65,
                    max_sgs_ratio: 0.5,
                    variant: crate::physics::turbulence::des::DESVariant::DES97,
                    rans_viscosity: 1e-5,
                    use_gpu: false,
                };
                let mut model = DetachedEddySimulation::new(nx, ny, 0.1, 0.1, config, &[]);

                let velocity_u = nalgebra::DMatrix::from_element(nx, ny, 1.0);
                let velocity_v = nalgebra::DMatrix::from_element(nx, ny, 0.0);
                let pressure = nalgebra::DMatrix::zeros(nx, ny);

                for _ in 0..iterations {
                    model
                        .update(
                            &velocity_u,
                            &velocity_v,
                            &pressure,
                            1.0,
                            1e-5,
                            0.001,
                            0.1,
                            0.1,
                        )
                        .unwrap();
                }
            }
        }

        let elapsed = start_time.elapsed();
        let operations_per_sec = (iterations * nx * ny) as f64 / elapsed.as_secs_f64();

        ValidationResult {
            test_name: format!("{} Performance", model.name()),
            passed: operations_per_sec > 1000.0,
            metric: format!("{operations_per_sec:.0} ops/sec"),
            details: format!("Grid: {nx}x{ny}, Iterations: {iterations}"),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn benchmark_model_parser_accepts_only_supported_models() {
        assert_eq!(
            BenchmarkModel::parse("smagorinsky-les"),
            Some(BenchmarkModel::SmagorinskyLes)
        );
        assert_eq!(
            BenchmarkModel::parse("des"),
            Some(BenchmarkModel::DetachedEddySimulation)
        );
        assert_eq!(BenchmarkModel::parse("unsupported"), None);
    }

    #[test]
    fn unsupported_benchmark_returns_typed_rejection_without_placeholder_text() {
        let validator = TurbulenceValidator::<f64>::new(1.0e-12);
        let result = validator.validate_model_performance("unsupported");

        assert!(!result.passed);
        assert_eq!(result.test_name, "unsupported Performance");
        assert_eq!(
            result.metric,
            "Unsupported turbulence benchmark model".to_string()
        );
        assert_eq!(result.details, "Supported models: smagorinsky-les, des");
        assert!(!result.details.contains("not implemented"));
    }
}
