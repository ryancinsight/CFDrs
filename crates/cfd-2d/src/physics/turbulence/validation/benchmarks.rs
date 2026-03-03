//! Performance benchmark validation methods for [`TurbulenceValidator`].

use super::ValidationResult;
use crate::physics::turbulence::traits::LESTurbulenceModel;
use crate::physics::turbulence::{DetachedEddySimulation, SmagorinskyLES};
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

use super::TurbulenceValidator;

impl<T: RealField + FromPrimitive + ToPrimitive + Copy> TurbulenceValidator<T> {
    /// Validate turbulence model performance
    pub fn validate_model_performance(&self, model_name: &str) -> ValidationResult {
        let start_time = std::time::Instant::now();

        let iterations = 100;
        let nx = 32;
        let ny = 32;

        match model_name {
            "smagorinsky-les" => {
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
            "des" => {
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
            _ => {
                return ValidationResult {
                    test_name: format!("{model_name} Performance"),
                    passed: false,
                    metric: "Unknown model".to_string(),
                    details: "Performance test not implemented".to_string(),
                }
            }
        }

        let elapsed = start_time.elapsed();
        let operations_per_sec = (iterations * nx * ny) as f64 / elapsed.as_secs_f64();

        ValidationResult {
            test_name: format!("{model_name} Performance"),
            passed: operations_per_sec > 1000.0,
            metric: format!("{operations_per_sec:.0} ops/sec"),
            details: format!("Grid: {nx}x{ny}, Iterations: {iterations}"),
        }
    }
}
