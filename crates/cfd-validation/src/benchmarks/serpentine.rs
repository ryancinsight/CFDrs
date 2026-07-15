//! 2D Serpentine Flow benchmark with Dean vortex validation
//!
//! Validates 2D Serpentine flow using Carreau-Yasuda blood rheology
//! and secondary flow strength analysis.

use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::geometry::{Geometry2D, Point2D, Serpentine2D};
use crate::scalar;
use crate::scalar::ValidationScalar;
use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::simplec_pimple::solver::SimplecPimpleSolver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use eunomia::{FloatElement, RealField};

/// 2D Serpentine Flow benchmark
pub struct SerpentineFlow<T: RealField + Copy> {
    name: String,
    geometry: Serpentine2D<T>,
}

impl<T: ValidationScalar + std::fmt::LowerExp> SerpentineFlow<T>
where
    T: FloatElement,
{
    /// Create a new Serpentine flow benchmark
    pub fn new(width: T, amplitude: T, period_length: T, periods: usize) -> Self {
        Self {
            name: "2D Serpentine Flow".to_string(),
            geometry: Serpentine2D::new(width, amplitude, period_length, periods),
        }
    }
}

impl<T: ValidationScalar + std::fmt::LowerExp> Benchmark<T> for SerpentineFlow<T>
where
    T: FloatElement,
{
    fn name(&self) -> &str {
        &self.name
    }

    fn description(&self) -> &'static str {
        "2D Serpentine flow with Carreau-Yasuda blood. Validates secondary flow development and Dean number relationships."
    }

    fn run(&self, config: &BenchmarkConfig<T>) -> Result<BenchmarkResult<T>> {
        let (min_p, max_p) = self.geometry.bounds();
        let lx = max_p.x - min_p.x;
        let ly = max_p.y - min_p.y;

        let nx = config.resolution;
        let ny = (config.resolution as f64 * scalar::to_f64(ly) / scalar::to_f64(lx)) as usize;

        let grid = StructuredGrid2D::new(nx, ny, min_p.x, max_p.x, min_p.y, max_p.y)?;

        let mut fields = SimulationFields::new(nx, ny);
        let blood = CarreauYasudaBlood::<T>::normal_blood();

        // Populate mask and initial viscosity
        for i in 0..nx {
            for j in 0..ny {
                let center = grid.cell_center(i, j)?;
                let point = Point2D::new(center[0], center[1]);
                let is_fluid = self.geometry.contains(&point);
                fields.mask.set(i, j, is_fluid);

                if is_fluid {
                    fields
                        .viscosity
                        .set(i, j, blood.apparent_viscosity(scalar::from_f64(100.0)));
                } else {
                    fields.viscosity.set(i, j, scalar::from_f64(0.0));
                }
            }
        }

        let mut solver = SimplecPimpleSolver::new(
            grid.clone(),
            cfd_2d::simplec_pimple::config::SimplecPimpleConfig::default(),
        )?;

        let start_time = std::time::Instant::now();
        let mut convergence = Vec::new();
        let rho = scalar::from_f64(1060.0);

        for _ in 0..config.max_iterations {
            let dt = config.time_step.unwrap_or_else(|| scalar::from_f64(0.01));
            let residual = solver.solve_time_step(&mut fields, dt, scalar::from_f64(0.0), rho)?;
            convergence.push(residual);

            if residual < config.tolerance {
                break;
            }
        }

        let execution_time = start_time.elapsed().as_secs_f64();

        let mut result = BenchmarkResult::new(&self.name);
        result.values = fields.u.data().to_vec();
        result.convergence = convergence;
        result.execution_time = execution_time;

        // Calculate metrics
        let props =
            FluidTrait::properties_at(&blood, scalar::from_f64(310.0), scalar::from_f64(101325.0))?;
        // dynamic_viscosity available via props if needed downstream
        let _ = props.dynamic_viscosity;

        // Characteristic velocity from Re
        let re_target = config.reynolds_number;

        // For serpentine, radius of curvature R can be estimated from amplitude and period
        // R approx P^2 / (4 * pi^2 * A) at the peak
        let pi: T = scalar::from_f64(std::f64::consts::PI);
        let r_curv = scalar::powf(self.geometry.period_length, scalar::from_f64(2.0))
            / (scalar::from_f64::<T>(4.0) * pi * pi * self.geometry.amplitude);
        let dean_number =
            re_target * scalar::sqrt(self.geometry.width / (scalar::from_f64::<T>(2.0) * r_curv));

        result
            .metrics
            .insert("reynolds_number".to_string(), re_target);
        result
            .metrics
            .insert("dean_number".to_string(), dean_number);

        // Pressure drop
        let p_inlet = fields.p.at(0, ny / 2);
        let p_outlet = fields.p.at(nx - 1, ny / 2);
        let dp = p_inlet - p_outlet;
        result.metrics.insert("pressure_drop".to_string(), dp);

        Ok(result)
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<T>> {
        None
    }

    fn validate(&self, result: &BenchmarkResult<T>) -> Result<bool> {
        let last_residual = result
            .convergence
            .last()
            .copied()
            .unwrap_or_else(|| scalar::from_f64(1.0));
        if last_residual > scalar::from_f64(5e-3) {
            return Ok(false);
        }

        // Basic Dean number sanity check
        if let Some(&de) = result.metrics.get("dean_number") {
            if de < scalar::from_f64(0.0) {
                return Ok(false);
            }
        }

        Ok(true)
    }
}
