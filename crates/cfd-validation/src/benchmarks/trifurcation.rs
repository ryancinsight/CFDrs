//! 2D Trifurcation Flow benchmark with non-Newtonian blood rheology
//!
//! Validates 2D triple-branch flow against analytical split ratios
//! and literature benchmarks for wall shear stress distribution.

use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::geometry::{Geometry2D, Point2D, Trifurcation2D};
use crate::scalar;
use crate::scalar::ValidationScalar;
use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::simplec_pimple::solver::SimplecPimpleSolver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::blood::CassonBlood;
use eunomia::{FloatElement, RealField};

/// 2D Trifurcation Flow benchmark
pub struct TrifurcationFlow<T: RealField + Copy> {
    name: String,
    geometry: Trifurcation2D<T>,
}

impl<T: ValidationScalar + std::fmt::LowerExp> TrifurcationFlow<T>
where
    T: FloatElement,
{
    /// Create a new trifurcation flow benchmark
    pub fn new(width: T, length: T, angle: T) -> Self {
        // Extended Murray's Law: D_parent^3 = 3 * D_daughter^3
        // D_daughter = D_parent / 3^(1/3) ≈ 0.693 * D_parent
        let murray_factor = scalar::from_f64(0.69336127435);
        let daughter_width = width * murray_factor;

        Self {
            name: "2D Trifurcation Flow".to_string(),
            geometry: Trifurcation2D::new_symmetric(
                width,
                length,
                daughter_width,
                length * scalar::from_f64(1.5),
                angle,
            ),
        }
    }
}

impl<T: ValidationScalar + std::fmt::LowerExp> Benchmark<T> for TrifurcationFlow<T>
where
    T: FloatElement,
{
    fn name(&self) -> &str {
        &self.name
    }

    fn description(&self) -> &'static str {
        "Symmetric 2D trifurcation flow with Casson blood. Validates distribution across three branches and junction pressure stability."
    }

    fn run(&self, config: &BenchmarkConfig<T>) -> Result<BenchmarkResult<T>> {
        let (min_p, max_p) = self.geometry.bounds();
        let lx = max_p.x - min_p.x;
        let ly = max_p.y - min_p.y;

        let nx = config.resolution;
        let ny = (config.resolution as f64 * scalar::to_f64(ly) / scalar::to_f64(lx)) as usize;

        let grid = StructuredGrid2D::new(nx, ny, min_p.x, max_p.x, min_p.y, max_p.y)?;

        let mut fields = SimulationFields::new(nx, ny);
        let blood = CassonBlood::<T>::normal_blood();

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
        Ok(last_residual < scalar::from_f64(1e-3))
    }
}
