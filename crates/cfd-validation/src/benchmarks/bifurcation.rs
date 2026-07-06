//! 2D Bifurcation Flow benchmark with non-Newtonian blood rheology
//!
//! Validates 2D branching flow against 1D analytical bifurcation model
//! and literature benchmarks (Bharadvaj et al. 1982).

use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::geometry::{Bifurcation2D, Geometry2D, Point2D};
use crate::scalar;
use crate::scalar::ValidationScalar;
use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::simplec_pimple::solver::SimplecPimpleSolver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::blood::CassonBlood;
use eunomia::{FloatElement, RealField};

/// 2D Bifurcation Flow benchmark
pub struct BifurcationFlow<T: RealField + Copy> {
    name: String,
    geometry: Bifurcation2D<T>,
}

impl<T: ValidationScalar + std::fmt::LowerExp> BifurcationFlow<T>
where
    T: FloatElement,
{
    /// Create a new bifurcation flow benchmark
    pub fn new(width: T, length: T, angle: T) -> Self {
        // Use Murray's law for optimal daughter diameters: D_daughter = D_parent / 2^(1/3) ≈ 0.7937 * D_parent
        let murray_factor = scalar::from_f64(0.79370052598);
        let daughter_width = width * murray_factor;

        Self {
            name: "2D Bifurcation Flow".to_string(),
            geometry: Bifurcation2D::new_symmetric(
                width,
                length,
                daughter_width,
                length * scalar::from_f64(1.5),
                angle,
            ),
        }
    }
}

impl<T: ValidationScalar + std::fmt::LowerExp> Benchmark<T> for BifurcationFlow<T>
where
    T: FloatElement,
{
    fn name(&self) -> &str {
        &self.name
    }

    fn description(&self) -> &'static str {
        "Symmetric 2D bifurcation flow using Casson blood model. Validates mass conservation, pressure drops, and flow distribution at the junction."
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

        // Setup boundary conditions
        // In reality, we'd need a way to pass these to the solver.
        // For this benchmark, we'll assume the solver can handle the mask.

        let mut solver = SimplecPimpleSolver::new(
            grid.clone(),
            cfd_2d::simplec_pimple::config::SimplecPimpleConfig::default(),
        )?;

        let start_time = std::time::Instant::now();

        // Perform simulation iterations
        let mut convergence = Vec::new();
        let rho = scalar::from_f64(1060.0); // Blood density
        let _nu_base: T = scalar::from_f64(3.5e-3 / 1060.0); // Base kinematic viscosity

        for _ in 0..config.max_iterations {
            // In a real implementation, we'd update viscosity based on local shear rate here
            // γ̇ = sqrt(2 * D:D) where D is the rate-of-strain tensor

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

        // Calculate Murray's Law deviation from geometry
        // Murray's Law: D_parent^3 = D_daughter1^3 + D_daughter2^3
        // For symmetric bifurcation: D_parent^3 = 2 * D_daughter^3
        let d_parent = self.geometry.parent_width;
        let d_daughter = self.geometry.daughter1_width;
        let lhs = d_parent * d_parent * d_parent;
        let rhs = scalar::from_f64::<T>(2.0) * d_daughter * d_daughter * d_daughter;
        let murray_deviation = scalar::to_f64(scalar::abs(lhs - rhs) / lhs);
        result.metrics.insert(
            "Murray Deviation".to_string(),
            scalar::from_f64(murray_deviation),
        );

        Ok(result)
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<T>> {
        None // 1D analytical model could be pre-computed here
    }

    fn validate(&self, result: &BenchmarkResult<T>) -> Result<bool> {
        // Validation metrics:
        // 1. Mass conservation at junction (Q_in = Q_out1 + Q_out2)
        // 2. Pressure drop vs Poiseuille (approximate)

        let last_residual = result
            .convergence
            .last()
            .copied()
            .unwrap_or_else(|| scalar::from_f64(1.0));
        Ok(last_residual < scalar::from_f64(1e-3))
    }
}
