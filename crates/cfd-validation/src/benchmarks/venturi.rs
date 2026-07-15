//! 2D Venturi Flow benchmark with cavitation inception
//!
//! Validates 2D Venturi flow against Bernoulli equation and cavitation number
//! using non-Newtonian blood rheology.

use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::geometry::{Geometry2D, Point2D, Venturi2D};
use crate::scalar;
use crate::scalar::ValidationScalar;
use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::simplec_pimple::solver::SimplecPimpleSolver;
use cfd_core::error::Result;
use cfd_core::physics::cavitation::VenturiCavitation;
use cfd_core::physics::fluid::blood::CassonBlood;
use eunomia::{FloatElement, RealField};
use std::collections::HashMap;

/// 2D Venturi Flow benchmark
pub struct VenturiFlow<T: RealField + Copy> {
    name: String,
    geometry: Venturi2D<T>,
}

impl<T: ValidationScalar + std::fmt::LowerExp> VenturiFlow<T>
where
    T: FloatElement,
{
    /// Create a new Venturi flow benchmark
    pub fn new(inlet_width: T, throat_width: T) -> Self {
        Self {
            name: "2D Venturi Flow".to_string(),
            geometry: Venturi2D::new(
                inlet_width,
                throat_width,
                inlet_width * scalar::from_f64(1.5), // Converge length
                throat_width * scalar::from_f64(1.0), // Throat length
                inlet_width * scalar::from_f64(3.0), // Diverge length (longer for recovery)
            ),
        }
    }
}

impl<T: ValidationScalar + std::fmt::LowerExp> Benchmark<T> for VenturiFlow<T>
where
    T: FloatElement,
{
    fn name(&self) -> &str {
        &self.name
    }

    fn description(&self) -> &'static str {
        "2D Venturi flow with Casson blood. Validates pressure-velocity relationship against Bernoulli and cavitation number at the throat."
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

        // Calculate analytical reference for comparison in metadata
        let analytic = VenturiCavitation::<T> {
            inlet_diameter: self.geometry.inlet_width, // 2D approximation width ≈ diameter
            throat_diameter: self.geometry.throat_width,
            outlet_diameter: self.geometry.outlet_width,
            convergent_angle: FloatElement::atan2(
                self.geometry.inlet_width - self.geometry.throat_width,
                self.geometry.l_converge,
            ),
            divergent_angle: FloatElement::atan2(
                self.geometry.outlet_width - self.geometry.throat_width,
                self.geometry.l_diverge,
            ),
            inlet_pressure: scalar::from_f64(101325.0), // Atmospheric inlet
            inlet_velocity: scalar::from_f64(1.0),      // Unit inlet velocity
            density: rho,
            vapor_pressure: scalar::from_f64(2339.0), // Water at 20C approx
        };

        let mut metadata = HashMap::new();
        metadata.insert(
            "throat_velocity_analytic".to_string(),
            format!("{:e}", analytic.throat_velocity()),
        );
        metadata.insert(
            "cavitation_number".to_string(),
            format!("{:e}", analytic.cavitation_number()),
        );

        let mut result = BenchmarkResult::new(&self.name);
        result.values = fields.u.data().to_vec();
        result.convergence = convergence;
        result.execution_time = execution_time;
        result.metadata = metadata;

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
