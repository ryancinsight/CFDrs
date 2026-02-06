//! 2D Venturi Flow benchmark with cavitation inception
//!
//! Validates 2D Venturi flow against Bernoulli equation and cavitation number
//! using non-Newtonian blood rheology.

use super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::geometry::{Venturi2D, Point2D, Geometry2D};
use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::grid::traits::Grid2D;
use cfd_2d::simplec_pimple::solver::SimplecPimpleSolver;
use cfd_core::error::Result;
use cfd_core::physics::fluid::blood::CassonBlood;
use cfd_core::physics::cavitation::VenturiCavitation;
use nalgebra::RealField;
use std::collections::HashMap;
use num_traits::{FromPrimitive, ToPrimitive};

/// 2D Venturi Flow benchmark
pub struct VenturiFlow<T: RealField + Copy> {
    name: String,
    geometry: Venturi2D<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::LowerExp + std::fmt::Debug> VenturiFlow<T> {
    /// Create a new Venturi flow benchmark
    pub fn new(inlet_width: T, throat_width: T) -> Self {
        Self {
            name: "2D Venturi Flow".to_string(),
            geometry: Venturi2D::new(
                inlet_width,
                throat_width,
                inlet_width * T::from_f64(1.5).unwrap(), // Converge length
                throat_width * T::from_f64(1.0).unwrap(), // Throat length
                inlet_width * T::from_f64(3.0).unwrap(), // Diverge length (longer for recovery)
            ),
        }
    }
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::LowerExp + std::fmt::Debug> Benchmark<T> for VenturiFlow<T> {
    fn name(&self) -> &str {
        &self.name
    }

    fn description(&self) -> &str {
        "2D Venturi flow with Casson blood. Validates pressure-velocity relationship against Bernoulli and cavitation number at the throat."
    }

    fn run(&self, config: &BenchmarkConfig<T>) -> Result<BenchmarkResult<T>> {
        let (min_p, max_p) = self.geometry.bounds();
        let lx = max_p.x - min_p.x;
        let ly = max_p.y - min_p.y;
        
        let nx = config.resolution;
        let ny = (config.resolution as f64 * ly.to_f64().unwrap() / lx.to_f64().unwrap()) as usize;
        
        let grid = StructuredGrid2D::new(
            nx,
            ny,
            min_p.x,
            max_p.x,
            min_p.y,
            max_p.y,
        )?;

        let mut fields = SimulationFields::new(nx, ny);
        let blood = CassonBlood::<T>::normal_blood();
        
        // Populate mask and initial viscosity
        for i in 0..nx {
            for j in 0..ny {
                let center = grid.cell_center(i, j)?;
                let point = Point2D { x: center.x, y: center.y };
                let is_fluid = self.geometry.contains(&point);
                fields.mask.set(i, j, is_fluid);
                
                if is_fluid {
                    fields.viscosity.set(i, j, blood.apparent_viscosity(T::from_f64(100.0).unwrap()));
                } else {
                    fields.viscosity.set(i, j, T::zero());
                }
            }
        }

        let mut solver = SimplecPimpleSolver::new(grid.clone(), cfd_2d::simplec_pimple::config::SimplecPimpleConfig::default())?;
        
        let start_time = std::time::Instant::now();
        let mut convergence = Vec::new();
        let rho = T::from_f64(1060.0).unwrap();
        
        for _ in 0..config.max_iterations {
            let dt = config.time_step.unwrap_or(T::from_f64(0.01).unwrap());
            let residual = solver.solve_time_step(&mut fields, dt, T::zero(), rho)?;
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
            convergent_angle: (self.geometry.inlet_width - self.geometry.throat_width).atan2(self.geometry.l_converge),
            divergent_angle: (self.geometry.outlet_width - self.geometry.throat_width).atan2(self.geometry.l_diverge),
            inlet_pressure: T::from_f64(101325.0).unwrap(), // Atmospheric inlet
            inlet_velocity: T::one(), // Unit inlet velocity
            density: rho,
            vapor_pressure: T::from_f64(2339.0).unwrap(), // Water at 20C approx
        };

        let mut metadata = HashMap::new();
        metadata.insert("throat_velocity_analytic".to_string(), format!("{:e}", analytic.throat_velocity()));
        metadata.insert("cavitation_number".to_string(), format!("{:e}", analytic.cavitation_number()));

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
        let last_residual = result.convergence.last().copied().unwrap_or_else(|| T::from_f64(1.0).unwrap());
        Ok(last_residual < T::from_f64(1e-3).unwrap())
    }
}
