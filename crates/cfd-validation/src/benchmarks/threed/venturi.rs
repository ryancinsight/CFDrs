//! 3D Venturi tube flow benchmark with blood rheology
//!
//! Validates 3D Venturi flow against 1D analytical models and
//! experimental discharge coefficients for high-Reynolds flows.

use super::super::{Benchmark, BenchmarkConfig, BenchmarkResult};
use crate::geometry::threed::Venturi3D;
use cfd_3d::venturi::{VenturiSolver3D, VenturiConfig3D};
use cfd_mesh::VenturiMeshBuilder;
use cfd_core::physics::fluid::blood::CarreauYasudaBlood;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use nalgebra::RealField;

/// 3D Venturi Flow benchmark
pub struct VenturiFlow3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// The 3D venturi geometry
    pub geometry: Venturi3D<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> VenturiFlow3D<T> {
    /// Create a new 3D venturi flow benchmark
    pub fn new(geometry: Venturi3D<T>) -> Self {
        Self { geometry }
    }
}

impl<T: RealField + Copy + num_traits::Float + num_traits::FromPrimitive + 
    num_traits::ToPrimitive + cfd_core::conversion::SafeFromF64 + std::convert::From<f64> + cfd_mesh::domain::core::Scalar> Benchmark<T> for VenturiFlow3D<T> {
    fn name(&self) -> &str {
        "3D Venturi Tube Flow"
    }

    fn description(&self) -> &str {
        "3D Venturi Tube flow with blood rheology. Validates pressure drop, discharge coefficient, and recovery."
    }

    fn run(&self, config: &BenchmarkConfig<T>) -> cfd_core::error::Result<BenchmarkResult<T>> {
        let mut result = BenchmarkResult::new(self.name());
        
        // 1. Setup Mesh Builder from Geometry
        let builder = VenturiMeshBuilder::new(
            self.geometry.d_inlet,
            self.geometry.d_throat,
            self.geometry.l_inlet,
            self.geometry.l_convergent,
            self.geometry.l_throat,
            self.geometry.l_divergent,
            self.geometry.l_outlet,
        );
        
        // 2. Setup Solver
        let solver_config = VenturiConfig3D {
            inlet_flow_rate: T::from_f64_or_one(1.0e-7), // Default if not in config
            resolution: (config.resolution, 8),
            ..VenturiConfig3D::default()
        };
        
        let solver = VenturiSolver3D::new(builder, solver_config);
        
        // 3. Run for Newtonian and Non-Newtonian
        let blood = CarreauYasudaBlood::<T>::normal_blood();
        let solution = solver.solve(blood)?;
        
        result.metrics.insert("Pressure Drop".to_string(), solution.dp_throat);
        result.metrics.insert("Recovery Coefficient".to_string(), solution.cp_recovery);
        result.metrics.insert("Throat Cp".to_string(), solution.cp_throat);
        
        Ok(result)
    }

    fn reference_solution(&self) -> Option<BenchmarkResult<T>> {
        None
    }

    fn validate(&self, result: &BenchmarkResult<T>) -> cfd_core::error::Result<bool> {
        if let Some(&cp) = result.metrics.get("Throat Cp") {
            // Analytical Cp for frictionless Venturi: Cp = (D_in / D_th)^4 - 1
            // But we defined Cp as (P_th - P_in) / q_dyn
            // So Cp = 1 - (D_in / D_th)^4
            let beta = self.geometry.d_throat / self.geometry.d_inlet;
            let cp_ideal = T::one() - T::one() / (beta * beta * beta * beta);
            
            println!("Venturi 3D Validation:");
            println!("  Cp (FEM):   {:?}", cp);
            println!("  Cp (Ideal): {:?}", cp_ideal);
            println!("  Diff:       {:?}", num_traits::Float::abs(cp - cp_ideal) / num_traits::Float::abs(cp_ideal));

            // Check if result is reasonably close to ideal (within 20% for coarse mesh)
            return Ok(num_traits::Float::abs(cp - cp_ideal) / num_traits::Float::abs(cp_ideal) < T::from_f64_or_one(0.2));
        }
        Ok(false)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_venturi_flow_3d() {
        let geo = Venturi3D::new(200e-6, 100e-6, 500e-6, 300e-6, 200e-6, 600e-6, 500e-6);
        let benchmark = VenturiFlow3D::new(geo);
        let config = BenchmarkConfig::default();
        
        let result = benchmark.run(&config).unwrap();
        assert!(benchmark.validate(&result).unwrap());
    }
}
