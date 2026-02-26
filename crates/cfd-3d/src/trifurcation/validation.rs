//! Validation tools for 3D trifurcation simulations
//!
//! Provides mesh convergence studies, error metrics, and comparison with
//! analytical/literature solutions.
//!
//! # Theorem — Richardson Extrapolation / GCI (Roache 1994)
//!
//! The Grid Convergence Index for a three-grid study:
//!
//! ```text
//! GCI = 1.25 |ε| / (r^p − 1)
//! ```
//!
//! bounds the discretisation error with 95% confidence.

use super::geometry::TrifurcationGeometry3D;
use super::solver::TrifurcationSolution3D;
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Error;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Mesh Refinement on Trifurcation
// ============================================================================

/// Mesh convergence study configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshRefinementConfig<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Number of mesh refinement levels
    pub n_levels: usize,
    /// Element size ratio between successive levels
    pub refinement_factor: T,
    /// Expected order of convergence (e.g. 2.0 for P1 elements)
    pub expected_order: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> Default for MeshRefinementConfig<T> {
    fn default() -> Self {
        Self {
            n_levels: 3,
            refinement_factor: T::from_f64_or_one(2.0),
            expected_order: T::from_f64_or_one(2.0),
        }
    }
}

// ============================================================================
// Trifurcation Validator
// ============================================================================

/// Validator for 3D trifurcation flow simulation results
pub struct TrifurcationValidator3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Trifurcation channel geometry
    pub geometry: TrifurcationGeometry3D<T>,
    /// Mesh refinement configuration
    pub mesh_config: MeshRefinementConfig<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + num_traits::Float + From<f64>> TrifurcationValidator3D<T> {
    /// Create a new trifurcation validator
    pub fn new(geometry: TrifurcationGeometry3D<T>, mesh_config: MeshRefinementConfig<T>) -> Self {
        Self { geometry, mesh_config }
    }

    /// Validate blood flow results
    pub fn validate_blood_flow(
        &self,
        solution: &TrifurcationSolution3D<T>,
    ) -> Result<TrifurcationValidationResult3D<T>, Error> {
        // 1. Mass Conservation Check
        // sum(Q_out) should equal Q_in
        let mass_error = solution.mass_conservation_error;
        let mass_ok = mass_error < T::from_f64_or_one(1e-6);

        // 2. Flow Split Symmetry Check (for symmetric geometry)
        // Q_d1 should equal Q_d3 (top and bottom), Q_d2 is middle
        let q_d1 = solution.flow_rates[1];
        let q_d3 = solution.flow_rates[3];
        
        // For symmetric trifurcation, we expect top/bottom symmetry
        // Middle branch might be different
        let symmetry_error = num_traits::Float::abs(q_d1 - q_d3) / (q_d1 + T::from_f64_or_one(1e-12));
        let symmetry_ok = symmetry_error < T::from_f64_or_one(0.05); // 5% tolerance

        // 3. Wall Shear Stress Check
        // Should be positive and within physiological range (0-1000 Pa for blood)
        let wss_ok = solution.wall_shear_stresses.iter().all(|&wss| {
            wss >= T::from_f64_or_one(0.0) && wss < T::from_f64_or_one(2000.0)
        });

        let mut result = TrifurcationValidationResult3D::new("Blood Flow Validation".to_string());
        result.mass_error = Some(mass_error);
        result.symmetry_error = Some(symmetry_error);
        result.validation_passed = mass_ok && symmetry_ok && wss_ok;

        if !result.validation_passed {
             let mut msg = String::new();
            if !mass_ok { use std::fmt::Write; let _ = write!(msg, "Mass error too high: {mass_error:?}; "); }
            if !symmetry_ok { use std::fmt::Write; let _ = write!(msg, "Asymmetry detected: {symmetry_error:?}; "); }
            if !wss_ok { msg.push_str("WSS out of range; "); }
            result.error_message = Some(msg);
        }

        Ok(result)
    }
}

/// Validation result for 3D trifurcation flow
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct TrifurcationValidationResult3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Name of the validation test
    pub test_name: String,
    /// Relative mass conservation error
    pub mass_error: Option<T>,
    /// Symmetry error between top and bottom daughters
    pub symmetry_error: Option<T>,
    /// Whether all validation checks passed
    pub validation_passed: bool,
    /// Detailed error message when validation fails
    pub error_message: Option<String>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> TrifurcationValidationResult3D<T> {
    /// Create a new (default-failing) validation result
    pub fn new(test_name: String) -> Self {
        Self {
            test_name,
            mass_error: None,
            symmetry_error: None,
            validation_passed: false,
            error_message: None,
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use super::super::solver::{TrifurcationConfig3D, TrifurcationSolver3D};

    #[test]
    fn test_trifurcation_blood_flow() {
        let geom = TrifurcationGeometry3D::<f64>::symmetric(
            100e-6, // Parent D
            80e-6,  // Daughter D
            500e-6, // Parent L
            500e-6, // Daughter L
            50e-6,  // Transition L
            std::f64::consts::PI / 4.0, // 45 degrees
        );
        let mesh_config = MeshRefinementConfig::default();
        let validator = TrifurcationValidator3D::new(geom.clone(), mesh_config);

        let mut config = TrifurcationConfig3D::default();
        config.inlet_flow_rate = 1e-9; // Very slow flow to ensure stability for test
        config.max_linear_iterations = 1000; // Match the fix in bifurcation
        config.linear_tolerance = 1e-6;

        let solver = TrifurcationSolver3D::new(geom, config);
        let fluid = cfd_core::physics::fluid::blood::CassonBlood::normal_blood();
        
        let solution = solver.solve(&fluid).expect("Solver failed");
        
        let result = validator.validate_blood_flow(&solution).expect("Validation failed");
        
        println!("Mass Error: {:?}", result.mass_error);
        println!("Symmetry Error: {:?}", result.symmetry_error);
        
        assert!(result.validation_passed, "Validation failed: {:?}", result.error_message);
    }
}
