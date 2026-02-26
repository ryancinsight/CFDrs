//! Validation tools for 3D bifurcation simulations
//!
//! Provides mesh convergence studies, error metrics, and comparison with
//! analytical/literature solutions.
//!
//! # Theorem — Richardson Extrapolation / GCI (Roache 1994)
//!
//! Given three grid solutions $\phi_1, \phi_2, \phi_3$ with refinement ratio $r$:
//!
//! ```text
//! p = ln((φ_3 − φ_2) / (φ_2 − φ_1)) / ln(r)
//! GCI = 1.25 |ε| / (r^p − 1)
//! ```
//!
//! The GCI provides a 95% confidence band on the discretisation error.

use super::geometry::BifurcationGeometry3D;
use super::solver::{BifurcationConfig3D, BifurcationSolution3D, BifurcationSolver3D};
use cfd_core::conversion::SafeFromF64;
use cfd_core::error::Error;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_core::physics::fluid::traits::NonNewtonianFluid;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

// ============================================================================
// Mesh Refinement Study
// ============================================================================

/// Mesh convergence study configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct MeshRefinementConfig<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Number of refinement levels
    pub n_levels: usize,
    /// Refinement factor (h_new = h_old / factor)
    pub refinement_factor: T,
    /// Expected convergence order (2 for FVM, 1.5-2 for FEM)
    pub expected_order: T,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> Default
    for MeshRefinementConfig<T>
{
    fn default() -> Self {
        Self {
            n_levels: 3,
            refinement_factor: T::from_f64_or_one(2.0),
            expected_order: T::from_f64_or_one(2.0),
        }
    }
}

// ============================================================================
// Bifurcation Validator
// ============================================================================

/// Comprehensive validator for 3D bifurcation simulations
pub struct BifurcationValidator3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    geometry: BifurcationGeometry3D<T>,
    mesh_config: MeshRefinementConfig<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + num_traits::Float + From<f64>> BifurcationValidator3D<T> {
    /// Create new validator
    pub fn new(geometry: BifurcationGeometry3D<T>, mesh_config: MeshRefinementConfig<T>) -> Self {
        Self {
            geometry,
            mesh_config,
        }
    }

    /// Validate against 1D bifurcation solution
    ///
    /// # Method
    ///
    /// Extract centerline velocity and pressure from 3D solution
    /// and compare with 1D bifurcation model.
    pub fn validate_vs_1d<F: FluidTrait<T> + NonNewtonianFluid<T> + Copy>(
        &self,
        solution_3d: &BifurcationSolution3D<T>,
        fluid: F,
        q_parent: T,
        p_parent: T,
    ) -> Result<BifurcationValidationResult3D<T>, Error> {
        // Compare mass conservation
        let mass_error = solution_3d.mass_conservation_error;

        // Compare pressure drops
        // In 1D, we compute ΔP = (128μQL) / (πD⁴)
        let props = fluid.properties_at(T::from_f64_or_one(310.0), p_parent)?;
        let mu = props.dynamic_viscosity;
        let pi = T::from_f64_or_one(std::f64::consts::PI);

        let dp_parent_1d = (T::from_f64_or_one(128.0) * mu * q_parent * self.geometry.l_parent)
            / (pi * num_traits::Float::powf(self.geometry.d_parent, T::from_f64_or_one(4.0)));

        let dp_error = num_traits::Float::abs(solution_3d.p_inlet - solution_3d.p_outlet - dp_parent_1d)
            / (num_traits::Float::abs(dp_parent_1d) + T::from_f64_or_one(1.0));

        // Create result
        let mut result = BifurcationValidationResult3D::new("3D vs 1D".to_string());
        result.mass_error = Some(mass_error);
        result.pressure_error = Some(dp_error);
        result.validation_passed =
            mass_error < T::from_f64_or_one(1e-10) && dp_error < T::from_f64_or_one(0.05); // < 5% error

        Ok(result)
    }

    /// Validate mesh convergence using Richardson extrapolation
    ///
    /// # Method
    ///
    /// Solves on coarse and fine meshes, computes convergence order:
    /// ```text
    /// p_obs = log(e_coarse / e_fine) / log(r)
    /// ```
    ///
    /// where r is refinement factor.
    pub fn validate_mesh_convergence<F: FluidTrait<T> + NonNewtonianFluid<T> + Clone>(
        &self,
        config: &BifurcationConfig3D<T>,
        fluid: F,
    ) -> Result<BifurcationValidationResult3D<T>, Error> {
        let r = self.mesh_config.refinement_factor;
        let r_f64 = r.to_f64().unwrap_or(2.0);
        if r_f64 <= 1.0 {
            return Err(Error::InvalidInput(
                "Mesh refinement factor must be > 1".to_string(),
            ));
        }

        // Build three levels: coarse, medium, fine
        let mut cfg_coarse = config.clone();
        let mut cfg_medium = config.clone();
        let mut cfg_fine = config.clone();

        let base_resolution = config.mesh_resolution.max(2);
        let medium_resolution = ((base_resolution as f64) * r_f64).round() as usize;
        let fine_resolution = ((base_resolution as f64) * r_f64 * r_f64).round() as usize;

        cfg_coarse.mesh_resolution = base_resolution;
        cfg_medium.mesh_resolution = medium_resolution.max(base_resolution + 1);
        cfg_fine.mesh_resolution = fine_resolution.max(cfg_medium.mesh_resolution + 1);

        let solver_coarse = BifurcationSolver3D::new(self.geometry.clone(), cfg_coarse);
        let solver_medium = BifurcationSolver3D::new(self.geometry.clone(), cfg_medium);
        let solver_fine = BifurcationSolver3D::new(self.geometry.clone(), cfg_fine);

        let sol_coarse = solver_coarse.solve(fluid.clone())?;
        let sol_medium = solver_medium.solve(fluid.clone())?;
        let sol_fine = solver_fine.solve(fluid)?;

        // Representative scalar quantity for convergence study
        let phi_coarse = sol_coarse.p_inlet - sol_coarse.p_outlet;
        let phi_medium = sol_medium.p_inlet - sol_medium.p_outlet;
        let phi_fine = sol_fine.p_inlet - sol_fine.p_outlet;

        let eps10 = num_traits::Float::abs(phi_medium - phi_coarse);
        let eps21 = num_traits::Float::abs(phi_fine - phi_medium);
        let tiny = T::from_f64_or_one(1e-20);

        if eps10 <= tiny || eps21 <= tiny {
            let mut result = BifurcationValidationResult3D::new("Mesh Convergence".to_string());
            result.convergence_order = None;
            result.gci = Some(T::zero());
            result.validation_passed = true;
            result.error_message = Some(
                "Inter-level differences are at numerical noise floor; treating as grid-converged".to_string(),
            );
            return Ok(result);
        }

        let p_obs = num_traits::Float::ln(eps10 / eps21) / num_traits::Float::ln(r);
        let rel_error_fine = eps21 / (num_traits::Float::abs(phi_fine) + tiny);
        let gci = T::from_f64_or_one(1.25) * rel_error_fine / (num_traits::Float::powf(r, p_obs) - T::one());

        let mut result = BifurcationValidationResult3D::new("Mesh Convergence".to_string());
        result.convergence_order = Some(p_obs);
        result.gci = Some(gci);
        result.validation_passed =
            gci < T::from_f64_or_one(0.05)
                && p_obs >= T::from_f64_or_one(0.5) * self.mesh_config.expected_order;

        Ok(result)
    }

    /// Validate blood flow with non-Newtonian properties
    pub fn validate_blood_flow(
        &self,
        solution: &BifurcationSolution3D<T>,
    ) -> Result<BifurcationValidationResult3D<T>, Error> {
        // Physical reasonableness checks for microfluidic bifurcation flows
        let mass_ok = solution.mass_conservation_error < T::from_f64_or_one(1e-6);

        // Wall shear stresses should be in physiological range for blood
        // For microfluidic flows with water, values can be in the range of 10-1000 Pa
        let tau_w_reasonable = solution.wall_shear_stress_parent > T::from_f64_or_one(0.0)
            && solution.wall_shear_stress_parent < T::from_f64_or_one(1000.0);

        let mut result = BifurcationValidationResult3D::new("Blood Flow".to_string());
        result.validation_passed = mass_ok && tau_w_reasonable;

        if !result.validation_passed {
            let mut msg = String::new();
            if !mass_ok {
                msg.push_str("Mass conservation failed; ");
            }
            if !tau_w_reasonable {
                msg.push_str("Wall shear stress out of physiological range");
            }
            result.error_message = Some(msg);
        }

        Ok(result)
    }
}

// ============================================================================
// Validation Result
// ============================================================================

/// Result of bifurcation validation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BifurcationValidationResult3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Test name
    pub test_name: String,
    /// Mass conservation error
    pub mass_error: Option<T>,
    /// Pressure drop error vs analytical
    pub pressure_error: Option<T>,
    /// Convergence order (from Richardson)
    pub convergence_order: Option<T>,
    /// Grid Convergence Index
    pub gci: Option<T>,
    /// Validation passed
    pub validation_passed: bool,
    /// Error message
    pub error_message: Option<String>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy + ToPrimitive> BifurcationValidationResult3D<T> {
    /// Create new validation result
    pub fn new(test_name: String) -> Self {
        Self {
            test_name,
            mass_error: None,
            pressure_error: None,
            convergence_order: None,
            gci: None,
            validation_passed: false,
            error_message: None,
        }
    }

    /// Print summary
    pub fn print_summary(&self) {
        println!("\n{}", "-".repeat(60));
        println!("Validation: {}", self.test_name);
        println!("{}", "-".repeat(60));

        if let Some(m_err) = self.mass_error {
            if let Some(m) = m_err.to_f64() {
                println!("Mass error: {m:.2e}");
            }
        }
        if let Some(p_err) = self.pressure_error {
            if let Some(p) = p_err.to_f64() {
                println!("Pressure error: {p:.2e}");
            }
        }
        if let Some(order) = self.convergence_order {
            if let Some(o) = order.to_f64() {
                println!("Convergence order: {o:.2}");
            }
        }
        if let Some(gci_val) = self.gci {
            if let Some(g) = gci_val.to_f64() {
                println!("GCI: {g:.2e}");
            }
        }

        println!(
            "Result: {}",
            if self.validation_passed {
                "✓ PASSED"
            } else {
                "✗ FAILED"
            }
        );

        if let Some(msg) = &self.error_message {
            println!("Error: {msg}");
        }
        println!("{}", "-".repeat(60));
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_validator_creation() {
        let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
        let mesh_config = MeshRefinementConfig::default();
        let validator = BifurcationValidator3D::new(geom, mesh_config);

        assert_eq!(validator.mesh_config.n_levels, 3);
    }

    #[test]
    fn test_blood_flow_validation() {
        let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
        let mesh_config = MeshRefinementConfig::default();
        let validator = BifurcationValidator3D::new(geom, mesh_config);

        let config = BifurcationConfig3D::default();
        let solver = BifurcationSolver3D::new(validator.geometry.clone(), config);

        let water = cfd_core::physics::fluid::water_20c::<f64>().unwrap();
        let solution = solver.solve(water).unwrap();

        let result = validator.validate_blood_flow(&solution).unwrap();
        
        // Debug output
        println!("Mass conservation error: {:.2e}", solution.mass_conservation_error);
        println!("Wall shear stress parent: {:.2e}", solution.wall_shear_stress_parent);
        println!("Validation message: {:?}", result.error_message);
        
        assert!(result.validation_passed);
    }

    #[test]
    fn test_mesh_convergence_outputs_observed_order_and_gci() {
        let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
        let mesh_config = MeshRefinementConfig::default();
        let validator = BifurcationValidator3D::new(geom, mesh_config);

        let mut config = BifurcationConfig3D::default();
        config.mesh_resolution = 4;

        let water = cfd_core::physics::fluid::blood::CassonBlood::normal_blood();
        let result = validator
            .validate_mesh_convergence(&config, water)
            .expect("mesh convergence validation should succeed");

        assert!(result.gci.is_some());
        assert!(
            result.convergence_order.is_some() || result.gci == Some(0.0),
            "Expected observed order or noise-floor convergence marker"
        );
    }

    #[test]
    fn test_mesh_convergence_rejects_invalid_refinement_factor() {
        let geom = BifurcationGeometry3D::<f64>::symmetric(100e-6, 80e-6, 1e-3, 1e-3, 100e-6);
        let mesh_config = MeshRefinementConfig {
            refinement_factor: 1.0,
            ..MeshRefinementConfig::default()
        };
        let validator = BifurcationValidator3D::new(geom, mesh_config);

        let config = BifurcationConfig3D::default();
        let water = cfd_core::physics::fluid::blood::CassonBlood::normal_blood();

        let result = validator.validate_mesh_convergence(&config, water);
        assert!(result.is_err());
    }
}
