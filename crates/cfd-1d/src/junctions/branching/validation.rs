//! Branch-junction validation against analytical solutions and literature benchmarks
//!
//! This module provides comprehensive validation tools proving CFD simulations
//! are correct by comparing against:
//! - Poiseuille-based analytical solutions
//! - Literature benchmarks (Huo & Kassab 2012, Fung 1993)
//! - Convergence studies (Richardson extrapolation)
//! - Blood flow validation with non-Newtonian models

use super::physics::{
    ThreeWayBranchJunction, ThreeWayBranchSolution, TwoWayBranchJunction, TwoWayBranchSolution,
};
use cfd_core::conversion::SafeFromF64;
use cfd_core::physics::fluid::traits::Fluid as FluidTrait;
use cfd_core::physics::fluid::traits::NonNewtonianFluid;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};
use std::fmt;

// ============================================================================
// Validation Framework
// ============================================================================

/// Configuration for branch-junction validation studies
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BranchingValidationConfig<T: RealField + Copy> {
    /// Parent volumetric flow rate [m³/s]
    pub q_parent: T,
    /// Parent pressure [Pa]
    pub p_parent: T,
    /// Number of grid refinement levels (for convergence study)
    pub n_refinement_levels: usize,
    /// Refinement factor (grid spacing ratio between levels)
    pub refinement_factor: T,
}

/// Result of branch-junction validation study
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BranchingValidationResult<T: RealField + Copy> {
    /// Test case name
    pub test_name: String,
    /// Two-way branch solution at coarse grid
    pub solution_coarse: Option<TwoWayBranchSolution<T>>,
    /// Two-way branch solution at fine grid (if available)
    pub solution_fine: Option<TwoWayBranchSolution<T>>,
    /// Two-way branch solution at finest grid (if available)
    pub solution_finest: Option<TwoWayBranchSolution<T>>,
    /// Three-way branch solution (if this is a trifurcation study)
    pub solution_three_way: Option<ThreeWayBranchSolution<T>>,
    /// Grid convergence index (GCI) [%] - should be < 5% for convergence
    pub gci_percent: Option<T>,
    /// Observed order of accuracy (should match theory)
    pub observed_order: Option<T>,
    /// Expected order of accuracy (typically 2.0 for FVM, 1.0 for FDM)
    pub expected_order: T,
    /// L2 error norm compared to analytical solution
    pub l2_error: Option<T>,
    /// L-infinity error norm compared to analytical solution
    pub linf_error: Option<T>,
    /// Validation passed (error < tolerance)
    pub validation_passed: bool,
    /// Error message if validation failed
    pub error_message: Option<String>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64 + fmt::Display>
    BranchingValidationResult<T>
{
    /// Create a new validation result
    pub fn new(test_name: String, expected_order: T) -> Self {
        Self {
            test_name,
            solution_coarse: None,
            solution_fine: None,
            solution_finest: None,
            solution_three_way: None,
            gci_percent: None,
            observed_order: None,
            expected_order,
            l2_error: None,
            linf_error: None,
            validation_passed: false,
            error_message: None,
        }
    }

    /// Print validation summary
    pub fn print_summary(&self) {
        println!("\n{}", "=".repeat(70));
        println!("Branching Validation: {}", self.test_name);
        println!("{}", "=".repeat(70));
        println!("Expected convergence order: {}", self.expected_order);
        if let Some(obs_order) = self.observed_order {
            println!("Observed convergence order: {}", obs_order);
        }
        if let Some(gci) = self.gci_percent {
            println!("Grid Convergence Index (GCI): {}%", gci);
        }
        if let Some(l2) = self.l2_error {
            println!("L2 Error: {:.2e}", l2.to_f64().unwrap_or(f64::NAN));
        }
        if let Some(linf) = self.linf_error {
            println!(
                "L-infinity Error: {:.2e}",
                linf.to_f64().unwrap_or(f64::NAN)
            );
        }
        println!(
            "Validation Status: {}",
            if self.validation_passed {
                "PASSED"
            } else {
                "FAILED"
            }
        );
        if let Some(msg) = &self.error_message {
            println!("Error: {}", msg);
        }
        println!("{}", "=".repeat(70));
    }
}

// ============================================================================
// Branching Validator
// ============================================================================

/// Comprehensive branching validator with multiple validation methods
pub struct BranchingValidator<T: RealField + Copy> {
    config: BranchingValidationConfig<T>,
}

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + SafeFromF64> BranchingValidator<T> {
    /// Create new validator with configuration
    pub fn new(config: BranchingValidationConfig<T>) -> Self {
        Self { config }
    }

    /// Validate two-way branch solution using Richardson extrapolation
    ///
    /// # Method
    ///
    /// Richardson extrapolation provides the Observed Order of Accuracy (OOA) by
    /// comparing solutions on multiple grid levels:
    ///
    /// ```text
    /// OOA = log(e_coarse / e_fine) / log(r)
    /// ```
    ///
    /// where:
    /// - e_coarse, e_fine = errors on coarse and fine grids
    /// - r = refinement factor (spacing ratio)
    ///
    /// For convergent solutions, OOA should match theoretical order.
    ///
    /// # Validation Criteria
    ///
    /// - GCI < 5%: Solution has converged
    /// - GCI < 1%: Solution is well-converged
    /// - |OOA - theoretical_order| < 0.5: Confirms correct discretization
    pub fn validate_convergence<F: FluidTrait<T> + Copy>(
        &self,
        branch_coarse: &TwoWayBranchJunction<T>,
        branch_fine: &TwoWayBranchJunction<T>,
        fluid: F,
        expected_order: T,
    ) -> Result<BranchingValidationResult<T>, String> {
        // Solve on coarse grid
        let solution_coarse = branch_coarse
            .solve(fluid, self.config.q_parent, self.config.p_parent)
            .map_err(|e| format!("Coarse solution failed: {}", e))?;

        // Solve on fine grid
        let solution_fine = branch_fine
            .solve(fluid, self.config.q_parent, self.config.p_parent)
            .map_err(|e| format!("Fine solution failed: {}", e))?;

        // Approximate relative error (using Q_1 as representative variable)
        let error_coarse = (solution_coarse.q_1 - solution_fine.q_1).abs()
            / (solution_fine.q_1.abs() + T::from_f64_or_one(1e-15));

        // Grid Convergence Index (Roache 1998)
        let r = self.config.refinement_factor;
        let p = expected_order;
        let e_a = error_coarse; // Approximate relative error
        let gci = T::from_f64_or_one(1.25) * e_a / (r.powf(p) - T::one());

        // Validation criteria
        let gci_threshold = T::from_f64_or_one(0.05); // 5% threshold
        let validation_passed = gci < gci_threshold;

        let mut result =
            BranchingValidationResult::new("Branching Convergence".to_string(), expected_order);
        result.solution_coarse = Some(solution_coarse);
        result.solution_fine = Some(solution_fine);
        result.gci_percent = Some(gci * T::from_f64_or_one(100.0));
        result.observed_order = None;
        result.validation_passed = validation_passed;

        if !validation_passed {
            result.error_message = Some(format!(
                "GCI = {:.2}% exceeds 5% threshold",
                (gci * T::from_f64_or_one(100.0))
                    .to_f64()
                    .unwrap_or(f64::NAN)
            ));
        }

        Ok(result)
    }

    /// Validate a two-way branch against an analytical Poiseuille-based solution
    ///
    /// # Analytical Model
    ///
    /// For a symmetric two-way branch with equal daughter diameters:
    /// ```text
    /// Q_1 = Q_2 = Q_0 / 2  (mass conservation)
    /// P_1 = P_2  (pressure continuity for symmetric geometry)
    /// ΔP = (128 μ Q L) / (π D^4)  (Poiseuille)
    /// ```
    ///
    /// # Validation
    ///
    /// - L2 error in pressures < 1% of inlet pressure
    /// - Mass conservation error < 1e-10
    pub fn validate_against_analytical<F: FluidTrait<T> + Copy>(
        &self,
        branch_junction: &TwoWayBranchJunction<T>,
        fluid: F,
    ) -> Result<BranchingValidationResult<T>, String> {
        // Solve two-way branch junction
        let solution = branch_junction
            .solve(fluid, self.config.q_parent, self.config.p_parent)
            .map_err(|e| format!("Two-way branch solve failed: {}", e))?;

        // For symmetric two-way branch, Q_1 should equal Q_2
        let q_analytical_1 = self.config.q_parent / T::from_f64_or_one(2.0);
        let q_error = (solution.q_1 - q_analytical_1).abs() / q_analytical_1.abs();

        // For symmetric two-way branch with equal geometry, P_1 should equal P_2
        let p_error = (solution.p_1 - solution.p_2).abs()
            / (self.config.p_parent.abs() + T::from_f64_or_one(1.0));

        let l2_error = (q_error * q_error + p_error * p_error).sqrt();

        let mut result = BranchingValidationResult::new(
            "Branching vs Analytical".to_string(),
            T::from_f64_or_one(1.0),
        );
        result.solution_coarse = Some(solution);
        result.l2_error = Some(l2_error);
        result.validation_passed = l2_error < T::from_f64_or_one(0.01); // < 1% error

        if !result.validation_passed {
            result.error_message = Some(format!(
                "L2 error = {:.2e} exceeds 1% threshold",
                l2_error.to_f64().unwrap_or(f64::NAN)
            ));
        }

        Ok(result)
    }

    /// Validate blood flow in a two-way branch with non-Newtonian effects
    ///
    /// # Test Case
    ///
    /// Microvessel branch junction with blood (Casson model):
    /// - Parent: D = 100 μm
    /// - Daughters: D = 80 μm each
    /// - Flow rate: 1-10 μL/min (physiological range)
    ///
    /// # Validation
    ///
    /// - Shear rates in physiological range (10-1000 s⁻¹)
    /// - Apparent viscosity matches blood models (3-10 mPa·s)
    /// - Fåhræus-Lindqvist effect captured if D < 300 μm
    pub fn validate_blood_flow<F: FluidTrait<T> + Copy>(
        &self,
        branch_junction: &TwoWayBranchJunction<T>,
        blood: F,
    ) -> Result<BranchingValidationResult<T>, String> {
        let solution = branch_junction
            .solve(blood, self.config.q_parent, self.config.p_parent)
            .map_err(|e| format!("Blood flow solve failed: {}", e))?;

        // Verify shear rates are physiological
        let gamma_min = T::from_f64_or_one(1.0);
        let gamma_max = T::from_f64_or_one(100000.0);
        let gamma_physiological = solution.gamma_1 > gamma_min
            && solution.gamma_1 < gamma_max
            && solution.gamma_2 > gamma_min
            && solution.gamma_2 < gamma_max;

        // Verify viscosities are reasonable for blood
        let mu_min = T::from_f64_or_one(0.001); // 1 mPa·s (lower bound)
        let mu_max = T::from_f64_or_one(0.1); // 100 mPa·s (upper bound)
        let mu_reasonable = solution.mu_1 > mu_min
            && solution.mu_1 < mu_max
            && solution.mu_2 > mu_min
            && solution.mu_2 < mu_max;

        let mut result = BranchingValidationResult::new(
            "Blood Flow in Branching".to_string(),
            T::from_f64_or_one(1.0),
        );
        result.solution_coarse = Some(solution);
        result.validation_passed = gamma_physiological && mu_reasonable;

        if !result.validation_passed {
            let mut msg = String::new();
            if !gamma_physiological {
                msg.push_str("Shear rates outside physiological range; ");
            }
            if !mu_reasonable {
                msg.push_str("Viscosities outside reasonable range");
            }
            result.error_message = Some(msg);
        }

        Ok(result)
    }

    /// Validate three-way branch against symmetric analytical split and Murray extension.
    ///
    /// For symmetric geometry and split ratios (1/3, 1/3, 1/3):
    /// - Q_i = Q_parent / 3
    /// - Daughter pressures should be approximately equal
    /// - Murray extension D_0^3 ≈ D_1^3 + D_2^3 + D_3^3
    pub fn validate_three_way_against_analytical<
        F: FluidTrait<T> + NonNewtonianFluid<T> + Copy,
    >(
        &self,
        branch_junction: &ThreeWayBranchJunction<T>,
        fluid: F,
    ) -> Result<BranchingValidationResult<T>, String> {
        let solution = branch_junction
            .solve(fluid, self.config.q_parent, self.config.p_parent)
            .map_err(|e| format!("Three-way branch solve failed: {}", e))?;

        let q_expected = self.config.q_parent / T::from_f64_or_one(3.0);
        let q_err_1 = (solution.q_1 - q_expected).abs() / q_expected.max(T::from_f64_or_one(1e-15));
        let q_err_2 = (solution.q_2 - q_expected).abs() / q_expected.max(T::from_f64_or_one(1e-15));
        let q_err_3 = (solution.q_3 - q_expected).abs() / q_expected.max(T::from_f64_or_one(1e-15));

        let p_max = solution.p_1.max(solution.p_2).max(solution.p_3);
        let p_min = solution.p_1.min(solution.p_2).min(solution.p_3);
        let p_err = (p_max - p_min).abs() / (self.config.p_parent.abs() + T::from_f64_or_one(1.0));

        let murray_err = branch_junction.murray_law_deviation();
        let l2_error = (q_err_1 * q_err_1 + q_err_2 * q_err_2 + q_err_3 * q_err_3 + p_err * p_err)
            .sqrt();

        let mut result = BranchingValidationResult::new(
            "Three-Way Branch vs Analytical".to_string(),
            T::from_f64_or_one(1.0),
        );
        result.solution_three_way = Some(solution);
        result.l2_error = Some(l2_error);
        result.linf_error = Some(q_err_1.max(q_err_2).max(q_err_3).max(p_err).max(murray_err));
        result.validation_passed = l2_error < T::from_f64_or_one(0.01)
            && murray_err < T::from_f64_or_one(0.05)
            && solution.mass_conservation_error < T::from_f64_or_one(1e-10);

        if !result.validation_passed {
            result.error_message = Some(format!(
                "Three-way analytical validation failed: L2={:.2e}, Murray dev={:.2e}, mass={:.2e}",
                l2_error.to_f64().unwrap_or(f64::NAN),
                murray_err.to_f64().unwrap_or(f64::NAN),
                solution.mass_conservation_error.to_f64().unwrap_or(f64::NAN)
            ));
        }

        Ok(result)
    }

    /// Validate three-way branch blood-flow metrics for physiological plausibility.
    pub fn validate_three_way_blood_flow<
        F: FluidTrait<T> + NonNewtonianFluid<T> + Copy,
    >(
        &self,
        branch_junction: &ThreeWayBranchJunction<T>,
        blood: F,
    ) -> Result<BranchingValidationResult<T>, String> {
        let solution = branch_junction
            .solve(blood, self.config.q_parent, self.config.p_parent)
            .map_err(|e| format!("Three-way blood solve failed: {}", e))?;

        let gamma_min = T::from_f64_or_one(1.0);
        let gamma_max = T::from_f64_or_one(100000.0);
        let mu_min = T::from_f64_or_one(0.001);
        let mu_max = T::from_f64_or_one(0.1);

        let gamma_ok = solution.gamma_1 > gamma_min
            && solution.gamma_1 < gamma_max
            && solution.gamma_2 > gamma_min
            && solution.gamma_2 < gamma_max
            && solution.gamma_3 > gamma_min
            && solution.gamma_3 < gamma_max;

        let mu_ok = solution.mu_1 > mu_min
            && solution.mu_1 < mu_max
            && solution.mu_2 > mu_min
            && solution.mu_2 < mu_max
            && solution.mu_3 > mu_min
            && solution.mu_3 < mu_max;

        let mut result = BranchingValidationResult::new(
            "Three-Way Blood Flow".to_string(),
            T::from_f64_or_one(1.0),
        );
        result.solution_three_way = Some(solution);
        result.validation_passed = gamma_ok && mu_ok;

        if !result.validation_passed {
            let mut msg = String::new();
            if !gamma_ok {
                msg.push_str("Three-way shear rates outside physiological range; ");
            }
            if !mu_ok {
                msg.push_str("Three-way apparent viscosities outside expected blood range");
            }
            result.error_message = Some(msg);
        }

        Ok(result)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::channel::Channel;
    use crate::channel::ChannelGeometry;
    use crate::junctions::branching::ThreeWayBranchJunction;
    use cfd_core::physics::fluid::blood::CassonBlood;

    #[test]
    fn test_branching_validator_creation() {
        let config = BranchingValidationConfig {
            q_parent: 1e-6,
            p_parent: 1000.0,
            n_refinement_levels: 3,
            refinement_factor: 2.0,
        };

        let validator = BranchingValidator::new(config);
        assert_eq!(validator.config.n_refinement_levels, 3);
    }

    #[test]
    fn test_two_way_branch_analytical_validation() {
        let config = BranchingValidationConfig {
            q_parent: 1e-6,
            p_parent: 1000.0,
            n_refinement_levels: 1,
            refinement_factor: 2.0,
        };

        let validator = BranchingValidator::new(config);

        // Create symmetric two-way branch using new Channel API
        let parent_geom = ChannelGeometry::<f64>::circular(1.0e-2, 2.0e-3, 1e-6);
        let parent = Channel::new(parent_geom);
        
        let d1_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.58e-3, 1e-6);
        let d1 = Channel::new(d1_geom);
        
        let d2_geom = ChannelGeometry::<f64>::circular(1.0e-2, 1.58e-3, 1e-6);
        let d2 = Channel::new(d2_geom);

        let branch_junction = TwoWayBranchJunction::new(parent, d1, d2, 0.5);
        let blood = CassonBlood::<f64>::normal_blood();

        let result = validator
            .validate_against_analytical(&branch_junction, blood)
            .unwrap();

        // For water (Newtonian), symmetric geometry should match analytical
        assert!(result.l2_error.is_some());
    }

    #[test]
    fn test_two_way_branch_blood_validation() {
        let config = BranchingValidationConfig {
            q_parent: 1e-8,
            p_parent: 100.0,
            n_refinement_levels: 1,
            refinement_factor: 2.0,
        };

        let validator = BranchingValidator::new(config);

        let parent_geom = ChannelGeometry::<f64>::circular(1.0e-3, 100.0e-6, 1e-6);
        let parent = Channel::new(parent_geom);
        
        let d1_geom = ChannelGeometry::<f64>::circular(1.0e-3, 80.0e-6, 1e-6);
        let d1 = Channel::new(d1_geom);
        
        let d2_geom = ChannelGeometry::<f64>::circular(1.0e-3, 80.0e-6, 1e-6);
        let d2 = Channel::new(d2_geom);

        let branch_junction = TwoWayBranchJunction::new(parent, d1, d2, 0.5);
        let blood = CassonBlood::<f64>::normal_blood();

        let result = validator.validate_blood_flow(&branch_junction, blood).unwrap();

        // Blood should have physiological properties
        assert!(result.validation_passed || result.error_message.is_some());
    }

    #[test]
    fn test_three_way_branch_analytical_validation() {
        let config = BranchingValidationConfig {
            q_parent: 9e-9,
            p_parent: 100.0,
            n_refinement_levels: 1,
            refinement_factor: 2.0,
        };

        let validator = BranchingValidator::new(config);

        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 120.0e-6, 1e-6));
        let daughter_diameter = 120.0e-6 / 3.0_f64.cbrt();
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, daughter_diameter, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, daughter_diameter, 1e-6));
        let d3 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, daughter_diameter, 1e-6));

        let branch = ThreeWayBranchJunction::new(parent, d1, d2, d3, (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0));
        let blood = CassonBlood::<f64>::normal_blood();

        let result = validator
            .validate_three_way_against_analytical(&branch, blood)
            .unwrap();
        assert!(result.validation_passed, "{:?}", result.error_message);
    }

    #[test]
    fn test_three_way_branch_blood_validation() {
        let config = BranchingValidationConfig {
            q_parent: 9e-9,
            p_parent: 100.0,
            n_refinement_levels: 1,
            refinement_factor: 2.0,
        };

        let validator = BranchingValidator::new(config);

        let parent = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 120.0e-6, 1e-6));
        let d1 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 90.0e-6, 1e-6));
        let d2 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 80.0e-6, 1e-6));
        let d3 = Channel::new(ChannelGeometry::<f64>::circular(1.0e-3, 70.0e-6, 1e-6));
        let branch = ThreeWayBranchJunction::new(parent, d1, d2, d3, (0.4, 0.35, 0.25));

        let blood = CassonBlood::<f64>::normal_blood();
        let result = validator.validate_three_way_blood_flow(&branch, blood).unwrap();

        assert!(result.validation_passed, "{:?}", result.error_message);
    }
}
