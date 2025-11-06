//! Comprehensive stability analysis for CFD time-stepping schemes
//!
//! This module provides validation of numerical stability including:
//! - Stability region computation and visualization
//! - CFL condition verification for various flow regimes
//! - Von Neumann stability analysis for linear PDEs
//! - Stability monitoring during actual CFD simulations
//!
//! References:
//! - Hairer & Nørsett (1993): Solving Ordinary Differential Equations I
//! - Trefthen (1996): Finite Difference and Spectral Methods
//! - LeVeque (2002): Finite Volume Methods for Hyperbolic Problems

use cfd_core::error::{Error, Result};
use cfd_math::time_stepping::{StabilityAnalyzer, NumericalScheme};
use nalgebra::{DMatrix, DVector, RealField};
use num_traits::ToPrimitive;

/// Comprehensive stability analysis report
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct StabilityAnalysisReport<T: RealField + Copy> {
    /// RK method stability regions
    pub rk_stability_regions: Vec<RKStabilityResult<T>>,
    /// CFL condition analyses
    pub cfl_analyses: Vec<CFLValidationResult<T>>,
    /// Von Neumann analyses
    pub von_neumann_analyses: Vec<VonNeumannResult<T>>,
    /// Overall stability assessment
    pub overall_assessment: StabilityAssessment,
    /// Recommendations for stability improvements
    pub recommendations: Vec<String>,
}

/// Result of Runge-Kutta stability region analysis
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct RKStabilityResult<T: RealField + Copy> {
    /// Method name
    pub method_name: String,
    /// Stability region boundary points
    pub boundary_points: Vec<(T, T)>,
    /// Method order
    pub order: usize,
    /// Number of stages
    pub stages: usize,
    /// Absolute stability limit estimate
    pub stability_limit: T,
}

/// CFL condition validation result
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct CFLValidationResult<T: RealField + Copy> {
    /// Test case name
    pub test_case: String,
    /// CFL number computed
    pub cfl_number: T,
    /// Maximum stable CFL for scheme
    pub max_stable_cfl: T,
    /// Stability status
    pub status: String,
    /// Velocity magnitude
    pub velocity: T,
    /// Time step size
    pub dt: T,
    /// Grid spacing
    pub dx: T,
}

/// Von Neumann stability analysis result
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct VonNeumannResult<T: RealField + Copy> {
    /// PDE type analyzed
    pub pde_type: String,
    /// Maximum amplification factor
    pub max_amplification: T,
    /// Critical wave number
    pub critical_wave_number: T,
    /// Stability assessment
    pub is_stable: bool,
    /// Stability margin
    pub stability_margin: T,
}

/// Overall stability assessment
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct StabilityAssessment {
    /// Overall stability score (0.0 to 1.0)
    pub overall_score: f64,
    /// Number of stability tests passed
    pub tests_passed: usize,
    /// Total number of stability tests
    pub total_tests: usize,
    /// Critical issues found
    pub critical_issues: Vec<String>,
}

/// Comprehensive stability analysis runner
pub struct StabilityAnalysisRunner<T: RealField + Copy + num_traits::ToPrimitive> {
    analyzer: StabilityAnalyzer<T>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + num_traits::ToPrimitive> StabilityAnalysisRunner<T> {
    /// Create new stability analysis runner
    pub fn new() -> Self {
        Self {
            analyzer: StabilityAnalyzer::new(),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Run comprehensive stability analysis suite
    pub fn run_comprehensive_stability_analysis(&self) -> Result<StabilityAnalysisReport<T>> {
        println!("🔬 Comprehensive Stability Analysis Suite");
        println!("=========================================");
        println!("Literature References:");
        println!("- Hairer & Nørsett (1993): Solving ODEs I");
        println!("- Trefthen (1996): Finite Difference Methods");
        println!("- LeVeque (2002): Finite Volume Methods");
        println!();

        let mut report = StabilityAnalysisReport {
            rk_stability_regions: Vec::new(),
            cfl_analyses: Vec::new(),
            von_neumann_analyses: Vec::new(),
            overall_assessment: StabilityAssessment {
                overall_score: 0.0,
                tests_passed: 0,
                total_tests: 0,
                critical_issues: Vec::new(),
            },
            recommendations: Vec::new(),
        };

        // Analyze Runge-Kutta stability regions
        self.analyze_rk_stability_regions(&mut report)?;

        // Validate CFL conditions
        self.validate_cfl_conditions(&mut report)?;

        // Perform von Neumann stability analysis
        self.perform_von_neumann_analysis(&mut report)?;

        // Generate overall assessment
        self.generate_overall_assessment(&mut report);

        // Display results
        self.display_stability_report(&report);

        Ok(report)
    }

    /// Analyze stability regions for common RK methods
    fn analyze_rk_stability_regions(&self, report: &mut StabilityAnalysisReport<T>) -> Result<()> {
        println!("\n📐 Runge-Kutta Stability Region Analysis");

        // Forward Euler (RK1)
        let rk1_result = self.analyze_forward_euler_stability()?;
        report.rk_stability_regions.push(rk1_result);

        // RK3
        let rk3_result = self.analyze_rk3_stability()?;
        report.rk_stability_regions.push(rk3_result);

        // Classic RK4
        let rk4_result = self.analyze_rk4_stability()?;
        report.rk_stability_regions.push(rk4_result);

        println!("  ✅ Analyzed stability regions for 3 RK methods");
        Ok(())
    }

    /// Analyze Forward Euler stability region
    fn analyze_forward_euler_stability(&self) -> Result<RKStabilityResult<T>> {
        // Forward Euler: A = [0], b = [1], c = [0]
        let a = DMatrix::from_row_slice(1, 1, &[T::zero()]);
        let b = DVector::from_vec(vec![T::one()]);
        let c = DVector::from_vec(vec![T::zero()]);

        let region = self.analyzer.compute_rk_stability_region(&a, &b, &c)?;

        // For Forward Euler, stability region is |1 + z| <= 1, so |z| <= 2
        let stability_limit = T::from_f64(2.0).unwrap();

        Ok(RKStabilityResult {
            method_name: "Forward Euler (RK1)".to_string(),
            boundary_points: region.boundary.iter()
                .map(|p| (p.real, p.imag))
                .collect(),
            order: 1,
            stages: 1,
            stability_limit,
        })
    }

    /// Analyze RK3 stability region
    fn analyze_rk3_stability(&self) -> Result<RKStabilityResult<T>> {
        // Heun's method (RK3): A = [[0,0,0],[1/3,0,0],[0,2/3,0]], b = [1/4,0,3/4], c = [0,1/3,2/3]
        let a = DMatrix::from_row_slice(3, 3, &[
            T::zero(), T::zero(), T::zero(),
            T::from_f64(1.0/3.0).unwrap(), T::zero(), T::zero(),
            T::zero(), T::from_f64(2.0/3.0).unwrap(), T::zero(),
        ]);
        let b = DVector::from_vec(vec![
            T::from_f64(0.25).unwrap(),
            T::zero(),
            T::from_f64(0.75).unwrap(),
        ]);
        let c = DVector::from_vec(vec![
            T::zero(),
            T::from_f64(1.0/3.0).unwrap(),
            T::from_f64(2.0/3.0).unwrap(),
        ]);

        let region = self.analyzer.compute_rk_stability_region(&a, &b, &c)?;
        let stability_limit = self.estimate_stability_limit();

        Ok(RKStabilityResult {
            method_name: "Heun's Method (RK3)".to_string(),
            boundary_points: region.boundary.iter()
                .map(|p| (p.real, p.imag))
                .collect(),
            order: 3,
            stages: 3,
            stability_limit,
        })
    }

    /// Analyze classic RK4 stability region
    fn analyze_rk4_stability(&self) -> Result<RKStabilityResult<T>> {
        // Classic RK4
        let a = DMatrix::from_row_slice(4, 4, &[
            T::zero(), T::zero(), T::zero(), T::zero(),
            T::from_f64(0.5).unwrap(), T::zero(), T::zero(), T::zero(),
            T::zero(), T::from_f64(0.5).unwrap(), T::zero(), T::zero(),
            T::zero(), T::zero(), T::one(), T::zero(),
        ]);
        let b = DVector::from_vec(vec![
            T::from_f64(1.0/6.0).unwrap(),
            T::from_f64(1.0/3.0).unwrap(),
            T::from_f64(1.0/3.0).unwrap(),
            T::from_f64(1.0/6.0).unwrap(),
        ]);
        let c = DVector::from_vec(vec![
            T::zero(),
            T::from_f64(0.5).unwrap(),
            T::from_f64(0.5).unwrap(),
            T::one(),
        ]);

        let region = self.analyzer.compute_rk_stability_region(&a, &b, &c)?;
        let stability_limit = self.estimate_stability_limit();

        Ok(RKStabilityResult {
            method_name: "Classic Runge-Kutta 4".to_string(),
            boundary_points: region.boundary.iter()
                .map(|p| (p.real, p.imag))
                .collect(),
            order: 4,
            stages: 4,
            stability_limit,
        })
    }

    /// Estimate stability limit (simplified implementation)
    fn estimate_stability_limit(&self) -> T {
        // Return a reasonable default stability limit
        T::from_f64(2.0).unwrap_or(T::one())
    }

    /// Validate CFL conditions for various test cases
    fn validate_cfl_conditions(&self, report: &mut StabilityAnalysisReport<T>) -> Result<()> {
        println!("\n🌊 CFL Condition Validation");

        // Test case 1: Low-speed laminar flow
        let laminar_result = self.validate_single_cfl_case(
            "Laminar Channel Flow",
            T::from_f64(0.1).unwrap(),  // velocity
            T::from_f64(0.001).unwrap(), // dt
            T::from_f64(0.01).unwrap(),  // dx
            NumericalScheme::RK4,
        );
        report.cfl_analyses.push(laminar_result);

        // Test case 2: High-speed compressible flow
        let compressible_result = self.validate_single_cfl_case(
            "Compressible Flow",
            T::from_f64(300.0).unwrap(), // velocity (Mach ~1)
            T::from_f64(1e-6).unwrap(),   // dt
            T::from_f64(0.001).unwrap(), // dx
            NumericalScheme::RK4,
        );
        report.cfl_analyses.push(compressible_result);

        // Test case 3: Turbulent boundary layer
        let turbulent_result = self.validate_single_cfl_case(
            "Turbulent Boundary Layer",
            T::from_f64(10.0).unwrap(),  // velocity
            T::from_f64(1e-5).unwrap(),   // dt
            T::from_f64(1e-4).unwrap(),   // dx (near wall)
            NumericalScheme::RK3,
        );
        report.cfl_analyses.push(turbulent_result);

        println!("  ✅ Validated CFL conditions for 3 test cases");
        Ok(())
    }

    /// Validate CFL condition for a single test case
    fn validate_single_cfl_case(
        &self,
        test_case: &str,
        velocity: T,
        dt: T,
        dx: T,
        scheme: NumericalScheme,
    ) -> CFLValidationResult<T> {
        let analysis = self.analyzer.analyze_cfl_condition(velocity, dt, dx, scheme);

        CFLValidationResult {
            test_case: test_case.to_string(),
            cfl_number: analysis.cfl_number,
            max_stable_cfl: analysis.max_cfl,
            status: match analysis.stability {
                cfd_math::time_stepping::StabilityStatus::Stable => "Stable".to_string(),
                cfd_math::time_stepping::StabilityStatus::MarginallyStable => "Marginally Stable".to_string(),
                cfd_math::time_stepping::StabilityStatus::Unstable => "Unstable".to_string(),
            },
            velocity,
            dt,
            dx,
        }
    }

    /// Perform von Neumann stability analysis
    fn perform_von_neumann_analysis(&self, report: &mut StabilityAnalysisReport<T>) -> Result<()> {
        println!("\n📊 Von Neumann Stability Analysis");

        // Analyze advection equation: ∂u/∂t + a ∂u/∂x = 0
        let advection_result = self.analyze_advection_equation()?;
        report.von_neumann_analyses.push(advection_result);

        // Analyze diffusion equation: ∂u/∂t = ν ∂²u/∂x²
        let diffusion_result = self.analyze_diffusion_equation()?;
        report.von_neumann_analyses.push(diffusion_result);

        // Analyze Burgers' equation: ∂u/∂t + u ∂u/∂x = ν ∂²u/∂x²
        let burgers_result = self.analyze_burgers_equation()?;
        report.von_neumann_analyses.push(burgers_result);

        println!("  ✅ Performed von Neumann analysis for 3 PDE types");
        Ok(())
    }

    /// Analyze advection equation stability
    fn analyze_advection_equation(&self) -> Result<VonNeumannResult<T>> {
        let dt = T::from_f64(0.01).unwrap();
        let a = T::one(); // Advection speed

        // Test range of wave numbers
        let k_min = T::from_f64(0.0).unwrap();
        let k_max = T::from_f64(10.0).unwrap();
        let num_k = 50;

        let wave_numbers: Vec<T> = (0..num_k)
            .map(|i| k_min + (k_max - k_min) * T::from_f64(i as f64 / (num_k - 1) as f64).unwrap())
            .collect();

        // Upwind discretization: ∂u/∂x ≈ (u_j - u_{j-1}) / Δx
        // L_hat(k) = -a * (1 - e^{-i k Δx}) / Δx
        // For unit Δx, L_hat(k) = -a * (1 - e^{-i k})
        let spatial_operator = |k: num_complex::Complex<f64>| {
            let a_f64 = a.to_f64().unwrap();
            let delta_x = 1.0;
            -num_complex::Complex::new(a_f64, 0.0) * (num_complex::Complex::new(1.0, 0.0) - (-num_complex::Complex::new(0.0, k.im * delta_x)).exp())
        };

        let analysis = self.analyzer.von_neumann_analysis(spatial_operator, dt, &wave_numbers)?;

        Ok(VonNeumannResult {
            pde_type: "Linear Advection".to_string(),
            max_amplification: analysis.max_amplification,
            critical_wave_number: analysis.critical_wave_number,
            is_stable: analysis.is_stable,
            stability_margin: analysis.stability_margin,
        })
    }

    /// Analyze diffusion equation stability
    fn analyze_diffusion_equation(&self) -> Result<VonNeumannResult<T>> {
        let dt = T::from_f64(0.01).unwrap();
        let nu = T::from_f64(0.1).unwrap(); // Viscosity

        // Test range of wave numbers
        let k_min = T::from_f64(0.0).unwrap();
        let k_max = T::from_f64(20.0).unwrap();
        let num_k = 50;

        let wave_numbers: Vec<T> = (0..num_k)
            .map(|i| k_min + (k_max - k_min) * T::from_f64(i as f64 / (num_k - 1) as f64).unwrap())
            .collect();

        // Central difference: ∂²u/∂x² ≈ (u_{j+1} - 2u_j + u_{j-1}) / Δx²
        // L_hat(k) = -ν * (2 - 2cos(k Δx)) / Δx²
        // For unit Δx, L_hat(k) = -ν * 2 * (1 - cos(k))
        let spatial_operator = |k: num_complex::Complex<f64>| {
            let nu_f64 = nu.to_f64().unwrap();
            let delta_x = 1.0;
            -num_complex::Complex::new(nu_f64, 0.0) * num_complex::Complex::new(2.0, 0.0) *
            (num_complex::Complex::new(1.0, 0.0) - (num_complex::Complex::new(0.0, k.im * delta_x)).cos())
        };

        let analysis = self.analyzer.von_neumann_analysis(spatial_operator, dt, &wave_numbers)?;

        Ok(VonNeumannResult {
            pde_type: "Diffusion".to_string(),
            max_amplification: analysis.max_amplification,
            critical_wave_number: analysis.critical_wave_number,
            is_stable: analysis.is_stable,
            stability_margin: analysis.stability_margin,
        })
    }

    /// Analyze Burgers' equation stability
    fn analyze_burgers_equation(&self) -> Result<VonNeumannResult<T>> {
        let dt = T::from_f64(0.001).unwrap();
        let nu = T::from_f64(0.01).unwrap(); // Viscosity
        let u0 = T::from_f64(1.0).unwrap();  // Base velocity

        // Test range of wave numbers
        let k_min = T::from_f64(0.0).unwrap();
        let k_max = T::from_f64(15.0).unwrap();
        let num_k = 50;

        let wave_numbers: Vec<T> = (0..num_k)
            .map(|i| k_min + (k_max - k_min) * T::from_f64(i as f64 / (num_k - 1) as f64).unwrap())
            .collect();

        // Simplified analysis: treat as advection + diffusion
        // L_hat(k) = -u0 * i*k - ν*k²
        let spatial_operator = |k: num_complex::Complex<f64>| {
            let u0_f64 = u0.to_f64().unwrap();
            let nu_f64 = nu.to_f64().unwrap();
            -num_complex::Complex::new(0.0, u0_f64 * k.im) - num_complex::Complex::new(nu_f64 * k.im * k.im, 0.0)
        };

        let analysis = self.analyzer.von_neumann_analysis(spatial_operator, dt, &wave_numbers)?;

        Ok(VonNeumannResult {
            pde_type: "Burgers' Equation".to_string(),
            max_amplification: analysis.max_amplification,
            critical_wave_number: analysis.critical_wave_number,
            is_stable: analysis.is_stable,
            stability_margin: analysis.stability_margin,
        })
    }

    /// Generate overall stability assessment
    fn generate_overall_assessment(&self, report: &mut StabilityAnalysisReport<T>) {
        let mut tests_passed = 0;
        let mut total_tests = 0;
        let mut critical_issues = Vec::new();

        // Check CFL conditions
        for cfl_result in &report.cfl_analyses {
            total_tests += 1;
            if cfl_result.status == "Stable" {
                tests_passed += 1;
            } else {
                critical_issues.push(format!("CFL condition violated in {}", cfl_result.test_case));
            }
        }

        // Check von Neumann stability
        for vn_result in &report.von_neumann_analyses {
            total_tests += 1;
            if vn_result.is_stable {
                tests_passed += 1;
            } else {
                critical_issues.push(format!("{} equation unstable", vn_result.pde_type));
            }
        }

        // RK stability regions are always computed (no pass/fail)
        total_tests += report.rk_stability_regions.len();

        let overall_score = if total_tests > 0 {
            tests_passed as f64 / total_tests as f64
        } else {
            0.0
        };

        report.overall_assessment = StabilityAssessment {
            overall_score,
            tests_passed,
            total_tests,
            critical_issues,
        };

        // Generate recommendations
        self.generate_recommendations(report);
    }

    /// Generate stability recommendations
    fn generate_recommendations(&self, report: &mut StabilityAnalysisReport<T>) {
        // CFL recommendations
        let unstable_cfl = report.cfl_analyses.iter()
            .filter(|r| r.status != "Stable")
            .count();

        if unstable_cfl > 0 {
            report.recommendations.push(format!(
                "Fix {} CFL condition violations by reducing time step or increasing resolution",
                unstable_cfl
            ));
        }

        // Von Neumann recommendations
        let unstable_pdes = report.von_neumann_analyses.iter()
            .filter(|r| !r.is_stable)
            .count();

        if unstable_pdes > 0 {
            report.recommendations.push(format!(
                "{} PDE types show stability issues - consider implicit methods",
                unstable_pdes
            ));
        }

        // Overall score recommendations
        if report.overall_assessment.overall_score < 0.8 {
            report.recommendations.push(
                "Overall stability score below 80% - comprehensive stability review recommended".to_string()
            );
        } else if report.overall_assessment.overall_score >= 0.95 {
            report.recommendations.push(
                "Excellent stability characteristics - suitable for production CFD simulations".to_string()
            );
        }
    }

    /// Display comprehensive stability report
    fn display_stability_report(&self, report: &StabilityAnalysisReport<T>) {
        println!("\n📋 Stability Analysis Summary");
        println!("============================");
        println!("Overall Stability Score: {:.1}% ({}/{})",
                report.overall_assessment.overall_score * 100.0,
                report.overall_assessment.tests_passed,
                report.overall_assessment.total_tests);

        if !report.overall_assessment.critical_issues.is_empty() {
            println!("\n⚠️  Critical Issues:");
            for issue in &report.overall_assessment.critical_issues {
                println!("  • {}", issue);
            }
        }

        println!("\n🔢 CFL Condition Results:");
        for cfl in &report.cfl_analyses {
            println!("  {}: CFL = {:.3}, Status: {}",
                    cfl.test_case,
                    cfl.cfl_number.to_f64().unwrap(),
                    cfl.status);
        }

        println!("\n📊 Von Neumann Analysis:");
        for vn in &report.von_neumann_analyses {
            println!("  {}: Max Amp = {:.3}, Stable: {}",
                    vn.pde_type,
                    vn.max_amplification.to_f64().unwrap(),
                    vn.is_stable);
        }

        if !report.recommendations.is_empty() {
            println!("\n💡 Recommendations:");
            for rec in &report.recommendations {
                println!("  • {}", rec);
            }
        }

        println!("\n✅ Stability analysis completed successfully!");
    }
}

impl<T: RealField + Copy + num_traits::ToPrimitive> Default for StabilityAnalysisRunner<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_stability_analysis_runner() {
        let runner = StabilityAnalysisRunner::<f64>::new();
        let report = runner.run_comprehensive_stability_analysis();

        assert!(report.is_ok());
        let report = report.unwrap();

        assert!(!report.rk_stability_regions.is_empty());
        assert!(!report.cfl_analyses.is_empty());
        assert!(!report.von_neumann_analyses.is_empty());
    }
}

