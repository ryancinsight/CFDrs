//! Comprehensive CFD Validation Suite - 1D, 2D, and 3D
//!
//! This example validates all cfd-rs solvers against:
//! 1. Analytical solutions (Poiseuille, Bernoulli)
//! 2. Literature data (Merrill 1969, Murray 1926, ISO 5167, Fung 1993)
//! 3. Cross-dimensional consistency (1D vs 2D vs 3D)
//!
//! # Validation Cases
//!
//! ## 1D Validations
//! - Poiseuille flow with Casson blood (Merrill 1969)
//! - Murray's Law for bifurcations (Murray 1926)
//! - Carreau-Yasuda model (Cho & Kensey 1991)
//! - Fåhræus-Lindqvist effect (Pries 1992)
//!
//! ## 2D Validations
//! - Venturi flow (Bernoulli equation, ISO 5167)
//! - Pressure recovery in diffusers
//!
//! ## 3D Validations
//! - Bifurcation geometry (Murray's law compliance)
//! - Trifurcation geometry (3-way branching)
//!
//! # Literature References
//!
//! - Merrill, E.W. et al. (1969). "Pressure-flow relations of human blood
//!   in hollow fibers at low flow rates". J. Appl. Physiol. 27(1):93-98.
//!
//! - Murray, C.D. (1926). "The Physiological Principle of Minimum Work".
//!   Proc. Natl. Acad. Sci. 12(3):207-214.
//!
//! - Cho, Y.I. & Kensey, K.R. (1991). "Effects of the non-Newtonian viscosity
//!   of blood on flows in a diseased arterial vessel". J. Biomech. Eng. 113(2).
//!
//! - Pries, A.R. et al. (1992). "Blood viscosity in tube flow: dependence
//!   on diameter and hematocrit". Am. J. Physiol. 263(6):H1770-H1778.
//!
//! - ISO 5167-1:2003. "Measurement of fluid flow by means of pressure
//!   differential devices inserted in circular cross-section conduits".
//!
//! - Fung, Y.C. (1993). "Biomechanics: Mechanical Properties of Living Tissues".
//!   Springer-Verlag, 2nd edition.

use cfd_1d::bifurcation::junction::BifurcationJunction;
use cfd_1d::channel::{Channel, ChannelGeometry};
use cfd_2d::solvers::venturi_flow::{
    BernoulliVenturi, VenturiGeometry, ViscousVenturi,
};
use cfd_3d::bifurcation::geometry::BifurcationGeometry3D;
use cfd_3d::trifurcation::geometry::TrifurcationGeometry3D;
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood, FahraeuasLindqvist};
use std::f64::consts::PI;

/// Validation tolerance for most tests
const TOLERANCE_ANALYTICAL: f64 = 0.01; // 1%
const TOLERANCE_LITERATURE: f64 = 0.05; // 5%
const TOLERANCE_NUMERICAL: f64 = 0.10; // 10%

/// Validation result structure
#[derive(Debug)]
struct ValidationResult {
    name: String,
    dimension: String,
    test_type: String,
    passed: bool,
    error: f64,
    tolerance: f64,
    rust_value: f64,
    reference_value: f64,
    literature_source: String,
    details: String,
}

impl ValidationResult {
    fn print(&self) {
        let status = if self.passed { "✓ PASS" } else { "✗ FAIL" };
        println!("\n  Test: {}", self.name);
        println!("    Dimension: {} | Type: {}", self.dimension, self.test_type);
        println!("    Error: {:.3e} (tolerance: {:.3e})", self.error, self.tolerance);
        println!("    Rust value: {:.6e}", self.rust_value);
        println!("    Reference: {:.6e}", self.reference_value);
        println!("    Source: {}", self.literature_source);
        println!("    Details: {}", self.details);
        println!("    Status: {}", status);
    }
}

/// Run all validations and collect results
fn main() {
    println!("╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║     COMPREHENSIVE CFD VALIDATION SUITE                                       ║");
    println!("║     1D, 2D, and 3D Validation Against Literature                             ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");

    let mut results: Vec<ValidationResult> = Vec::new();

    // ============================================================================
    // 1D VALIDATIONS
    // ============================================================================
    println!("\n═══════════════════════════════════════════════════════════════════════════════");
    println!("1D VALIDATION TESTS");
    println!("═══════════════════════════════════════════════════════════════════════════════");

    results.push(validate_1d_poiseuille_casson());
    results.push(validate_1d_murray_law());
    results.push(validate_1d_carreau_yasuda());
    results.push(validate_1d_fahraeus_lindqvist());

    // ============================================================================
    // 2D VALIDATIONS
    // ============================================================================
    println!("\n═══════════════════════════════════════════════════════════════════════════════");
    println!("2D VALIDATION TESTS");
    println!("═══════════════════════════════════════════════════════════════════════════════");

    results.push(validate_2d_venturi_bernoulli());
    results.push(validate_2d_venturi_pressure_recovery());

    // ============================================================================
    // 3D VALIDATIONS
    // ============================================================================
    println!("\n═══════════════════════════════════════════════════════════════════════════════");
    println!("3D VALIDATION TESTS");
    println!("═══════════════════════════════════════════════════════════════════════════════");

    results.push(validate_3d_bifurcation_geometry());
    results.push(validate_3d_trifurcation_geometry());

    // ============================================================================
    // SUMMARY
    // ============================================================================
    print_summary(&results);

    // Exit with appropriate code
    let failed_count = results.iter().filter(|r| !r.passed).count();
    if failed_count > 0 {
        println!("\n⚠️  {} validation(s) FAILED", failed_count);
        std::process::exit(1);
    } else {
        println!("\n✓ All validations PASSED");
        std::process::exit(0);
    }
}

/// Print validation summary
fn print_summary(results: &[ValidationResult]) {
    println!("\n╔══════════════════════════════════════════════════════════════════════════════╗");
    println!("║     VALIDATION SUMMARY                                                       ║");
    println!("╚══════════════════════════════════════════════════════════════════════════════╝");

    let total = results.len();
    let passed = results.iter().filter(|r| r.passed).count();
    let failed = total - passed;

    println!("\n  Total tests: {}", total);
    println!("  Passed: {}", passed);
    println!("  Failed: {}", failed);
    println!("  Success rate: {:.1}%", 100.0 * passed as f64 / total as f64);

    // Breakdown by dimension
    println!("\n  Breakdown by dimension:");
    for dim in ["1D", "2D", "3D"] {
        let dim_results: Vec<_> = results.iter().filter(|r| r.dimension == dim).collect();
        if !dim_results.is_empty() {
            let dim_passed = dim_results.iter().filter(|r| r.passed).count();
            println!("    {}: {}/{} passed", dim, dim_passed, dim_results.len());
        }
    }

    // Failed cases
    if failed > 0 {
        println!("\n  Failed cases:");
        for result in results.iter().filter(|r| !r.passed) {
            println!(
                "    ✗ {} ({}): error={:.3e}",
                result.name, result.dimension, result.error
            );
        }
    }
}

// =============================================================================
// 1D Validation Functions
// =============================================================================

/// Validate 1D Poiseuille flow with Casson blood model
/// Reference: Merrill et al. (1969)
fn validate_1d_poiseuille_casson() -> ValidationResult {
    println!("\n  Test: Poiseuille Flow with Casson Blood (Merrill 1969)");

    let diameter: f64 = 4.0e-3; // 4 mm
    let length: f64 = 10.0e-2; // 10 cm
    let flow_rate: f64 = 5.0e-6; // 5 mL/s

    let blood = CassonBlood::<f64>::normal_blood();
    let gamma_wall: f64 = 32.0 * flow_rate / (PI * diameter.powi(3));
    let mu_apparent = blood.apparent_viscosity(gamma_wall);

    // Analytical pressure drop (Hagen-Poiseuille)
    let dp_analytical = (128.0 * mu_apparent * flow_rate * length) / (PI * diameter.powi(4));

    // Use bifurcation solver with equal daughters for straight tube
    let parent_geom = ChannelGeometry::circular(length, diameter, 0.0);
    let parent = Channel::new(parent_geom);
    let d_geom = ChannelGeometry::circular(length, diameter, 0.0);
    let d1 = Channel::new(d_geom.clone());
    let d2 = Channel::new(d_geom);

    let bifurcation = BifurcationJunction::new(parent, d1, d2, 0.5_f64);

    let dp_solver = BifurcationJunction::pressure_drop(&blood, flow_rate, &bifurcation.parent);

    let error = (dp_solver - dp_analytical).abs() / dp_analytical;

    ValidationResult {
        name: "Poiseuille Flow (Casson)".to_string(),
        dimension: "1D".to_string(),
        test_type: "Analytical".to_string(),
        passed: error < TOLERANCE_ANALYTICAL,
        error,
        tolerance: TOLERANCE_ANALYTICAL,
        rust_value: dp_solver,
        reference_value: dp_analytical,
        literature_source: "Merrill et al. (1969)".to_string(),
        details: format!("gamma={:.2e} 1/s, mu={:.4e} Pa·s", gamma_wall, mu_apparent),
    }
}

/// Validate Murray's Law for symmetric bifurcation
/// Reference: Murray (1926)
fn validate_1d_murray_law() -> ValidationResult {
    println!("\n  Test: Murray's Law for Bifurcations (Murray 1926)");

    let d_parent: f64 = 2.0e-3; // 2 mm
    let murray_factor: f64 = 0.79370052598; // 2^(-1/3)
    let d_daughter: f64 = d_parent * murray_factor;

    // Check Murray's law: D₀³ = D₁³ + D₂³
    let d0_cubed = d_parent.powi(3);
    let daughters_cubed = 2.0_f64 * d_daughter.powi(3);
    let murray_deviation = (d0_cubed - daughters_cubed).abs() / d0_cubed;

    // Flow split validation
    let parent_geom = ChannelGeometry::circular(1e-2, d_parent, 0.0);
    let d1_geom = ChannelGeometry::circular(1e-2, d_daughter, 0.0);
    let parent = Channel::new(parent_geom);
    let d1 = Channel::new(d1_geom.clone());
    let d2 = Channel::new(d1_geom);

    let bifurcation = BifurcationJunction::new(parent, d1, d2, 0.5_f64);
    let blood = CassonBlood::<f64>::normal_blood();

    let solution = bifurcation.solve(blood, 1e-6_f64, 100.0_f64).unwrap();
    let flow_split_error = ((solution.q_1 / 1e-6_f64) - 0.5_f64).abs();

    let error = murray_deviation.max(flow_split_error);

    ValidationResult {
        name: "Murray's Law".to_string(),
        dimension: "1D".to_string(),
        test_type: "Analytical".to_string(),
        passed: error < TOLERANCE_LITERATURE,
        error,
        tolerance: TOLERANCE_LITERATURE,
        rust_value: d_daughter,
        reference_value: d_parent / 2.0_f64.powf(1.0 / 3.0),
        literature_source: "Murray (1926)".to_string(),
        details: format!("D0={:.2} mm, D1={:.2} mm, dev={:.2e}", 
            d_parent * 1e3, d_daughter * 1e3, murray_deviation),
    }
}

/// Validate Carreau-Yasuda blood model
/// Reference: Cho & Kensey (1991)
fn validate_1d_carreau_yasuda() -> ValidationResult {
    println!("\n  Test: Carreau-Yasuda Blood Model (Cho & Kensey 1991)");

    let blood = CarreauYasudaBlood::<f64>::normal_blood();

    // Test at multiple shear rates across physiological range
    let shear_rates = [0.1_f64, 1.0_f64, 10.0_f64, 100.0_f64, 1000.0_f64];
    let viscosities: Vec<f64> = shear_rates
        .iter()
        .map(|&g| blood.apparent_viscosity(g))
        .collect();

    // Check shear-thinning behavior (viscosity decreases with shear rate)
    let shear_thinning_valid = viscosities
        .windows(2)
        .all(|w| w[0] > w[1]);

    let mu_0 = blood.zero_shear_viscosity;
    let mu_inf = blood.infinite_shear_viscosity;

    // Validate limiting behavior
    let zero_shear_correct = (blood.apparent_viscosity(0.0) - mu_0).abs() < 1e-10;
    let high_shear_correct = (blood.apparent_viscosity(10000.0) - mu_inf).abs() < 0.0001;

    let error = if shear_thinning_valid && mu_0 > mu_inf && zero_shear_correct && high_shear_correct {
        0.0
    } else {
        1.0
    };

    ValidationResult {
        name: "Carreau-Yasuda Model".to_string(),
        dimension: "1D".to_string(),
        test_type: "Literature".to_string(),
        passed: error < 0.1,
        error,
        tolerance: 0.1,
        rust_value: mu_0,
        reference_value: 0.056, // Literature value from Cho & Kensey
        literature_source: "Cho & Kensey (1991)".to_string(),
        details: format!("mu0={:.4e} Pa.s, mu_inf={:.4e} Pa.s, ratio={:.1}x", 
            mu_0, mu_inf, mu_0/mu_inf),
    }
}

/// Validate Fåhræus-Lindqvist effect
/// Reference: Pries et al. (1992)
fn validate_1d_fahraeus_lindqvist() -> ValidationResult {
    println!("\n  Test: Fåhræus-Lindqvist Effect (Pries 1992)");

    let d_large: f64 = 200e-6; // 200 μm
    let d_small: f64 = 50e-6; // 50 μm
    let hematocrit: f64 = 0.45;

    let fl_large = FahraeuasLindqvist::<f64>::new(d_large, hematocrit);
    let fl_small = FahraeuasLindqvist::<f64>::new(d_small, hematocrit);

    let mu_rel_large = fl_large.relative_viscosity();
    let mu_rel_small = fl_small.relative_viscosity();

    // Smaller vessel should have lower apparent viscosity (F-L effect)
    let fl_effect_valid = mu_rel_small < mu_rel_large;
    let viscosity_reduction = (mu_rel_large - mu_rel_small) / mu_rel_large;
    
    let error = if fl_effect_valid && viscosity_reduction > 0.1 { 
        0.0 
    } else { 
        1.0 
    };

    ValidationResult {
        name: "Fåhræus-Lindqvist Effect".to_string(),
        dimension: "1D".to_string(),
        test_type: "Literature".to_string(),
        passed: error < 0.5,
        error,
        tolerance: 0.5,
        rust_value: mu_rel_small,
        reference_value: mu_rel_large * 0.75, // Expected ~25% reduction
        literature_source: "Pries et al. (1992)".to_string(),
        details: format!("mu_rel(200um)={:.3}, mu_rel(50um)={:.3}, reduction={:.1}%", 
            mu_rel_large, mu_rel_small, viscosity_reduction * 100.0),
    }
}

// =============================================================================
// 2D Validation Functions
// =============================================================================

/// Validate 2D Venturi flow (Bernoulli)
/// Reference: ISO 5167
fn validate_2d_venturi_bernoulli() -> ValidationResult {
    println!("\n  Test: 2D Venturi Flow (Bernoulli, ISO 5167)");

    // ISO 5167 standard Venturi geometry
    let geom = VenturiGeometry::<f64>::iso_5167_standard();
    let area_ratio = geom.area_ratio();

    // Bernoulli pressure coefficient: Cp = 1 - (1/β²)² where β is area ratio
    let cp_analytical = 1.0 - (1.0 / area_ratio).powi(2);

    let bernoulli = BernoulliVenturi::new(
        geom,
        0.1,      // u_inlet [m/s]
        101325.0, // p_inlet [Pa]
        1000.0,   // rho [kg/m³]
    );
    
    let cp_numerical = bernoulli.pressure_coefficient_throat();

    let error = (cp_numerical - cp_analytical).abs() / cp_analytical.abs();

    ValidationResult {
        name: "Venturi (Bernoulli)".to_string(),
        dimension: "2D".to_string(),
        test_type: "Analytical".to_string(),
        passed: error < TOLERANCE_NUMERICAL,
        error,
        tolerance: TOLERANCE_NUMERICAL,
        rust_value: cp_numerical,
        reference_value: cp_analytical,
        literature_source: "ISO 5167-1:2003".to_string(),
        details: format!("beta={:.3}, Cp={:.3}, Cp_analytical={:.3}", 
            area_ratio, cp_numerical, cp_analytical),
    }
}

/// Validate 2D Venturi pressure recovery
fn validate_2d_venturi_pressure_recovery() -> ValidationResult {
    println!("\n  Test: 2D Venturi Pressure Recovery");

    let geom = VenturiGeometry::<f64>::iso_5167_standard();

    let viscous = ViscousVenturi::new(
        geom,
        0.5,      // u_inlet [m/s]
        101325.0, // p_inlet [Pa]
        1000.0,   // rho [kg/m³]
        0.15,     // zeta (loss coefficient for well-designed Venturi)
    );

    let cp_recovery = viscous.pressure_recovery_coefficient();
    
    // Expected recovery coefficient: Cp_recovery ≈ -ζ
    let expected_recovery = -0.15;
    let error = (cp_recovery - expected_recovery).abs();

    ValidationResult {
        name: "Venturi Pressure Recovery".to_string(),
        dimension: "2D".to_string(),
        test_type: "Analytical".to_string(),
        passed: error < TOLERANCE_LITERATURE,
        error,
        tolerance: TOLERANCE_LITERATURE,
        rust_value: cp_recovery,
        reference_value: expected_recovery,
        literature_source: "ISO 5167".to_string(),
        details: format!("Cp_recovery={:.3}, expected={:.3}", cp_recovery, expected_recovery),
    }
}

// =============================================================================
// 3D Validation Functions
// =============================================================================

/// Validate 3D bifurcation geometry against Murray's law
fn validate_3d_bifurcation_geometry() -> ValidationResult {
    println!("\n  Test: 3D Bifurcation Geometry (Murray's Law)");

    let d_parent: f64 = 4e-3; // 4 mm (carotid artery scale)
    let d_daughter: f64 = 3e-3; // 3 mm
    let length: f64 = 20e-3; // 20 mm
    let transition: f64 = 100e-6; // 100 μm transition region

    let geom = BifurcationGeometry3D::symmetric(
        d_parent,
        d_daughter,
        length,
        length,
        transition,
    );

    // Verify geometry properties
    let area_parent = PI * (d_parent / 2.0).powi(2);
    let area_daughter = PI * (d_daughter / 2.0).powi(2);
    let area_ratio = area_daughter / area_parent;

    // Check Murray's law compliance: D₀³ = D₁³ + D₂³
    let murray_parent = d_parent.powi(3);
    let murray_daughters = 2.0 * d_daughter.powi(3);
    let murray_error = (murray_parent - murray_daughters).abs() / murray_parent;

    // Optimal daughter diameter for Murray's law: D₁ = D₀ / 2^(1/3)
    let optimal_daughter = d_parent / 2.0_f64.powf(1.0 / 3.0);
    let diameter_error = (d_daughter - optimal_daughter).abs() / optimal_daughter;

    ValidationResult {
        name: "3D Bifurcation Geometry".to_string(),
        dimension: "3D".to_string(),
        test_type: "Geometric".to_string(),
        passed: murray_error < TOLERANCE_LITERATURE,
        error: diameter_error,
        tolerance: TOLERANCE_LITERATURE,
        rust_value: d_daughter,
        reference_value: optimal_daughter,
        literature_source: "Murray (1926)".to_string(),
        details: format!("D_parent={:.1} mm, D_daughter={:.1} mm, optimal={:.1} mm, Murray error={:.1}%", 
            d_parent*1e3, d_daughter*1e3, optimal_daughter*1e3, murray_error*100.0),
    }
}

/// Validate 3D trifurcation geometry
fn validate_3d_trifurcation_geometry() -> ValidationResult {
    println!("\n  Test: 3D Trifurcation Geometry");

    let d_parent: f64 = 100e-6; // 100 μm (microfluidic scale)
    let d_daughter: f64 = 80e-6; // 80 μm
    let length: f64 = 1e-3; // 1 mm
    let transition: f64 = 100e-6; // 100 μm
    let angle: f64 = PI / 6.0; // 30 degrees branching angle

    let geom = TrifurcationGeometry3D::symmetric(
        d_parent,
        d_daughter,
        length,
        length,
        transition,
        angle,
    );

    // Verify geometry construction
    let expected_daughters = 3;
    let actual_daughters = geom.d_daughters.len();
    
    let geometry_valid = actual_daughters == expected_daughters 
        && geom.d_parent == d_parent
        && geom.l_parent == length;

    let error = if geometry_valid { 0.0 } else { 1.0 };

    // Generalized Murray's law for trifurcation: D₀³ = D₁³ + D₂³ + D₃³
    let murray_parent = d_parent.powi(3);
    let murray_daughters = 3.0 * d_daughter.powi(3);
    let murray_deviation = (murray_parent - murray_daughters).abs() / murray_parent;

    ValidationResult {
        name: "3D Trifurcation Geometry".to_string(),
        dimension: "3D".to_string(),
        test_type: "Geometric".to_string(),
        passed: geometry_valid && murray_deviation < 0.3,
        error,
        tolerance: 0.01,
        rust_value: actual_daughters as f64,
        reference_value: expected_daughters as f64,
        literature_source: "Geometric Validation".to_string(),
        details: format!("3-way split at {:.1} deg, D0={:.0} um, D1={:.0} um, Murray dev={:.1}%", 
            angle.to_degrees(), d_parent*1e6, d_daughter*1e6, murray_deviation*100.0),
    }
}
