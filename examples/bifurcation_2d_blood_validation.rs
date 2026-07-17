//! 2D bifurcation flow validation with blood rheology.
//!
//! This example drives the public [`BifurcationSolver2D`] path for a symmetric
//! Casson case and an asymmetric Carreau-Yasuda case.  The reported flow rates
//! and mass-balance error are integrated from the solved staggered field; no
//! result is prescribed by the example.
//!
//! The validation limits match the established `cfd-2d` regression contract.
//! The Cartesian representation of angled branches contributes up to five
//! percent geometric outlet-measurement error at convergence.  The asymmetric
//! case uses the ten percent cross-fidelity bound because its narrower branch
//! occupies fewer cells at the validated resolution.

use std::error::Error;
use std::io;

use cfd_2d::solvers::bifurcation_flow::{BifurcationGeometry, BifurcationSolver2D};
use cfd_2d::solvers::ns_fvm::{BloodModel, SIMPLEConfig};
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};

/// Relative mass-balance limit for the converged symmetric grid.
const SYMMETRIC_MASS_BALANCE_LIMIT: f64 = 0.05;
/// Relative difference limit for equal geometric daughter branches.
const SYMMETRY_FLOW_LIMIT: f64 = 0.05;
/// Relative mass-balance limit for the asymmetric cross-fidelity grid.
const ASYMMETRIC_MASS_BALANCE_LIMIT: f64 = 0.10;
/// Blood density for both solved cases, in kg/m³.
const BLOOD_DENSITY: f64 = 1060.0;
/// Prescribed inlet velocity for both solved cases, in m/s.
///
/// This is the operating point exercised by the canonical `cfd-2d`
/// bifurcation regressions.
const INLET_VELOCITY: f64 = 0.1;

/// Solved flow quantities extracted from a bifurcation simulation.
#[derive(Debug)]
struct ValidationResult {
    name: &'static str,
    parent_flow: f64,
    daughter1_flow: f64,
    daughter2_flow: f64,
    mass_balance_error: f64,
}

impl ValidationResult {
    fn daughter_flow_difference(&self) -> f64 {
        (self.daughter1_flow - self.daughter2_flow).abs() / self.parent_flow
    }

    fn print(&self) {
        println!("\n{}", "=".repeat(70));
        println!("{}", self.name);
        println!("{}", "=".repeat(70));
        println!("Q_parent    = {:.6e} m²/s", self.parent_flow);
        println!("Q_daughter1 = {:.6e} m²/s", self.daughter1_flow);
        println!("Q_daughter2 = {:.6e} m²/s", self.daughter2_flow);
        println!(
            "Mass-balance error = {:.3}%",
            self.mass_balance_error * 100.0
        );
    }
}

fn solve(
    name: &'static str,
    geometry: BifurcationGeometry<f64>,
    blood: BloodModel<f64>,
    nx: usize,
    ny: usize,
    config: SIMPLEConfig<f64>,
) -> Result<ValidationResult, Box<dyn Error>> {
    let mut solver = BifurcationSolver2D::new(geometry, blood, BLOOD_DENSITY, nx, ny, config);
    let solution = solver.solve(INLET_VELOCITY)?;

    Ok(ValidationResult {
        name,
        parent_flow: solution.q_parent,
        daughter1_flow: solution.q_daughter1,
        daughter2_flow: solution.q_daughter2,
        mass_balance_error: solution.mass_balance_error,
    })
}

fn require(condition: bool, message: impl Into<String>) -> Result<(), Box<dyn Error>> {
    if condition {
        Ok(())
    } else {
        Err(io::Error::other(message.into()).into())
    }
}

fn validate_symmetric_casson() -> Result<ValidationResult, Box<dyn Error>> {
    let geometry = BifurcationGeometry::new_symmetric(1.0e-3, 2.0e-3, 0.7e-3, 2.0e-3, 0.5);
    let config = SIMPLEConfig {
        max_iterations: 5_000,
        tolerance: 1e-5,
        alpha_u: 0.5,
        alpha_p: 0.2,
        alpha_mu: 0.5,
        ..SIMPLEConfig::default()
    };
    let result = solve(
        "Case 1: symmetric Casson bifurcation",
        geometry,
        BloodModel::Casson(CassonBlood::normal_blood()),
        50,
        30,
        config,
    )?;

    require(
        result.parent_flow.is_finite() && result.parent_flow > 0.0,
        format!(
            "symmetric solve produced invalid parent flow {}",
            result.parent_flow
        ),
    )?;
    require(
        result.mass_balance_error <= SYMMETRIC_MASS_BALANCE_LIMIT,
        format!(
            "symmetric mass-balance error {} exceeds the {} limit",
            result.mass_balance_error, SYMMETRIC_MASS_BALANCE_LIMIT
        ),
    )?;
    require(
        result.daughter_flow_difference() <= SYMMETRY_FLOW_LIMIT,
        format!(
            "symmetric daughter-flow difference {} exceeds the {} limit",
            result.daughter_flow_difference(),
            SYMMETRY_FLOW_LIMIT
        ),
    )?;

    Ok(result)
}

fn validate_asymmetric_carreau_yasuda() -> Result<ValidationResult, Box<dyn Error>> {
    let geometry = BifurcationGeometry {
        parent_width: 2.0e-3,
        parent_length: 10.0e-3,
        daughter1_width: 1.5e-3,
        daughter1_length: 20.0e-3,
        daughter1_angle: 0.2,
        daughter2_width: 0.75e-3,
        daughter2_length: 20.0e-3,
        daughter2_angle: -0.2,
    };
    let config = SIMPLEConfig {
        max_iterations: 1_000,
        ..SIMPLEConfig::default()
    };
    let result = solve(
        "Case 2: asymmetric Carreau-Yasuda bifurcation",
        geometry,
        BloodModel::CarreauYasuda(CarreauYasudaBlood::normal_blood()),
        100,
        60,
        config,
    )?;

    require(
        result.parent_flow.is_finite() && result.parent_flow > 0.0,
        format!(
            "asymmetric solve produced invalid parent flow {}",
            result.parent_flow
        ),
    )?;
    require(
        result.mass_balance_error <= ASYMMETRIC_MASS_BALANCE_LIMIT,
        format!(
            "asymmetric mass-balance error {} exceeds the {} limit",
            result.mass_balance_error, ASYMMETRIC_MASS_BALANCE_LIMIT
        ),
    )?;
    require(
        result.daughter1_flow > result.daughter2_flow,
        format!(
            "wider daughter must receive more flow, got Q1={} and Q2={}",
            result.daughter1_flow, result.daughter2_flow
        ),
    )?;

    Ok(result)
}

fn main() -> Result<(), Box<dyn Error>> {
    println!("2D bifurcation blood-flow validation");

    let results = [
        validate_symmetric_casson()?,
        validate_asymmetric_carreau_yasuda()?,
    ];
    for result in &results {
        result.print();
    }

    println!("\nAll solved bifurcation validations passed.");
    Ok(())
}
