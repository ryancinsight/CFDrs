//! 2D Poiseuille flow validation with Casson blood rheology
//!
//! Solves steady non-Newtonian Poiseuille flow in a 2D channel with Casson
//! blood model and prints the velocity profile, comparing against the
//! Newtonian analytical solution for the same pressure gradient.

use cfd_2d::solvers::{BloodModel, PoiseuilleConfig, PoiseuilleFlow2D};
use cfd_core::physics::fluid::blood::CassonBlood;

fn main() {
    println!("CFD Validation: 2D Poiseuille Flow with Casson Blood");
    println!("====================================================\n");

    // Channel geometry: 200 µm height, 1 mm length
    let mut config = PoiseuilleConfig::<f64>::default();
    config.height = 200e-6;
    config.width = 200e-6;
    config.length = 1e-3;
    config.ny = 64;
    config.pressure_gradient = 500.0; // Pa/m
    config.tolerance = 1e-8;
    config.max_iterations = 200;

    let blood = BloodModel::Casson(CassonBlood::normal_blood());

    let mut solver = PoiseuilleFlow2D::new(config, blood);
    let iterations = solver.solve().expect("Poiseuille solve failed");
    println!("Converged in {iterations} iterations");

    // Newtonian reference viscosity for analytical comparison (blood ~3.5 mPa.s)
    let mu_newtonian: f64 = 3.5e-3;
    let analytical_profile = solver.analytical_solution(mu_newtonian);
    let u_max_analytical = analytical_profile.iter().cloned().fold(0.0_f64, f64::max);
    let u_max_numerical = solver.max_velocity();

    println!("Max velocity (numerical) : {u_max_numerical:.6e} m/s");
    println!("Max velocity (Newtonian) : {u_max_analytical:.6e} m/s");

    // Non-Newtonian Casson model produces a blunted profile compared to the
    // Newtonian parabola, so a non-zero L2 difference is expected.
    let profile = solver.velocity_profile();
    let mut err_sq = 0.0;
    let mut norm_sq = 0.0;
    for (i, &u_num) in profile.iter().enumerate() {
        let u_ana = analytical_profile[i];
        err_sq += (u_num - u_ana).powi(2);
        norm_sq += u_ana.powi(2);
    }
    let l2_error = (err_sq / norm_sq).sqrt();

    println!("Relative L2 error vs Newtonian analytical: {l2_error:.4e}");
}
