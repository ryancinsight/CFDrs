//! Matrix-Free Linear Solvers Demo
//!
//! This example demonstrates the use of matrix-free linear solvers
//! for solving systems without explicitly storing the coefficient matrix.
//!
//! Run with: `cargo run --example matrix_free_demo`

use cfd_math::error::Result;
use cfd_math::linear_solver::{
    ConjugateGradient, IterativeLinearSolver, IterativeSolverConfig, LaplacianOperator2D,
    LinearOperator,
};
use leto::Array1;

/// Simple 1D diffusion operator for demonstration
struct DiffusionOperator1D {
    n: usize,
    dx: f64,
}

impl DiffusionOperator1D {
    fn new(n: usize, dx: f64) -> Self {
        Self { n, dx }
    }
}

impl LinearOperator<f64> for DiffusionOperator1D {
    fn apply(&self, x: &Array1<f64>, y: &mut Array1<f64>) -> Result<()> {
        if x.shape()[0] != self.size() || y.shape()[0] != self.size() {
            return Err(cfd_core::error::Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let dx2_inv = 1.0 / (self.dx * self.dx);

        // Interior points: d²u/dx²
        for i in 1..(self.n - 1) {
            y[i] = (x[i - 1] + x[i + 1] - 2.0 * x[i]) * dx2_inv;
        }

        // Boundary conditions (homogeneous Dirichlet)
        y[0] = 0.0;
        y[self.n - 1] = 0.0;

        Ok(())
    }

    fn size(&self) -> usize {
        self.n
    }

    fn is_symmetric(&self) -> bool {
        true
    }
}

fn main() -> Result<()> {
    println!("🚀 Matrix-Free Linear Solvers Demo");
    println!("===================================");
    println!();

    // Example 1: 1D Diffusion Problem
    println!("📐 Example 1: 1D Diffusion Equation");
    println!("-----------------------------------");

    let n = 10;
    let dx = 0.1;
    let operator = DiffusionOperator1D::new(n, dx);

    // Create a manufactured solution: u(x) = sin(πx)
    // RHS: -d²u/dx² = π² sin(πx)
    let mut b_data = vec![0.0; n];
    let pi = std::f64::consts::PI;

    for (i, b_i) in b_data.iter_mut().enumerate() {
        let x = i as f64 * dx;
        *b_i = pi * pi * x.sin();
    }
    let b = Array1::from_shape_vec([n], b_data).expect("shape matches data length");

    // Solve using matrix-free CG
    let config = IterativeSolverConfig::new(1e-10).with_max_iterations(1000);
    let solver = ConjugateGradient::new(config);

    let mut x = Array1::zeros([n]);
    solver.solve(
        &operator,
        &b,
        &mut x,
        None::<&cfd_math::linear_solver::preconditioners::IdentityPreconditioner>,
    )?;

    // Check solution
    let mut max_error: f64 = 0.0;
    for i in 0..n {
        let x_pos = i as f64 * dx;
        let exact = (pi * x_pos).sin();
        let error = (x[i] - exact).abs();
        max_error = max_error.max(error);
    }

    println!("  Grid points: {}", n);
    println!("  Grid spacing: {:.3}", dx);
    println!("  Max error: {:.2e}", max_error);
    println!("  ✅ Solution computed successfully");
    println!();

    // Example 2: 2D Laplacian (CFD Pressure Solver)
    println!("🌊 Example 2: 2D Laplacian (CFD Pressure)");
    println!("----------------------------------------");

    let nx = 8;
    let ny = 8;
    let dx = 0.125;
    let dy = 0.125;

    let laplacian = LaplacianOperator2D::new(nx, ny, dx, dy);

    // Simple test: solve -∇²p = 1 with homogeneous Neumann BC
    let size = laplacian.size();
    let b_laplace = Array1::from_shape_fn([size], |_| 1.0);

    let mut p = Array1::zeros([size]);
    solver.solve(
        &laplacian,
        &b_laplace,
        &mut p,
        None::<&cfd_math::linear_solver::preconditioners::IdentityPreconditioner>,
    )?;

    println!("  Grid: {}x{}", nx, ny);
    println!("  Total DOF: {}", size);
    println!("  ✅ Pressure field computed");
    println!("  Sample pressure values:");
    for j in 0..3 {
        print!("    ");
        for i in 0..3 {
            let idx = j * nx + i;
            print!("{:.3} ", p[idx]);
        }
        println!();
    }
    println!();

    println!("🎉 Matrix-free solvers demonstration complete!");
    println!("💡 Key benefits:");
    println!("   • No matrix storage required");
    println!("   • Memory efficient for large problems");
    println!("   • Natural for physics-based discretizations");
    println!("   • Easy parallelization opportunities");

    Ok(())
}
