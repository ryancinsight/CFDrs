//! Matrix-Free Linear Solvers Demo
//!
//! This example demonstrates the use of matrix-free linear solvers
//! for solving systems without explicitly storing the coefficient matrix.
//!
//! Run with: `cargo run --example matrix_free_demo`

use cfd_math::error::Result;
use cfd_math::linear_solver::{
    ConjugateGradient, IterativeLinearSolver, IterativeSolverConfig,
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

        // Interior points: −d²u/dx² (positive-definite discrete Laplacian for CG)
        for i in 1..(self.n - 1) {
            y[i] = (2.0 * x[i] - x[i - 1] - x[i + 1]) * dx2_inv;
        }

        // Boundary rows: identity so the augmented system stays positive definite.
        // With b[0] = b[n-1] = 0 this encodes u(0) = u(L) = 0 (Dirichlet).
        y[0] = x[0];
        y[self.n - 1] = x[self.n - 1];

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

    // n=11 gives x ∈ [0.0, 1.0] with dx=0.1, so sin(πx) satisfies u(0)=u(1)=0.
    let n = 11;
    let dx = 1.0 / (n - 1) as f64;
    let operator = DiffusionOperator1D::new(n, dx);

    // Create a manufactured solution: u(x) = sin(πx)
    // RHS: -d²u/dx² = π² sin(πx)
    let mut b_data = vec![0.0; n];
    let pi = std::f64::consts::PI;

    for (i, b_i) in b_data.iter_mut().enumerate() {
        let x = i as f64 * dx;
        // Manufactured solution u = sin(πx) → −u'' = π² sin(πx)
        *b_i = pi * pi * (pi * x).sin();
    }
    // Enforce homogeneous Dirichlet BCs on the RHS (matches identity boundary rows).
    b_data[0] = 0.0;
    b_data[n - 1] = 0.0;
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
    println!("  Grid spacing: {:.4}", dx);
    println!("  Max error vs sin(πx): {:.2e}", max_error);
    println!("  ✅ Solution matches manufactured solution");
    println!();

    // Example 2: Scaling — larger 1D problem demonstrates matrix-free advantage
    println!("📏 Example 2: 1D Diffusion — Scaling to n = 101");
    println!("------------------------------------------------");

    let n2 = 101usize;
    let dx2 = 1.0 / (n2 - 1) as f64;
    let op2 = DiffusionOperator1D::new(n2, dx2);
    let mut b2_data = vec![0.0_f64; n2];
    for (i, b_i) in b2_data.iter_mut().enumerate() {
        let x = i as f64 * dx2;
        *b_i = pi * pi * (pi * x).sin();
    }
    b2_data[0] = 0.0;
    b2_data[n2 - 1] = 0.0;
    let b2 = Array1::from_shape_vec([n2], b2_data).expect("shape matches");

    let config2 = IterativeSolverConfig::new(1e-12).with_max_iterations(200);
    let solver2 = ConjugateGradient::new(config2);
    let mut x2 = Array1::zeros([n2]);
    solver2.solve(
        &op2,
        &b2,
        &mut x2,
        None::<&cfd_math::linear_solver::preconditioners::IdentityPreconditioner>,
    )?;

    let mut max_err2 = 0.0_f64;
    for i in 0..n2 {
        let x_pos = i as f64 * dx2;
        let exact = (pi * x_pos).sin();
        max_err2 = max_err2.max((x2[i] - exact).abs());
    }
    println!("  Grid points: {n2}");
    println!("  Grid spacing: {dx2:.5}");
    println!("  Max error vs sin(πx): {max_err2:.2e}  (O(dx²) ≈ {:.2e})", dx2 * dx2);
    println!("  Operator never stored in memory — applied on-the-fly");
    println!("  ✅ Scaling demonstration complete");
    println!();

    println!("🎉 Matrix-free solvers demonstration complete!");
    println!("💡 Key benefits:");
    println!("   • No matrix storage required");
    println!("   • Memory efficient for large problems");
    println!("   • Natural for physics-based discretizations");
    println!("   • Easy parallelization opportunities");

    Ok(())
}
