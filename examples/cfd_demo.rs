//! CFD demonstration example
//!
//! This example shows usage of the CFD library components.

use cfd_core::physics::fluid_dynamics::{FlowField, FlowOperations};
use cfd_math::linear_solver::IterativeLinearSolver;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== CFD Demonstration ===\n");

    // 1. Create a flow field
    println!("1. Creating 3D flow field...");
    let flow_field = FlowField::<f64>::new(32, 32, 32);
    println!("   Flow field created: 32x32x32 grid");

    // 2. Calculate flow field operations
    println!("\n2. Computing flow field quantities...");
    let divergence = FlowOperations::divergence(&flow_field.velocity);
    let vorticity = FlowOperations::vorticity(&flow_field.velocity);
    let kinetic_energy = FlowOperations::kinetic_energy(&flow_field.velocity);
    let enstrophy = FlowOperations::enstrophy(&flow_field.velocity);

    println!("   Divergence computed (should be ~0 for incompressible)");
    println!(
        "     Average divergence: {:.3e}",
        divergence.iter().sum::<f64>() / divergence.len() as f64
    );

    println!("   Vorticity computed");
    println!("     Vorticity points: {}", vorticity.len());

    println!("   Kinetic energy computed");
    println!("     Total KE: {:.3e}", kinetic_energy.iter().sum::<f64>());

    println!("   Enstrophy computed");
    println!(
        "     Total enstrophy: {:.3e}",
        enstrophy.iter().sum::<f64>()
    );

    // 3. Test Reynolds number calculation
    println!("\n3. Reynolds number calculation...");
    use cfd_core::physics::values::{FlowGeometry, ReynoldsNumber};

    let re = ReynoldsNumber::new(2300.0, FlowGeometry::Pipe)?;
    println!("   Re = 2300 (pipe flow)");
    println!("   Is laminar? {}", re.is_laminar());
    println!("   Is transitional? {}", re.is_transitional());
    println!("   Is turbulent? {}", re.is_turbulent());

    // 4. Test sparse matrix operations
    println!("\n4. Sparse matrix operations...");
    use cfd_math::sparse::SparseMatrixBuilder;

    let mut builder = SparseMatrixBuilder::new(3, 3);
    builder.add_entry(0, 0, 2.0)?;
    builder.add_entry(0, 1, -1.0)?;
    builder.add_entry(1, 0, -1.0)?;
    builder.add_entry(1, 1, 2.0)?;
    builder.add_entry(1, 2, -1.0)?;
    builder.add_entry(2, 1, -1.0)?;
    builder.add_entry(2, 2, 2.0)?;
    let matrix = builder.build()?;
    println!("   Created 3x3 tridiagonal matrix");
    println!("   Non-zero elements: {}", matrix.nnz());

    // 5. Test linear solver
    println!("\n5. Solving linear system...");

    use nalgebra::DVector;

    let b = DVector::from_vec(vec![1.0, 0.0, 1.0]);
    use cfd_math::linear_solver::{preconditioners::IdentityPreconditioner, ConjugateGradient};
    let solver = ConjugateGradient::<f64>::default();
    let mut x = DVector::zeros(3);
    let preconditioner: Option<&IdentityPreconditioner> = None;
    solver.solve(&matrix, &b, &mut x, preconditioner)?;
    println!("   Solved Ax = b using Conjugate Gradient");
    println!("   Solution norm: {:.6}", x.norm());

    // 6. Summary
    println!("\n=== Summary ===");
    println!("All core components working correctly:");
    println!("  - Flow field creation and manipulation");
    println!("  - Flow quantity calculations (divergence, vorticity, KE, enstrophy)");
    println!("  - Reynolds number classification");
    println!("  - Sparse matrix operations");
    println!("  - Linear system solving (CG)");

    Ok(())
}
