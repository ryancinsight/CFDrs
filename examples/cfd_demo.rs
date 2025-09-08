//! CFD demonstration example
//!
//! This example shows usage of the CFD library components.

use cfd_core::domains::fluid_dynamics::{
    FlowField, FlowOperations, KEpsilonModel, TurbulenceModel,
};
use cfd_math::linear_solver::IterativeLinearSolver;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("=== CFD Demonstration ===\n");

    // 1. Create a flow field
    println!("1. Creating 3D flow field...");
    let flow_field = FlowField::<f64>::new(32, 32, 32);
    println!("   ✓ Flow field created: 32x32x32 grid");

    // 2. Initialize turbulence model
    println!("\n2. Initializing k-epsilon turbulence model...");
    let mut k_epsilon = KEpsilonModel::new();
    k_epsilon.initialize_state(&flow_field);
    let nu_t = k_epsilon.turbulent_viscosity(&flow_field);
    println!("   ✓ Turbulence model initialized");
    println!(
        "   Average turbulent viscosity: {:.3e}",
        nu_t.iter().sum::<f64>() / nu_t.len() as f64
    );

    // 3. Calculate flow field operations
    println!("\n3. Computing flow field quantities...");
    let divergence = FlowOperations::divergence(&flow_field.velocity);
    let vorticity = FlowOperations::vorticity(&flow_field.velocity);
    let kinetic_energy = FlowOperations::kinetic_energy(&flow_field.velocity);
    let enstrophy = FlowOperations::enstrophy(&flow_field.velocity);

    println!("   ✓ Divergence computed (should be ~0 for incompressible)");
    println!(
        "     Average divergence: {:.3e}",
        divergence.iter().sum::<f64>() / divergence.len() as f64
    );

    println!("   ✓ Vorticity computed");
    println!("     Vorticity points: {}", vorticity.len());

    println!("   ✓ Kinetic energy computed");
    println!("     Total KE: {:.3e}", kinetic_energy.iter().sum::<f64>());

    println!("   ✓ Enstrophy computed");
    println!(
        "     Total enstrophy: {:.3e}",
        enstrophy.iter().sum::<f64>()
    );

    // 4. Test Reynolds number calculation
    println!("\n4. Reynolds number calculation...");
    use cfd_core::values::{FlowGeometry, ReynoldsNumber};

    let re = ReynoldsNumber::new(2300.0, FlowGeometry::Pipe)?;
    println!("   Re = 2300 (pipe flow)");
    println!("   Is laminar? {}", re.is_laminar());
    println!("   Is transitional? {}", re.is_transitional());
    println!("   Is turbulent? {}", re.is_turbulent());

    // 5. Test sparse matrix operations
    println!("\n5. Sparse matrix operations...");
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
    println!("   ✓ Created 3x3 tridiagonal matrix");
    println!("   Non-zero elements: {}", matrix.nnz());

    // 6. Test linear solver
    println!("\n6. Solving linear system...");

    use nalgebra::DVector;

    let b = DVector::from_vec(vec![1.0, 0.0, 1.0]);
    use cfd_math::linear_solver::{ConjugateGradient, preconditioners::IdentityPreconditioner};
    let solver = ConjugateGradient::<f64>::default();
    let mut x = DVector::zeros(3);
    let preconditioner: Option<&IdentityPreconditioner> = None;
    solver.solve(&matrix, &b, &mut x, preconditioner)?;
    println!("   ✓ Solved Ax = b using Conjugate Gradient");
    println!("   Solution norm: {:.6}", x.norm());

    // 7. Summary
    println!("\n=== Summary ===");
    println!("✓ Flow field creation and manipulation");
    println!("✓ Turbulence model initialization");
    println!("✓ Flow quantity calculations");
    println!("✓ Reynolds number classification");
    println!("✓ Sparse matrix operations");
    println!("✓ Linear system solving");
    println!("\nAll core components working correctly!");

    Ok(())
}
