//! Demonstration of MUSCL (Monotonic Upstream-Centered Scheme for Conservation Laws) schemes
//!
//! This example shows how to use higher-order MUSCL schemes for improved spatial discretization
//! in CFD simulations, providing 2nd/3rd order accuracy while maintaining monotonicity.

use cfd_2d::fields::SimulationFields;
use cfd_2d::grid::StructuredGrid2D;
use cfd_2d::physics::momentum::{
    DiscretizationScheme, MusclDiscretization, MusclOrder, MusclScheme, Superbee, VanLeer,
};
use nalgebra::Vector2;

fn main() -> Result<(), Box<dyn std::error::Error>> {
    println!("MUSCL Schemes Demonstration");
    println!("===========================");

    // Create a simple 2D grid (10x10 cells)
    let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)?;
    println!("Grid: {}x{} cells", grid.nx, grid.ny);

    // Initialize simulation fields
    let mut fields = SimulationFields::new(10, 10);

    // Set up test velocity field (simple shear flow)
    for i in 0..10 {
        for j in 0..10 {
            let y_pos = j as f64 * 0.1; // Normalized y position
                                        // Linear velocity profile: u = y, v = 0
            fields.set_velocity_at(i, j, &Vector2::new(y_pos, 0.0));
            fields.p[(i, j)] = 0.0;
        }
    }

    println!("Initial conditions: Linear shear flow (u = y, v = 0)");

    // Demonstrate MUSCL2 with Superbee limiter
    demonstrate_muscl2_superbee(&fields)?;

    // Demonstrate MUSCL2 with van Leer limiter
    demonstrate_muscl2_vanleer(&fields)?;

    println!("\nDemonstration completed successfully!");
    println!("MUSCL schemes provide higher-order accuracy while preserving monotonicity.");
    println!("- MUSCL2: 2nd order accurate with TVD limiting");
    println!("- MUSCL3: 3rd order accurate (QUICK-like) with limiter blending");
    println!("- TVD limiters ensure monotonicity preservation at discontinuities");

    Ok(())
}

fn demonstrate_muscl2_superbee(
    fields: &SimulationFields<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n--- MUSCL2 with Superbee Limiter ---");

    let limiter = Superbee;
    let muscl = MusclScheme::new(limiter, MusclOrder::SecondOrder);
    let scheme = MusclDiscretization::new(muscl);

    println!("Testing convective flux calculation:");
    println!("Position: i=5, j=5 (interior cell)");
    println!(
        "Upwind: {:.3}, Central: {:.3}, Downwind: {:.3}",
        fields.u.at(4, 5),
        fields.u.at(5, 5),
        fields.u.at(6, 5)
    );

    // Test with positive velocity (flow from left to right)
    let flux_pos = scheme.convective_flux(
        fields.u.at(4, 5), // upwind
        fields.u.at(5, 5), // central
        fields.u.at(6, 5), // downwind
        1.0,               // positive velocity
    );

    // Test with negative velocity (flow from right to left)
    let flux_neg = scheme.convective_flux(
        fields.u.at(4, 5), // upwind
        fields.u.at(5, 5), // central
        fields.u.at(6, 5), // downwind
        -1.0,              // negative velocity
    );

    println!("Superbee limiter results:");
    println!("  Positive velocity flux: {:.3}", flux_pos);
    println!("  Negative velocity flux: {:.3}", flux_neg);

    Ok(())
}

fn demonstrate_muscl2_vanleer(
    fields: &SimulationFields<f64>,
) -> Result<(), Box<dyn std::error::Error>> {
    println!("\n--- MUSCL2 with van Leer Limiter ---");

    let limiter = VanLeer;
    let muscl = MusclScheme::new(limiter, MusclOrder::SecondOrder);
    let scheme = MusclDiscretization::new(muscl);

    // Test with positive velocity
    let flux_pos = scheme.convective_flux(
        fields.u.at(4, 5), // upwind
        fields.u.at(5, 5), // central
        fields.u.at(6, 5), // downwind
        1.0,               // positive velocity
    );

    println!("van Leer limiter results:");
    println!("  Positive velocity flux: {:.3}", flux_pos);
    println!("  (Balanced accuracy/stability compromise)");

    Ok(())
}
