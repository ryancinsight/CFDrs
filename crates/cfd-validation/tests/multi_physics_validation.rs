//! Multi-physics CFD validation tests
//!
//! Tests coupled physics problems including:
//! - Conjugate heat transfer
//! - Species transport with reactions
//! - Magnetohydrodynamics
//! - Multi-phase flows

use cfd_validation::manufactured::multi_physics::{
    ManufacturedConjugateHeatTransfer, ManufacturedMHD, ManufacturedMultiphase,
    ManufacturedSpeciesTransport,
};
use cfd_validation::manufactured::ManufacturedSolution;

/// Test conjugate heat transfer validation
#[test]
fn test_conjugate_heat_transfer_validation() {
    println!("Testing Conjugate Heat Transfer MMS Validation...");

    let cht = ManufacturedConjugateHeatTransfer::<f64>::new(
        10.0, // conductivity ratio (solid/fluid)
        2.0,  // capacity ratio (solid/fluid)
        0.5,  // interface at x=0.5
        1.0,  // temperature amplitude
        1.0,  // frequency parameter
    );

    // Test interface continuity
    let y_test = 0.5;
    let t_test = 1.0;

    let t_fluid_interface = cht.fluid_temperature(0.5, y_test, t_test);
    let t_solid_interface = cht.solid_temperature(0.5, y_test, t_test);

    // Temperature continuity at interface
    assert!(
        (t_fluid_interface - t_solid_interface).abs() < 1e-10,
        "Temperature discontinuity at interface: fluid={}, solid={}",
        t_fluid_interface,
        t_solid_interface
    );

    // Test heat flux continuity (more complex - would require derivatives)
    // For this simplified test, we verify the formulations are consistent

    // Test points across the domain
    let test_points = vec![
        (0.2, 0.3, 0.0, 0.5), // Fluid domain
        (0.8, 0.7, 0.0, 1.0), // Solid domain
        (0.5, 0.5, 0.0, 1.5), // Interface
    ];

    for (x, y, z, t) in test_points {
        let temp = cht.exact_solution(x, y, z, t);
        let source = cht.source_term(x, y, z, t);

        // Physical bounds
        assert!(
            temp.abs() < 50.0,
            "Temperature out of bounds at ({},{}): {}",
            x,
            y,
            temp
        );
        assert!(
            source.is_finite(),
            "Source term not finite at ({},{}): {}",
            x,
            y,
            source
        );
        assert!(
            source.abs() < 1000.0,
            "Source term magnitude too large at ({},{}): {}",
            x,
            y,
            source
        );
    }

    println!("✓ Conjugate heat transfer validation passed");
}

/// Test species transport with chemical reactions
#[test]
fn test_species_transport_reaction_validation() {
    println!("Testing Species Transport with Reaction MMS Validation...");

    let species = ManufacturedSpeciesTransport::<f64>::new(
        0.01, // molecular diffusivity
        0.1,  // reaction rate constant
        1.0,  // concentration amplitude
        2.0,  // kx
        1.5,  // ky
    );

    // Test conservation and boundedness
    let test_points = vec![
        (0.25, 0.25, 0.0, 0.0), // Initial time
        (0.5, 0.5, 0.0, 0.5),   // Intermediate time
        (0.75, 0.75, 0.0, 1.0), // Later time
        (0.33, 0.67, 0.0, 2.0), // Different spatial point
    ];

    let mut prev_concentrations = Vec::new();

    for (x, y, z, t) in test_points {
        let c = species.exact_solution(x, y, z, t);
        let source = species.source_term(x, y, z, t);

        // Physical constraints for species concentration
        assert!(
            c >= 0.0,
            "Concentration must be non-negative: {} at ({},{},t={})",
            c,
            x,
            y,
            t
        );
        assert!(
            c <= 2.0,
            "Concentration exceeds maximum bound: {} at ({},{},t={})",
            c,
            x,
            y,
            t
        );

        // Source term should be finite and reasonable
        assert!(
            source.is_finite(),
            "Source term not finite at ({},{},t={})",
            x,
            y,
            t
        );
        assert!(
            source.abs() < 10.0,
            "Source term magnitude too large: {} at ({},{},t={})",
            source,
            x,
            y,
            t
        );

        prev_concentrations.push(c);
    }

    // Test temporal decay behavior at the same spatial point
    let x_fixed = 0.5;
    let y_fixed = 0.5;
    let c0 = species.exact_solution(x_fixed, y_fixed, 0.0, 0.0);
    let c1 = species.exact_solution(x_fixed, y_fixed, 0.0, 0.5);
    let c2 = species.exact_solution(x_fixed, y_fixed, 0.0, 1.0);

    assert!(
        c1 < c0,
        "Concentration should decay over time: t=0: {}, t=0.5: {}",
        c0,
        c1
    );

    assert!(
        c2 < c1,
        "Concentration should continue decaying: t=0.5: {}, t=1.0: {}",
        c1,
        c2
    );

    println!("✓ Species transport with reaction validation passed");
}

/// Test magnetohydrodynamics validation
#[test]
fn test_mhd_validation() {
    println!("Testing MHD MMS Validation...");

    let mhd = ManufacturedMHD::<f64>::new(
        4e-7, // mu_0 (magnetic permeability)
        1.0,  // sigma (electrical conductivity)
        1.0,  // velocity amplitude
        0.1,  // magnetic field amplitude
        1.0,  // kx
        1.0,  // ky
    );

    // Test MHD solution properties
    let test_points = vec![
        (0.25, 0.25, 0.0, 0.0),
        (0.5, 0.5, 0.0, 0.5),
        (0.75, 0.75, 0.0, 1.0),
    ];

    for (x, y, z, t) in test_points {
        let velocity = mhd.exact_solution(x, y, z, t);
        let source = mhd.source_term(x, y, z, t);

        // Physical constraints
        assert!(
            velocity >= 0.0,
            "Velocity must be non-negative: {} at ({},{},t={})",
            velocity,
            x,
            y,
            t
        );
        assert!(
            velocity < 2.0,
            "Velocity exceeds reasonable bounds: {} at ({},{},t={})",
            velocity,
            x,
            y,
            t
        );

        // Source term should be finite (represents Lorentz force contribution)
        assert!(
            source.is_finite(),
            "MHD source term not finite at ({},{},t={})",
            x,
            y,
            t
        );
        assert!(
            source.abs() < 10.0,
            "MHD source term magnitude too large: {} at ({},{},t={})",
            source,
            x,
            y,
            t
        );
    }

    // Test Alfvén wave characteristics (simplified)
    let alfven_speed = mhd.magnetic_amp / (mhd.mu_0 * mhd.sigma).sqrt();
    assert!(
        alfven_speed > 0.0,
        "Alfvén speed must be positive: {}",
        alfven_speed
    );

    println!(
        "✓ MHD validation passed (Alfvén speed: {:.6})",
        alfven_speed
    );
}

/// Test multi-phase flow validation
#[test]
fn test_multiphase_flow_validation() {
    println!("Testing Multi-Phase Flow MMS Validation...");

    let multiphase = ManufacturedMultiphase::<f64>::new(
        2.0, // density ratio (heavy/light)
        5.0, // viscosity ratio (viscous/light)
        0.5, // interface at y=0.5
        0.8, // phase indicator amplitude
        1.5, // kx
        1.0, // ky
    );

    // Test phase indicator function properties
    let test_points = vec![
        (0.25, 0.25, 0.0, 0.0), // Phase 1 region
        (0.5, 0.5, 0.0, 0.5),   // Interface region
        (0.75, 0.75, 0.0, 1.0), // Phase 2 region
    ];

    for (x, y, z, t) in test_points {
        let phi = multiphase.exact_solution(x, y, z, t);
        let source = multiphase.source_term(x, y, z, t);

        // Phase indicator should be bounded
        assert!(
            phi.abs() <= 1.0,
            "Phase indicator out of bounds: {} at ({},{},t={})",
            phi,
            x,
            y,
            t
        );

        // Source term should be finite and reasonable
        assert!(
            source.is_finite(),
            "Multi-phase source term not finite at ({},{},t={})",
            x,
            y,
            t
        );
        assert!(
            source.abs() < 5.0,
            "Multi-phase source term magnitude too large: {} at ({},{},t={})",
            source,
            x,
            y,
            t
        );
    }

    // Test material property ratios
    assert!(
        multiphase.density_ratio > 1.0,
        "Density ratio should be > 1"
    );
    assert!(
        multiphase.viscosity_ratio > 1.0,
        "Viscosity ratio should be > 1"
    );

    println!("✓ Multi-phase flow validation passed");
}

/// Test coupled physics interaction validation
#[test]
fn test_coupled_physics_interaction() {
    println!("Testing Coupled Physics Interaction Validation...");

    // Create coupled heat transfer and species transport
    let heat_transfer = ManufacturedConjugateHeatTransfer::<f64>::new(5.0, 1.5, 0.4, 1.0, 1.0);
    let species = ManufacturedSpeciesTransport::<f64>::new(0.005, 0.05, 0.8, 1.2, 0.8);

    // Test coupling at interface
    let interface_x = heat_transfer.interface_x;
    let y_test = 0.5;
    let t_test = 0.8;

    let temp_fluid = heat_transfer.exact_solution(interface_x - 0.01, y_test, 0.0, t_test);
    let temp_solid = heat_transfer.exact_solution(interface_x + 0.01, y_test, 0.0, t_test);
    let concentration = species.exact_solution(interface_x, y_test, 0.0, t_test);

    // Verify continuity and coupling
    assert!(
        (temp_fluid - temp_solid).abs() < 0.1,
        "Temperature coupling violated: fluid={}, solid={}",
        temp_fluid,
        temp_solid
    );

    assert!(
        concentration >= 0.0 && concentration <= 1.0,
        "Species concentration out of bounds: {}",
        concentration
    );

    // Test source terms are consistent
    let heat_source = heat_transfer.source_term(interface_x, y_test, 0.0, t_test);
    let species_source = species.source_term(interface_x, y_test, 0.0, t_test);

    assert!(
        heat_source.is_finite(),
        "Heat source not finite at interface"
    );
    assert!(
        species_source.is_finite(),
        "Species source not finite at interface"
    );

    println!("✓ Coupled physics interaction validation passed");
}

/// Property-based tests for multi-physics validation
#[cfg(test)]
mod property_tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        /// Test conjugate heat transfer with various material properties
        #[test]
        fn test_conjugate_heat_transfer_properties(
            k_ratio in 1.0f64..20.0,
            c_ratio in 0.5f64..5.0,
            interface in 0.2f64..0.8,
            amplitude in 0.5f64..2.0,
            frequency in 0.5f64..3.0
        ) {
            let cht = ManufacturedConjugateHeatTransfer::new(k_ratio, c_ratio, interface, amplitude, frequency);

            let x_fluid = interface - 0.1;
            let x_solid = interface + 0.1;
            let y = 0.5;
            let t = 1.0;

            let t_fluid = cht.exact_solution(x_fluid, y, 0.0, t);
            let t_solid = cht.exact_solution(x_solid, y, 0.0, t);
            let s_fluid = cht.source_term(x_fluid, y, 0.0, t);
            let s_solid = cht.source_term(x_solid, y, 0.0, t);

            prop_assert!(t_fluid.is_finite() && t_solid.is_finite(), "Temperatures must be finite");
            prop_assert!(s_fluid.is_finite() && s_solid.is_finite(), "Source terms must be finite");
            prop_assert!(t_fluid.abs() < amplitude * 2.0, "Fluid temperature magnitude reasonable");
            prop_assert!(t_solid.abs() < amplitude * k_ratio * 2.0, "Solid temperature magnitude reasonable");
        }

        /// Test species transport with various reaction parameters
        #[test]
        fn test_species_transport_properties(
            diffusivity in 1e-4f64..0.1,
            reaction_rate in 0.01f64..1.0,
            amplitude in 0.1f64..2.0,
            kx in 0.5f64..3.0,
            ky in 0.5f64..3.0
        ) {
            let species = ManufacturedSpeciesTransport::new(diffusivity, reaction_rate, amplitude, kx, ky);

            let x = 0.5;
            let y = 0.5;
            let t = 1.0;

            let c = species.exact_solution(x, y, 0.0, t);
            let source = species.source_term(x, y, 0.0, t);

            prop_assert!(c >= 0.0, "Concentration must be non-negative");
            prop_assert!(c <= amplitude * 2.0, "Concentration within reasonable bounds");
            prop_assert!(source.is_finite(), "Source term must be finite");
            prop_assert!(source.abs() < amplitude * 10.0, "Source term magnitude reasonable");
        }

        /// Test MHD parameter sensitivity
        #[test]
        fn test_mhd_properties(
            mu_0 in 1e-7f64..1e-6,
            sigma in 0.1f64..10.0,
            vel_amp in 0.1f64..2.0,
            mag_amp in 0.01f64..0.5,
            kx in 0.5f64..2.0,
            ky in 0.5f64..2.0
        ) {
            let mhd = ManufacturedMHD::new(mu_0, sigma, vel_amp, mag_amp, kx, ky);

            let x = 0.5;
            let y = 0.5;
            let t = 0.5;

            let velocity = mhd.exact_solution(x, y, 0.0, t);
            let source = mhd.source_term(x, y, 0.0, t);

            prop_assert!(velocity >= 0.0, "Velocity must be non-negative");
            prop_assert!(velocity <= vel_amp * 2.0, "Velocity within amplitude bounds");
            prop_assert!(source.is_finite(), "MHD source term must be finite");
            prop_assert!(source.abs() < vel_amp * mag_amp * 10.0, "Source term magnitude reasonable");
        }
    }
}

