//! Physics validation tests against literature benchmarks

use cfd_validation::literature::{
    patankar_1980::PatankarLidDrivenCavity,
    LiteratureValidation,
};

#[test]
fn validate_lid_driven_cavity() {
    let test = PatankarLidDrivenCavity::new(100.0, 32);
    let result = test.validate().expect("Validation should succeed");
    
    assert!(result.passed, "Lid-driven cavity validation failed");
    assert!(result.max_error < test.expected_accuracy(), 
            "Error {} exceeds expected accuracy {}", 
            result.max_error, test.expected_accuracy());
}

#[test]
fn validate_rhie_chow_interpolation() {
    // Test that Rhie-Chow prevents checkerboard pressure
    use cfd_2d::pressure_velocity::rhie_chow_complete::RhieChowInterpolation;
    use cfd_2d::grid::StructuredGrid2D;
    use cfd_2d::fields::Field2D;
    use nalgebra::Vector2;
    
    let grid = StructuredGrid2D::<f64>::new(10, 10, 1.0, 1.0);
    let mut interpolator = RhieChowInterpolation::new(&grid);
    
    // Create checkerboard pressure field
    let mut pressure = Field2D::new(10, 10, 0.0);
    for i in 0..10 {
        for j in 0..10 {
            pressure.set(i, j, if (i + j) % 2 == 0 { 1.0 } else { -1.0 });
        }
    }
    
    // Create uniform velocity field
    let velocity = Field2D::new(10, 10, Vector2::new(1.0, 0.0));
    
    // Rhie-Chow should smooth out the checkerboard
    let face_velocity = interpolator.face_velocity_x(&velocity, &pressure, 0.1, 5, 5);
    
    // Face velocity should be close to uniform (1.0)
    assert!((face_velocity - 1.0).abs() < 0.5, 
            "Rhie-Chow failed to prevent checkerboard: face_velocity = {}", face_velocity);
}

#[test]
fn validate_supg_stabilization() {
    use cfd_3d::fem::stabilization::StabilizationParameters;
    use nalgebra::Vector3;
    
    // Test case: High Peclet number flow
    let h = 0.1;  // Element size
    let nu = 0.001;  // Low viscosity (high Pe)
    let velocity = Vector3::new(1.0, 0.0, 0.0);
    
    let params = StabilizationParameters::new(h, nu, velocity, None);
    
    // Check Peclet number
    let pe = params.peclet_number();
    assert!(pe > 10.0, "Peclet number should be high: {}", pe);
    
    // Check that tau is reasonable
    let tau = params.tau_supg();
    assert!(tau > 0.0 && tau < h, "SUPG tau out of range: {}", tau);
    
    // For high Pe, optimal tau should be close to tau_supg
    let optimal = params.optimal_tau();
    assert!((optimal - tau).abs() < 0.1 * tau, 
            "Optimal tau {} differs too much from SUPG tau {}", optimal, tau);
}

#[test]
fn validate_weno_accuracy() {
    use cfd_2d::schemes::{WENO5, SpatialDiscretization, Grid2D as SchemeGrid2D};
    
    // Test WENO5 on smooth function: sin(x)
    let mut grid = SchemeGrid2D::new(100, 1, 0.01, 1.0, 3);
    
    // Initialize with sin(x)
    for i in 0..106 {
        let x = (i as f64 - 3.0) * 0.01;
        grid.data[(i, 3)] = x.sin();
    }
    
    let weno = WENO5::new();
    
    // Compute derivative at center
    let derivative = weno.compute_derivative(&grid, 50, 3);
    
    // Analytical derivative of sin(x) is cos(x)
    let x = 0.47;  // x at i=50
    let exact = x.cos();
    
    // WENO5 should be 5th order accurate for smooth functions
    let error = (derivative - exact).abs();
    assert!(error < 1e-4, "WENO5 error {} too large for smooth function", error);
}

#[test]
fn validate_energy_conservation() {
    use cfd_2d::energy::EnergyEquationSolver;
    use cfd_2d::grid::StructuredGrid2D;
    use cfd_2d::fields::SimulationFields;
    use cfd_core::Fluid;
    
    let grid = StructuredGrid2D::new(10, 10, 0.1, 0.1);
    let fluid = Fluid::water().expect("Should create water");
    let mut fields = SimulationFields::new(10, 10);
    
    // Set initial uniform temperature
    for i in 0..10 {
        for j in 0..10 {
            fields.temperature.set(i, j, 300.0);
        }
    }
    
    let solver = EnergyEquationSolver::new(&grid, &fluid);
    
    // With no heat sources and insulated boundaries, total energy should be conserved
    let initial_energy: f64 = fields.temperature.data.iter().sum();
    
    // Advance one time step with no sources
    solver.solve(&mut fields, 0.01).expect("Energy solve should succeed");
    
    let final_energy: f64 = fields.temperature.data.iter().sum();
    
    // Energy should be conserved (within numerical precision)
    let energy_change = (final_energy - initial_energy).abs() / initial_energy;
    assert!(energy_change < 1e-10, 
            "Energy not conserved: change = {}", energy_change);
}