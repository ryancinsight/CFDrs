//! Hierarchical integration harness for 3D solver contracts.
//!
//! The leaf modules retain their original contract assertions untouched;
//! one Cargo target replaces the previous flat target-per-file topology
//! (18 binaries linking the stack per build). Nextest still isolates per
//! test, and the committed timeout/serial-group profile applies unchanged:
//! test-name filters match the trailing test name with or without the
//! harness module prefix.

#[path = "main/blueprint_integration.rs"]
mod blueprint_integration;
#[path = "main/cascade_solver_test.rs"]
mod cascade_solver_test;
#[path = "main/cavitation_solver_validation.rs"]
mod cavitation_solver_validation;
#[path = "main/domain_solver_validation.rs"]
mod domain_solver_validation;
#[path = "main/fem_boundary_conditions.rs"]
mod fem_boundary_conditions;
#[path = "main/fem_tests.rs"]
mod fem_tests;
#[path = "main/fourier_validation.rs"]
mod fourier_validation;
#[path = "main/ibm_tests.rs"]
mod ibm_tests;
#[path = "main/level_set_tests.rs"]
mod level_set_tests;
#[path = "main/poiseuille_test.rs"]
mod poiseuille_test;
#[path = "main/poisson_validation.rs"]
mod poisson_validation;
#[path = "main/projection_solver_validation.rs"]
mod projection_solver_validation;
#[path = "main/property_tests.rs"]
mod property_tests;
#[path = "main/robustness_tests.rs"]
mod robustness_tests;
#[path = "main/smagorinsky_test.rs"]
mod smagorinsky_test;
#[path = "main/turbulence_math_tests.rs"]
mod turbulence_math_tests;
#[path = "main/venturi_blood_test.rs"]
mod venturi_blood_test;
#[path = "main/vof_tests.rs"]
mod vof_tests;
