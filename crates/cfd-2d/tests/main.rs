//! Hierarchical integration harness for 2D solver contracts.
//!
//! The leaf modules retain their original contract assertions untouched;
//! one Cargo target replaces the previous flat target-per-file topology
//! (14 binaries linking the stack per build). Nextest still isolates per
//! test, and the committed timeout/serial-group profile applies unchanged:
//! test-name filters match the trailing test name with or without the
//! harness module prefix.

#[path = "main/blueprint_validation.rs"]
mod blueprint_validation;
#[path = "main/cif_junction_flow_split.rs"]
mod cif_junction_flow_split;
#[path = "main/energy_equation_validation.rs"]
mod energy_equation_validation;
#[path = "main/extended_boundary_conditions_validation.rs"]
mod extended_boundary_conditions_validation;
#[path = "main/ghia_cavity_simplec_validation.rs"]
mod ghia_cavity_simplec_validation;
#[path = "main/momentum_solver_simple_tests.rs"]
mod momentum_solver_simple_tests;
#[path = "main/muscl_validation.rs"]
mod muscl_validation;
#[path = "main/n_furcation_flow_split.rs"]
mod n_furcation_flow_split;
#[path = "main/poisson_fdm_validation.rs"]
mod poisson_fdm_validation;
#[path = "main/primitive_selective_split_tree_cross_fidelity.rs"]
mod primitive_selective_split_tree_cross_fidelity;
#[path = "main/reynolds_stress_comprehensive_tests.rs"]
mod reynolds_stress_comprehensive_tests;
#[path = "main/simplec_pimple_validation.rs"]
mod simplec_pimple_validation;
#[path = "main/spalart_allmaras_comprehensive_tests.rs"]
mod spalart_allmaras_comprehensive_tests;
#[path = "main/tvd_scheme_validation.rs"]
mod tvd_scheme_validation;
