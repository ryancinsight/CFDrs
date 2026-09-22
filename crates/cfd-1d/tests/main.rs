//! Hierarchical integration harness for 1D network contracts.
//!
//! The leaf modules retain their original contract assertions untouched;
//! one Cargo target replaces the previous flat target-per-file topology
//! (25 binaries linking the stack per build). Nextest still isolates per
//! test, and the committed timeout/serial-group profile applies unchanged:
//! test-name filters match the trailing test name with or without the
//! harness module prefix.

#[path = "main/adversarial_tests.rs"]
mod adversarial_tests;
#[path = "main/analysis_tests.rs"]
mod analysis_tests;
#[path = "main/blueprint_metadata_physics.rs"]
mod blueprint_metadata_physics;
#[path = "main/blueprint_solve_trace.rs"]
mod blueprint_solve_trace;
#[path = "main/blueprint_validation.rs"]
mod blueprint_validation;
#[path = "main/branch_reverse_flow_orientation.rs"]
mod branch_reverse_flow_orientation;
#[path = "main/cell_separation_tests.rs"]
mod cell_separation_tests;
#[path = "main/cell_separation_validation.rs"]
mod cell_separation_validation;
#[path = "main/channel_solver_tests.rs"]
mod channel_solver_tests;
#[path = "main/component_validation.rs"]
mod component_validation;
#[path = "main/components_tests.rs"]
mod components_tests;
#[path = "main/integration_schematics.rs"]
mod integration_schematics;
#[path = "main/manufactured_network.rs"]
mod manufactured_network;
#[path = "main/membrane_component_parity.rs"]
mod membrane_component_parity;
#[path = "main/millifluidics_tests.rs"]
mod millifluidics_tests;
#[path = "main/network_analysis_validation.rs"]
mod network_analysis_validation;
#[path = "main/nonlinear_linearization_and_dirichlet.rs"]
mod nonlinear_linearization_and_dirichlet;
#[path = "main/primary_solve_reliability.rs"]
mod primary_solve_reliability;
#[path = "main/property_tests.rs"]
mod property_tests;
#[path = "main/resistance_model_validation.rs"]
mod resistance_model_validation;
#[path = "main/solver_core_tests.rs"]
mod solver_core_tests;
#[path = "main/spd_and_invariants.rs"]
mod spd_and_invariants;
#[path = "main/transient_composition_parity.rs"]
mod transient_composition_parity;
#[path = "main/transient_droplet_parity.rs"]
mod transient_droplet_parity;
#[path = "main/transient_literature_validation.rs"]
mod transient_literature_validation;
