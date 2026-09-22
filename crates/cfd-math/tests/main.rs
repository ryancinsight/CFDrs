//! Hierarchical integration harness for math-kernel contracts.
//!
//! The leaf modules retain their original contract assertions untouched;
//! one Cargo target replaces the previous flat target-per-file topology
//! (7 binaries linking the stack per build). Nextest still isolates per
//! test, and the committed timeout/serial-group profile applies unchanged:
//! test-name filters match the trailing test name with or without the
//! harness module prefix.

#[path = "main/amg_coarsening_tests.rs"]
mod amg_coarsening_tests;
#[path = "main/amg_integration_test.rs"]
mod amg_integration_test;
#[path = "main/core_solver_tests.rs"]
mod core_solver_tests;
#[path = "main/preconditioner_edge_cases.rs"]
mod preconditioner_edge_cases;
#[path = "main/simd_tests.rs"]
mod simd_tests;
#[path = "main/simple_gmres_tests.rs"]
mod simple_gmres_tests;
#[path = "main/window_iterator_test.rs"]
mod window_iterator_test;
