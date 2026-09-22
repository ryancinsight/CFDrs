//! Hierarchical integration harness for core substrate contracts.
//!
//! The leaf modules retain their original contract assertions untouched;
//! one Cargo target replaces the previous flat target-per-file topology
//! (2 binaries linking the stack per build). Nextest still isolates per
//! test, and the committed timeout/serial-group profile applies unchanged:
//! test-name filters match the trailing test name with or without the
//! harness module prefix.

#[path = "main/backend_validation.rs"]
mod backend_validation;
#[path = "main/gpu_integration_test.rs"]
mod gpu_integration_test;
