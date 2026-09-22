//! Hierarchical integration harness for optimization-pipeline contracts.
//!
//! The leaf modules retain their original contract assertions untouched;
//! one Cargo target replaces the previous flat target-per-file topology
//! (4 binaries linking the stack per build). Nextest still isolates per
//! test, and the committed timeout/serial-group profile applies unchanged:
//! test-name filters match the trailing test name with or without the
//! harness module prefix.

#[path = "main/milestone12_ga.rs"]
mod milestone12_ga;
#[path = "main/milestone12_option1.rs"]
mod milestone12_option1;
#[path = "main/milestone12_option2.rs"]
mod milestone12_option2;
#[path = "main/schematic_export_annotations.rs"]
mod schematic_export_annotations;
