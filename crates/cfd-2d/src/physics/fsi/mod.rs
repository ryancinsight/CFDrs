//! Fluid-structure interaction: the fluid half of the coupling.
//!
//! Atlas ADR 0059 splits Phase 0 coupling so CFDrs produces the interface
//! traction from its own flow state and a structural solver elsewhere consumes
//! it. Neither depends on the other; a coupling driver joins them.

mod traction;

pub use traction::{interface_traction, FaceTraction, InterfaceFace};
