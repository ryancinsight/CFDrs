//! Fluid dynamics domain - Core fluid mechanics concepts and operations.
//!
//! This module encapsulates all fluid dynamics-specific knowledge following DDD principles.
//! Organized into focused submodules adhering to SLAP and SOC.

pub mod fields;
pub mod turbulence;
pub mod rans;
pub mod flow_regimes;
pub mod operations;

// Re-export core types for convenience
pub use fields::{FlowField, VelocityField, PressureField, ScalarField};
pub use turbulence::{TurbulenceModel, SmagorinskyModel, MixingLengthModel};
pub use rans::{RANSModel, KEpsilonModel, KEpsilonConstants};
pub use flow_regimes::{FlowRegime, FlowClassifier};
pub use operations::FlowOperations;