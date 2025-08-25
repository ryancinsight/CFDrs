//! Fluid dynamics domain - Core fluid mechanics concepts and operations.
//!
//! This module encapsulates all fluid dynamics-specific knowledge following DDD principles.
//! Organized into focused submodules adhering to SLAP and SOC.

pub mod fields;
pub mod flow_regimes;
pub mod operations;
pub mod rans;
pub mod turbulence;

// Re-export core types for convenience
pub use fields::{FlowField, PressureField, ScalarField, VelocityField};
pub use flow_regimes::{FlowClassifier, FlowRegime};
pub use operations::FlowOperations;
pub use rans::{KEpsilonConstants, KEpsilonModel, RANSModel};
pub use turbulence::{MixingLengthModel, SmagorinskyModel, TurbulenceModel};
