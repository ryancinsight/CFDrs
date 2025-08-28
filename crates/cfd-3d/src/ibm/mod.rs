//! Immersed Boundary Method (IBM) for complex geometries in 3D
//!
//! The IBM allows simulation of flow around complex objects without
//! body-fitted meshes by using forcing terms in the momentum equations.

pub mod config;
pub mod forcing;
pub mod interpolation;
pub mod lagrangian;
pub mod solver;

pub use config::IbmConfig;
pub use forcing::{DirectForcing, FeedbackForcing, ForcingMethod};
pub use interpolation::{DeltaFunction, InterpolationKernel};
pub use lagrangian::LagrangianPoint;
pub use solver::IbmSolver;
