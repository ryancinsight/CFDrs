//! Time-stepping methods for ODE integration in CFD simulations.
//!
//! This module provides various time integration schemes optimized for
//! CFD applications, including explicit and implicit methods, adaptive
//! time stepping, and specialized schemes for stiff equations.

pub mod runge_kutta;
pub mod adaptive;
pub mod imex;
pub mod traits;

// Re-export main types for convenience
pub use runge_kutta::{RungeKutta4, RungeKutta3, LowStorageRK4};
pub use adaptive::AdaptiveTimeStepper;
pub use imex::IMEXTimeStepper;
pub use traits::{TimeStepper, TimeStepController};

// Type aliases for common schemes
pub type RK4<T> = RungeKutta4<T>;
pub type RK3<T> = RungeKutta3<T>;
