//! Time-stepping methods for ODE integration in CFD simulations.
//!
//! This module provides various time integration schemes optimized for
//! CFD applications, including explicit and implicit methods, adaptive
//! time stepping, and specialized schemes for stiff equations.

pub mod adaptive;
pub mod exponential;
pub mod imex;
pub mod rk_chebyshev;
pub mod runge_kutta;
pub mod stability;
pub mod traits;

// Re-export main types for convenience
pub use adaptive::AdaptiveTimeStepper;
pub use exponential::{ExponentialConfig, ExponentialRungeKutta4, ExponentialTimeDifferencing};
pub use imex::IMEXTimeStepper;
pub use rk_chebyshev::{RhsFunction, RkcConfig, RungeKuttaChebyshev};
pub use runge_kutta::{LowStorageRK4, RungeKutta3, RungeKutta4};
pub use stability::{
    CFLAnalysis, NumericalScheme, StabilityAnalyzer, StabilityStatus, VonNeumannAnalysis,
};
pub use traits::{TimeStepController, TimeStepper};

// Type aliases for common schemes

/// Fourth-order Runge-Kutta method (RK4)
pub type RK4<T> = RungeKutta4<T>;

/// Third-order Runge-Kutta method (RK3)
pub type RK3<T> = RungeKutta3<T>;
