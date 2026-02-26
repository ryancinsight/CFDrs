//! Time-stepping methods for ODE integration in CFD simulations.
//!
//! This module provides various time integration schemes optimized for
//! CFD applications, including explicit and implicit methods, adaptive
//! time stepping, and specialized schemes for stiff equations.
//!
//! ## Theorem — Dahlquist Test Equation & A-Stability (Dahlquist 1963)
//!
//! **Theorem**: A time integration method is *A-stable* if and only if, for the
//! Dahlquist test equation dy/dt = λy with Re(λ) < 0, the method gives |y^{n+1}| ≤ |y^n|
//! for all Δt > 0. Among all explicit Runge-Kutta methods of order p:
//!
//! - p = 1 (Forward Euler): A-stability boundary |1 + λΔt| ≤ 1
//! - p = 4 (RK4): stability region extends to Re(λΔt) ≈ -2.79
//!
//! Implicit methods (IMEX, exponential, BDF) can be unconditionally A-stable.
//!
//! ## Theorem — RK4 Stability Region
//!
//! **Theorem**: The classical 4th-order Runge-Kutta stability function is:
//!
//! ```text
//! R(z) = 1 + z + z²/2 + z³/6 + z⁴/24   where z = λΔt
//! ```
//!
//! The stability region {z ∈ ℂ : |R(z)| ≤ 1} on the real axis extends from
//! z ∈ [-2.785, 0], which defines the maximum allowable Δt for a given eigenvalue λ.
//!
//! ## Theorem — CFL Stability Condition
//!
//! **Theorem (Courant-Friedrichs-Lewy)**: For hyperbolic PDEs discretized with
//! explicit time integration, a necessary condition for stability is:
//!
//! ```text
//! CFL = |u| · Δt / Δx ≤ C_max
//! ```
//!
//! where u is the wave speed, Δx the grid spacing, and C_max depends on the
//! spatial scheme: C_max = 1 for upwind, C_max = 1/√d in d dimensions for
//! the D2Q9 LBM scheme, C_max ≈ 2.78/|ρ(A)| for an explicit RK method.
//!
//! This is enforced as a hard invariant: the `AdaptiveTimeStepper` will
//! never advance with a time step violating the configured CFL limit.
//!
//! ## Theorem — Global Error Accumulation
//!
//! **Theorem**: For a one-step p-th order method applied over fixed time T = N·Δt:
//!
//! ```text
//! ‖y(T) - y_N‖ ≤ C · T · Δt^p · max_{t∈[0,T]} ‖y^{(p+1)}(t)‖ / (p+1)!
//! ```
//!
//! The factor T (not Δt^p alone) shows that long integrations require either
//! small Δt or higher-order methods to maintain accuracy.
//!
//! ## Theorem — Runge-Kutta-Chebyshev (RKC) Extended Stability
//!
//! **Theorem (Verwer et al. 2004)**: The s-stage RKC method achieves stability along
//! the real axis up to:
//!
//! ```text
//! |λΔt| ≤ 0.653 · s²
//! ```
//!
//! making it superior for parabolic (diffusion-dominated) problems where
//! standard RK4 requires Δt ~ Δx²/2ν.
//!
//! ## Invariants
//! - Δt > 0 strictly; zero or negative step size is rejected at construction.
//! - CFL number is recorded for each step and must satisfy CFL ≤ CFL_max.
//! - Adaptive steppers never increase Δt by more than a safety factor ≤ 5.

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
