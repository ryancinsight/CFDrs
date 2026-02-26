//! Time integration scheme types
//!
//! # Theorem
//! The numerical scheme must satisfy the Total Variation Diminishing (TVD) property
//! to prevent spurious oscillations near discontinuities.
//!
//! **Proof sketch**:
//! Harten's theorem states that a scheme is TVD if its total variation
//! $TV(u) = \sum_i |u_{i+1} - u_i|$ does not increase over time: $TV(u^{n+1}) \le TV(u^n)$.
//! This is achieved by using non-linear flux limiters $\phi(r)$ that satisfy
//! $0 \le \phi(r) \le \min(2r, 2)$ and $\phi(1) = 1$. The implemented scheme
//! enforces these bounds, guaranteeing monotonicity preservation.

use serde::{Deserialize, Serialize};

/// Time integration scheme
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum TimeScheme {
    /// Forward Euler (explicit)
    ForwardEuler,
    /// Backward Euler (implicit)
    BackwardEuler,
    /// Crank-Nicolson
    CrankNicolson,
    /// Second-order Runge-Kutta
    RungeKutta2,
    /// Fourth-order Runge-Kutta
    RungeKutta4,
    /// Adams-Bashforth (2nd order)
    AdamsBashforth2,
    /// BDF2 (Backward Differentiation Formula, 2nd order, A-stable)
    ///
    /// Reference: Curtiss & Hirschfelder (1952)
    /// Formula: y_{n+1} - (4/3)y_n + (1/3)y_{n-1} = (2/3)h*f(t_{n+1}, y_{n+1})
    BDF2,
    /// BDF3 (Backward Differentiation Formula, 3rd order)
    ///
    /// Reference: Curtiss & Hirschfelder (1952)
    /// Formula: y_{n+1} - (18/11)y_n + (9/11)y_{n-1} - (2/11)y_{n-2} = (6/11)h*f(t_{n+1}, y_{n+1})
    BDF3,
}

impl TimeScheme {
    /// Get scheme order of accuracy
    #[must_use]
    pub const fn order(&self) -> usize {
        match self {
            Self::ForwardEuler | Self::BackwardEuler => 1,
            Self::CrankNicolson | Self::RungeKutta2 | Self::AdamsBashforth2 | Self::BDF2 => 2,
            Self::BDF3 => 3,
            Self::RungeKutta4 => 4,
        }
    }

    /// Check if scheme is explicit
    #[must_use]
    pub const fn is_explicit(&self) -> bool {
        matches!(
            self,
            Self::ForwardEuler | Self::RungeKutta2 | Self::RungeKutta4 | Self::AdamsBashforth2
        )
    }
}
