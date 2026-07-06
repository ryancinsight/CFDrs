//! Discretization coefficients for STANDARD algorithm
//!
//! # Invariant (Coefficient Positivity)
//!
//! All neighbour coefficients $a_{nb} \ge 0$ under upwind convection,
//! and $a_P = \sum_{nb} a_{nb} + S_P$ with $S_P \ge 0$. This ensures
//! diagonal dominance and the discrete maximum principle.

use crate::scalar;
use crate::scalar::Cfd2dScalar;
use eunomia::NumericElement;
use serde::{Deserialize, Serialize};

/// Cell coefficients for discretized momentum equation
/// Following Patankar notation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellCoefficients<T: Cfd2dScalar + Copy> {
    /// Central coefficient (diagonal)
    pub ap: T,
    /// East neighbor coefficient
    pub ae: T,
    /// West neighbor coefficient
    pub aw: T,
    /// North neighbor coefficient
    pub an: T,
    /// South neighbor coefficient
    pub as_: T,
    /// Source term
    pub su: T,
    /// Pressure gradient coefficient
    pub d: T,
}

impl<T: Cfd2dScalar + Copy + NumericElement> CellCoefficients<T> {
    /// Create zero coefficients
    #[must_use]
    pub fn zero() -> Self {
        Self {
            ap: scalar::zero::<T>(),
            ae: scalar::zero::<T>(),
            aw: scalar::zero::<T>(),
            an: scalar::zero::<T>(),
            as_: scalar::zero::<T>(),
            su: scalar::zero::<T>(),
            d: scalar::zero::<T>(),
        }
    }

    /// Calculate central coefficient from neighbors
    pub fn calculate_ap(&mut self) {
        self.ap = self.ae + self.aw + self.an + self.as_;
    }

    /// Apply under-relaxation
    pub fn apply_relaxation(&mut self, alpha: T) {
        let one = scalar::one::<T>();
        self.ap /= alpha;
        self.su += ((one - alpha) * self.ap) * self.d;
    }
}
