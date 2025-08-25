//! Discretization coefficients for STANDARD algorithm

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Cell coefficients for discretized momentum equation
/// Following Patankar notation
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellCoefficients<T: RealField + Copy> {
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

impl<T: RealField + Copy> CellCoefficients<T> {
    /// Create zero coefficients
    #[must_use]
    pub fn zero() -> Self {
        Self {
            ap: T::zero(),
            ae: T::zero(),
            aw: T::zero(),
            an: T::zero(),
            as_: T::zero(),
            su: T::zero(),
            d: T::zero(),
        }
    }

    /// Calculate central coefficient from neighbors
    pub fn calculate_ap(&mut self) {
        self.ap = self.ae + self.aw + self.an + self.as_;
    }

    /// Apply under-relaxation
    pub fn apply_relaxation(&mut self, alpha: T) {
        let one = T::one();
        self.ap /= alpha;
        self.su += ((one - alpha) * self.ap) * self.d;
    }
}
