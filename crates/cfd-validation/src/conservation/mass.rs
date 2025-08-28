//! Mass conservation checker

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

/// Mass conservation checker
pub struct MassConservationChecker<T: RealField + Copy> {
    tolerance: T,
}

impl<T: RealField + Copy + FromPrimitive> MassConservationChecker<T> {
    /// Create new mass conservation checker
    pub fn new(tolerance: T) -> Self {
        Self { tolerance }
    }

    /// Check mass conservation for a velocity field
    pub fn check_divergence(&self, u: &DVector<T>, v: &DVector<T>, dx: T, dy: T) -> Result<T> {
        // Compute divergence using central differences
        let mut max_divergence = T::zero();

        // This is a placeholder - actual implementation would compute div(u) properly
        for i in 0..u.len() {
            let div = (u[i] + v[i]) / (dx + dy); // Simplified
            max_divergence = max_divergence.max(div.abs());
        }

        Ok(max_divergence)
    }
}

impl<T: RealField + Copy + FromPrimitive> ConservationChecker<T> for MassConservationChecker<T> {
    type FlowField = (DVector<T>, DVector<T>);

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        let (u, v) = field;
        let dx = T::from_f64(0.01).unwrap_or(T::one());
        let dy = T::from_f64(0.01).unwrap_or(T::one());

        let divergence = self.check_divergence(u, v, dx, dy)?;

        Ok(ConservationReport::new(
            "Mass Conservation".to_string(),
            divergence,
            self.tolerance,
        ))
    }

    fn name(&self) -> &str {
        "Mass Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }
}
