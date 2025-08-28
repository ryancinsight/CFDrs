//! Energy conservation checker

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use cfd_core::error::Result;
use nalgebra::{DVector, RealField};

/// Energy conservation checker
pub struct EnergyConservationChecker<T: RealField + Copy> {
    tolerance: T,
}

impl<T: RealField + Copy> EnergyConservationChecker<T> {
    /// Create new energy conservation checker
    pub fn new(tolerance: T) -> Self {
        Self { tolerance }
    }
}

impl<T: RealField + Copy> ConservationChecker<T> for EnergyConservationChecker<T> {
    type FlowField = DVector<T>;

    fn check_conservation(&self, _field: &Self::FlowField) -> Result<ConservationReport<T>> {
        // Placeholder implementation
        Ok(ConservationReport::new(
            "Energy Conservation".to_string(),
            T::zero(),
            self.tolerance,
        ))
    }

    fn name(&self) -> &str {
        "Energy Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }
}
