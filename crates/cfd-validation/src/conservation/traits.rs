//! Core traits for conservation checking

use super::report::ConservationReport;
use cfd_core::error::Result;
use nalgebra::RealField;

/// Trait for conservation checking
pub trait ConservationChecker<T: RealField + Copy> {
    /// Type representing the flow field
    type FlowField;

    /// Check conservation for the given flow field
    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>>;

    /// Get the name of this conservation check
    fn name(&self) -> &str;

    /// Get the tolerance for this conservation check
    fn tolerance(&self) -> T;
}
