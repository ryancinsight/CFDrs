//! Core traits for network analysis

use cfd_core::Result;
use nalgebra::RealField;

/// Trait for domain-specific network analyzers
pub trait NetworkAnalyzer<T: RealField + Copy> {
    /// Analysis result type
    type Result;

    /// Perform analysis on the network
    fn analyze(&mut self, network: &crate::network::Network<T>) -> Result<Self::Result>;

    /// Get analyzer name
    fn name(&self) -> &str;
}
