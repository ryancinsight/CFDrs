//! Conservation report structures

use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Report on conservation properties
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationReport<T: RealField + Copy> {
    /// Name of the conservation check
    pub check_name: String,
    /// Whether conservation is satisfied within tolerance
    pub is_conserved: bool,
    /// Conservation error magnitude
    pub error: T,
    /// Tolerance used for the check
    pub tolerance: T,
    /// Additional details about the conservation check
    pub details: HashMap<String, T>,
}

impl<T: RealField + Copy> ConservationReport<T> {
    /// Create a new conservation report
    pub fn new(check_name: String, error: T, tolerance: T) -> Self {
        Self {
            check_name,
            is_conserved: error <= tolerance,
            error,
            tolerance,
            details: HashMap::new(),
        }
    }

    /// Add a detail to the report
    pub fn add_detail(&mut self, key: String, value: T) {
        self.details.insert(key, value);
    }
}
