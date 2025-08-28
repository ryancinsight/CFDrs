//! Conservation history tracking

use nalgebra::RealField;
use serde::{Deserialize, Serialize};

/// Track conservation properties over time
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ConservationHistory<T: RealField + Copy> {
    /// Time points
    pub times: Vec<T>,
    /// Conservation errors at each time
    pub errors: Vec<T>,
    /// Whether conservation was satisfied at each time
    pub satisfied: Vec<bool>,
}

impl<T: RealField + Copy> ConservationHistory<T> {
    /// Create new conservation history
    pub fn new() -> Self {
        Self {
            times: Vec::new(),
            errors: Vec::new(),
            satisfied: Vec::new(),
        }
    }

    /// Add a conservation check result
    pub fn add_check(&mut self, time: T, error: T, tolerance: T) {
        self.times.push(time);
        self.errors.push(error);
        self.satisfied.push(error <= tolerance);
    }

    /// Get the maximum error in the history
    pub fn max_error(&self) -> Option<T> {
        self.errors
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal))
            .copied()
    }

    /// Get the fraction of time points where conservation was satisfied
    pub fn satisfaction_rate(&self) -> f64 {
        if self.satisfied.is_empty() {
            1.0
        } else {
            let count = self.satisfied.iter().filter(|&&s| s).count();
            count as f64 / self.satisfied.len() as f64
        }
    }
}

impl<T: RealField + Copy> Default for ConservationHistory<T> {
    fn default() -> Self {
        Self::new()
    }
}
