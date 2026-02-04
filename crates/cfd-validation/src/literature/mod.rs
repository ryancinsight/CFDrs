//! Literature-based validation tests for CFD algorithms
//!
//! This module provides validation against published benchmark results
//! to ensure physics accuracy and numerical correctness.

pub mod blood_flow_1d;
pub mod chapman_enskog;
pub mod patankar_1980;

use cfd_core::error::Result;
use nalgebra::RealField;

/// Literature validation test trait
pub trait LiteratureValidation<T: RealField + Copy> {
    /// Run validation test
    fn validate(&self) -> Result<ValidationReport<T>>;

    /// Get reference citation
    fn citation(&self) -> &str;

    /// Get expected accuracy
    fn expected_accuracy(&self) -> T;
}

/// Validation report
#[derive(Debug, Clone)]
pub struct ValidationReport<T: RealField + Copy> {
    /// Test name
    pub test_name: String,
    /// Citation
    pub citation: String,
    /// Maximum error
    pub max_error: T,
    /// Average error
    pub avg_error: T,
    /// Pass/fail status
    pub passed: bool,
    /// Detailed results
    pub details: String,
}
