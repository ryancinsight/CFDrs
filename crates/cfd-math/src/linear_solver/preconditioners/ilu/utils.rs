//! Utility functions for ILU factorization

use cfd_core::error::{Error, Result};

/// Find index of diagonal element in CSR row
///
/// # Arguments
///
/// * `offsets` - CSR row offsets array
/// * `indices` - CSR column indices array
/// * `row` - Row index to search
///
/// # Returns
///
/// Index into values array where diagonal element is stored
pub fn find_diagonal_index(offsets: &[usize], indices: &[usize], row: usize) -> Result<usize> {
    let row_start = offsets[row];
    let row_end = offsets[row + 1];

    for idx in row_start..row_end {
        if indices[idx] == row {
            return Ok(idx);
        }
    }

    Err(Error::InvalidInput(
        "Diagonal element not found in matrix row".to_string(),
    ))
}
