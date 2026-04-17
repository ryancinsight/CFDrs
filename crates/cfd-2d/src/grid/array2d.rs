//! `Array2D<T>` — flat, contiguous 2D array for cache-efficient field storage.
//!
//! # Motivation
//!
//! The SIMPLE solver iterates over 2D fields millions of times per solve.
//! `Vec<Vec<T>>` (jagged arrays) stores each row as a separate heap allocation,
//! causing:
//! - O(rows) allocations instead of 1
//! - Cache misses on row transitions (non-contiguous memory)
//! - Poor hardware prefetcher utilization
//!
//! # Theorem (Spatial Locality)
//!
//! For a row-major contiguous array of size R×C, sequential access to elements
//! `a[i][j], a[i][j+1], ..., a[i][C-1], a[i+1][0], ...` produces a stride-1
//! memory access pattern. Modern CPUs prefetch L1 cache lines of 64 bytes
//! (8 `f64` values), so every 8th element access triggers a cache miss in the
//! best case. For jagged arrays, *every row transition* is a potential cache
//! miss because rows reside at unrelated heap addresses.
//!
//! **Proof**: Let `base` be the array's starting address. Element `(i, j)` is
//! stored at address `base + (i * C + j) * sizeof(T)`. For consecutive elements
//! `(i, j)` and `(i, j+1)`, the address difference is `sizeof(T)` — exactly
//! stride-1. At the row boundary, `(i, C-1)` to `(i+1, 0)` also has stride-1
//! because `(i * C + C - 1) + 1 = (i+1) * C + 0`. Thus the entire array is a
//! single contiguous block with no gaps, maximizing prefetcher effectiveness.
//!
//! **Expected impact**: For an N×M grid, reduces heap allocations from
//! N × fields + fields to just `fields` (one per Array2D). Eliminates all
//! row-transition cache misses.

use std::ops::{Index, IndexMut};

/// A flat, row-major 2D array.
///
/// Elements are stored in a single contiguous `Vec<T>` with row-major layout:
/// element `(i, j)` is at index `i * cols + j`.
#[derive(Debug, Clone)]
pub struct Array2D<T> {
    data: Vec<T>,
    rows: usize,
    cols: usize,
}

impl<T: Clone> Array2D<T> {
    /// Create an `rows × cols` array filled with `val`.
    pub fn new(rows: usize, cols: usize, val: T) -> Self {
        Self {
            data: vec![val; rows * cols],
            rows,
            cols,
        }
    }

    /// Number of rows.
    #[inline]
    pub fn rows(&self) -> usize {
        self.rows
    }

    /// Number of columns.
    #[inline]
    pub fn cols(&self) -> usize {
        self.cols
    }

    /// Fill all elements with `val`.
    #[inline]
    pub fn fill(&mut self, val: T) {
        self.data.fill(val);
    }

    /// Copy data from another `Array2D` of the same dimensions.
    ///
    /// Uses `clone_from_slice` which compiles to `memcpy` for `Copy` types —
    /// zero allocation, O(rows × cols) time.
    ///
    /// # Panics
    /// Panics if dimensions differ.
    pub fn copy_from(&mut self, other: &Self) {
        assert_eq!(self.rows, other.rows, "Array2D::copy_from: row mismatch");
        assert_eq!(self.cols, other.cols, "Array2D::copy_from: col mismatch");
        self.data.clone_from_slice(&other.data);
    }

    /// Get an immutable reference to element `(i, j)`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> &T {
        debug_assert!(i < self.rows && j < self.cols);
        &self.data[i * self.cols + j]
    }

    /// Get a mutable reference to element `(i, j)`.
    #[inline]
    pub fn get_mut(&mut self, i: usize, j: usize) -> &mut T {
        debug_assert!(i < self.rows && j < self.cols);
        &mut self.data[i * self.cols + j]
    }

    /// View the underlying flat slice.
    #[inline]
    pub fn as_slice(&self) -> &[T] {
        &self.data
    }

    /// View the underlying flat mutable slice.
    #[inline]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        &mut self.data
    }

    /// Iterate over all elements (row-major order).
    #[inline]
    pub fn iter(&self) -> impl Iterator<Item = &T> {
        self.data.iter()
    }
}

impl<T> Index<(usize, usize)> for Array2D<T> {
    type Output = T;
    #[inline]
    fn index(&self, (i, j): (usize, usize)) -> &T {
        debug_assert!(i < self.rows && j < self.cols);
        &self.data[i * self.cols + j]
    }
}

impl<T> IndexMut<(usize, usize)> for Array2D<T> {
    #[inline]
    fn index_mut(&mut self, (i, j): (usize, usize)) -> &mut T {
        debug_assert!(i < self.rows && j < self.cols);
        &mut self.data[i * self.cols + j]
    }
}

/// A flat, row-major 2D boolean mask.
///
/// Stored as `Vec<bool>` for direct compatibility with Array2D indexing
/// patterns while keeping the type distinct from numeric arrays.
pub type Mask2D = Array2D<bool>;

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_new_and_dimensions() {
        let a: Array2D<f64> = Array2D::new(3, 4, 0.0);
        assert_eq!(a.rows(), 3);
        assert_eq!(a.cols(), 4);
        assert_eq!(a.as_slice().len(), 12);
    }

    #[test]
    fn test_indexing() {
        let mut a = Array2D::new(3, 4, 0.0_f64);
        a[(1, 2)] = 42.0;
        assert_eq!(a[(1, 2)], 42.0);
        assert_eq!(*a.get(1, 2), 42.0);
    }

    #[test]
    fn test_fill() {
        let mut a = Array2D::new(2, 3, 0.0_f64);
        a[(0, 0)] = 1.0;
        a.fill(5.0);
        for val in a.iter() {
            assert_eq!(*val, 5.0);
        }
    }

    #[test]
    fn test_copy_from() {
        let src = Array2D::new(2, 3, 7.0_f64);
        let mut dst = Array2D::new(2, 3, 0.0);
        dst.copy_from(&src);
        for val in dst.iter() {
            assert_eq!(*val, 7.0);
        }
    }

    #[test]
    #[should_panic(expected = "row mismatch")]
    fn test_copy_from_dimension_mismatch() {
        let src = Array2D::new(2, 3, 0.0_f64);
        let mut dst = Array2D::new(3, 3, 0.0);
        dst.copy_from(&src);
    }

    #[test]
    fn test_contiguous_layout() {
        let mut a = Array2D::new(3, 4, 0.0_f64);
        // Write in row-major order and verify flat layout
        for i in 0..3 {
            for j in 0..4 {
                a[(i, j)] = (i * 4 + j) as f64;
            }
        }
        let expected: Vec<f64> = (0..12).map(|x| x as f64).collect();
        assert_eq!(a.as_slice(), &expected[..]);
    }

    #[test]
    fn test_mask() {
        let mut mask: Mask2D = Array2D::new(5, 5, true);
        mask[(2, 3)] = false;
        assert!(!mask[(2, 3)]);
        assert!(mask[(0, 0)]);
    }
}
