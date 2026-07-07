use std::ops::{Index, IndexMut};

use eunomia::NumericElement;
use leto::Array2;

/// Minimal 2D matrix wrapper used to decouple cfd-validation internals from legacy dense backends.
#[derive(Debug, Clone)]
pub struct DMatrix<T> {
    data: Array2<T>,
}

impl<T> DMatrix<T> {
    /// Create a zero-initialized matrix with `(rows, cols)` shape.
    #[inline]
    pub fn zeros(rows: usize, cols: usize) -> Self
    where
        T: NumericElement,
    {
        Self {
            data: Array2::zeros([rows, cols]),
        }
    }

    #[cfg(test)]
    /// Create a matrix by evaluating `f(row, col)` for each entry.
    #[inline]
    pub(crate) fn from_fn<F>(rows: usize, cols: usize, mut f: F) -> Self
    where
        F: FnMut(usize, usize) -> T,
    {
        Self {
            data: Array2::from_shape_fn([rows, cols], |[i, j]| f(i, j)),
        }
    }

    /// Return `(rows, cols)` for this matrix.
    #[inline]
    pub fn shape(&self) -> (usize, usize) {
        let shape = self.data.shape();
        (shape[0], shape[1])
    }

    /// Return the number of rows.
    #[inline]
    pub fn nrows(&self) -> usize {
        self.data.shape()[0]
    }

    /// Return the number of columns.
    #[inline]
    pub fn ncols(&self) -> usize {
        self.data.shape()[1]
    }

    /// Fill every matrix entry with `value`.
    #[inline]
    pub fn fill(&mut self, value: T)
    where
        T: Copy,
    {
        let (rows, cols) = self.shape();
        for i in 0..rows {
            for j in 0..cols {
                self[(i, j)] = value;
            }
        }
    }

    /// Copy all values from `other` into `self`.
    #[inline]
    pub fn copy_from(&mut self, other: &Self)
    where
        T: Clone,
    {
        self.data.clone_from(&other.data);
    }
}

impl<T> Index<(usize, usize)> for DMatrix<T> {
    type Output = T;

    #[inline]
    fn index(&self, index: (usize, usize)) -> &Self::Output {
        &self.data[[index.0, index.1]]
    }
}

impl<T> IndexMut<(usize, usize)> for DMatrix<T> {
    #[inline]
    fn index_mut(&mut self, index: (usize, usize)) -> &mut Self::Output {
        &mut self.data[[index.0, index.1]]
    }
}
