//! Broadcasting operations for efficient array manipulation
//!
//! Provides zero-copy broadcasting operations similar to `NumPy` broadcasting,
//! enabling efficient element-wise operations on arrays of different shapes.

use nalgebra::RealField;
use std::borrow::Cow;

/// Broadcasting trait for array operations
pub trait Broadcast<T: RealField + Copy> {
    /// Broadcast a scalar to match array shape
    fn broadcast_scalar(&self, scalar: T) -> Cow<'_, [T]>;

    /// Broadcast operation between arrays
    fn broadcast_with<F>(&self, other: &[T], op: F) -> Vec<T>
    where
        F: Fn(T, T) -> T;
}

/// Broadcasting operations for slices
impl<T: RealField + Copy> Broadcast<T> for [T] {
    fn broadcast_scalar(&self, scalar: T) -> Cow<'_, [T]> {
        // If all elements are the same as scalar, return borrowed
        if self.iter().all(|&x| x == scalar) {
            Cow::Borrowed(self)
        } else {
            // Otherwise create new array
            Cow::Owned(vec![scalar; self.len()])
        }
    }

    fn broadcast_with<F>(&self, other: &[T], op: F) -> Vec<T>
    where
        F: Fn(T, T) -> T,
    {
        match (self.len(), other.len()) {
            (n, m) if n == m => {
                // Same size: element-wise operation
                self.iter()
                    .zip(other.iter())
                    .map(|(&a, &b)| op(a, b))
                    .collect()
            }
            (_n, 1) => {
                // Broadcast single element
                let b = other[0];
                self.iter().map(|&a| op(a, b)).collect()
            }
            (1, _m) => {
                // Broadcast self
                let a = self[0];
                other.iter().map(|&b| op(a, b)).collect()
            }
            _ => {
                // Incompatible shapes - use cycling iterator
                self.iter()
                    .cycle()
                    .zip(other.iter().cycle())
                    .take(self.len().max(other.len()))
                    .map(|(&a, &b)| op(a, b))
                    .collect()
            }
        }
    }
}

/// Broadcasting view for efficient operations
pub struct BroadcastView<'a, T: RealField + Copy> {
    data: Cow<'a, [T]>,
    shape: Vec<usize>,
    strides: Vec<usize>,
}

impl<'a, T: RealField + Copy> BroadcastView<'a, T> {
    /// Create a new broadcast view
    pub fn new(data: &'a [T], shape: Vec<usize>) -> Self {
        let strides = Self::compute_strides(&shape);
        Self {
            data: Cow::Borrowed(data),
            shape,
            strides,
        }
    }

    /// Create from owned data
    #[must_use]
    pub fn from_owned(data: Vec<T>, shape: Vec<usize>) -> Self {
        let strides = Self::compute_strides(&shape);
        Self {
            data: Cow::Owned(data),
            shape,
            strides,
        }
    }

    /// Compute strides for given shape
    fn compute_strides(shape: &[usize]) -> Vec<usize> {
        let mut strides = vec![1; shape.len()];
        for i in (0..shape.len() - 1).rev() {
            strides[i] = strides[i + 1] * shape[i + 1];
        }
        strides
    }

    /// Get element at multi-dimensional index
    #[must_use]
    pub fn get(&self, indices: &[usize]) -> Option<T> {
        if indices.len() != self.shape.len() {
            return None;
        }

        let flat_idx = indices
            .iter()
            .zip(&self.strides)
            .map(|(i, s)| i * s)
            .sum::<usize>();

        self.data.get(flat_idx).copied()
    }

    /// Apply operation with broadcasting
    pub fn apply<F>(&self, other: &BroadcastView<'a, T>, op: F) -> Vec<T>
    where
        F: Fn(T, T) -> T,
    {
        // Determine output shape
        let out_shape: Vec<usize> = self
            .shape
            .iter()
            .zip(&other.shape)
            .map(|(&a, &b)| a.max(b))
            .collect();

        let total_size: usize = out_shape.iter().product();
        let mut result = Vec::with_capacity(total_size);

        // Iterate over output indices
        for i in 0..total_size {
            let mut indices = vec![0; out_shape.len()];
            let mut temp = i;
            for (j, &dim) in out_shape.iter().enumerate().rev() {
                indices[j] = temp % dim;
                temp /= dim;
            }

            // Map indices for broadcasting
            let a_indices: Vec<usize> = indices
                .iter()
                .zip(&self.shape)
                .map(|(&i, &s)| if s == 1 { 0 } else { i })
                .collect();

            let b_indices: Vec<usize> = indices
                .iter()
                .zip(&other.shape)
                .map(|(&i, &s)| if s == 1 { 0 } else { i })
                .collect();

            let a = self.get(&a_indices).unwrap_or(T::zero());
            let b = other.get(&b_indices).unwrap_or(T::zero());
            result.push(op(a, b));
        }

        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_broadcast_scalar() {
        let arr = [1.0, 2.0, 3.0];
        let broadcast = arr.broadcast_scalar(5.0);
        
        // Scalar broadcast should always return an owned vector
        assert!(matches!(broadcast, Cow::Owned(_)), 
                "Expected owned vector for scalar broadcast, got borrowed");
        
        if let Cow::Owned(v) = broadcast {
            assert_eq!(v, vec![5.0, 5.0, 5.0]);
        }
    }

    #[test]
    fn test_broadcast_with() {
        let a = [1.0, 2.0, 3.0];
        let b = [2.0];
        let result = a.broadcast_with(&b, |x, y| x + y);
        assert_eq!(result, vec![3.0, 4.0, 5.0]);
    }

    #[test]
    fn test_broadcast_view() {
        let data = vec![1.0, 2.0];
        let view = BroadcastView::new(&data, vec![2, 1]);
        assert_eq!(view.get(&[0, 0]), Some(1.0));
        assert_eq!(view.get(&[1, 0]), Some(2.0));
    }
}
