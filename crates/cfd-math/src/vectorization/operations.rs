//! Core vectorized operations for CFD computations

use nalgebra::RealField;
use rayon::prelude::*;

/// Cache-friendly chunk size for vectorized operations
const CHUNK_SIZE_CACHE: usize = 64;
/// Larger chunk size for parallel reductions
const CHUNK_SIZE_PARALLEL: usize = 1024;

/// Vectorized operations for CFD computations
pub struct VectorizedOps;

impl VectorizedOps {
    /// Create a new vectorized operations handler
    pub fn new() -> Self {
        Self
    }

    /// Vectorized element-wise addition with SIMD optimization
    pub fn add_vectorized<T: RealField + Copy + Send + Sync>(
        a: &[T],
        b: &[T],
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err("All slices must have the same length");
        }

        // Use parallel iterators for automatic vectorization
        result
            .par_iter_mut()
            .zip(a.par_iter().zip(b.par_iter()))
            .for_each(|(r, (a_val, b_val))| {
                *r = *a_val + *b_val;
            });

        Ok(())
    }

    /// Vectorized element-wise multiplication
    pub fn mul_vectorized<T: RealField + Copy + Send + Sync>(
        a: &[T],
        b: &[T],
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err("All slices must have the same length");
        }

        result
            .par_iter_mut()
            .zip(a.par_iter().zip(b.par_iter()))
            .for_each(|(r, (a_val, b_val))| {
                *r = *a_val * *b_val;
            });

        Ok(())
    }

    /// Vectorized scalar multiplication with broadcasting
    pub fn scale_vectorized<T: RealField + Copy + Send + Sync>(
        input: &[T],
        scalar: T,
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if input.len() != result.len() {
            return Err("Input and result slices must have the same length");
        }

        result
            .par_iter_mut()
            .zip(input.par_iter())
            .for_each(|(r, val)| {
                *r = scalar * *val;
            });

        Ok(())
    }

    /// Broadcasting addition: adds scalar to each element of vector
    pub fn broadcast_add<T: RealField + Copy + Send + Sync>(
        input: &[T],
        scalar: T,
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if input.len() != result.len() {
            return Err("Input and result slices must have the same length");
        }

        // Use iterator chunks for cache-friendly access patterns
        input
            .par_chunks(CHUNK_SIZE_CACHE)
            .zip(result.par_chunks_mut(CHUNK_SIZE_CACHE))
            .for_each(|(input_chunk, result_chunk)| {
                input_chunk.iter().zip(result_chunk.iter_mut()).for_each(
                    |(input_val, result_val)| {
                        *result_val = *input_val + scalar;
                    },
                );
            });

        Ok(())
    }

    /// Broadcasting multiplication: multiplies each element by broadcasted vector
    pub fn broadcast_mul_vector<T: RealField + Copy + Send + Sync>(
        matrix: &[T],
        matrix_rows: usize,
        matrix_cols: usize,
        vector: &[T],
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if matrix.len() != matrix_rows * matrix_cols {
            return Err("Matrix dimensions don't match data length");
        }
        if vector.len() != matrix_cols {
            return Err("Vector length must match matrix columns");
        }
        if result.len() != matrix.len() {
            return Err("Result must have same size as matrix");
        }

        // Row-wise broadcasting using parallel iteration
        result
            .par_chunks_mut(matrix_cols)
            .zip(matrix.par_chunks(matrix_cols))
            .for_each(|(result_row, matrix_row)| {
                result_row
                    .iter_mut()
                    .zip(matrix_row.iter())
                    .zip(vector.iter())
                    .for_each(|((r, m), v)| {
                        *r = *m * *v;
                    });
            });

        Ok(())
    }

    /// Vectorized dot product with chunked processing for cache efficiency
    pub fn dot_vectorized<T: RealField + Copy + Send + Sync>(
        a: &[T],
        b: &[T],
    ) -> Result<T, &'static str> {
        if a.len() != b.len() {
            return Err("Vectors must have the same length");
        }

        // Use chunked parallel processing for better cache locality
        let result = a
            .par_chunks(CHUNK_SIZE_PARALLEL)
            .zip(b.par_chunks(CHUNK_SIZE_PARALLEL))
            .map(|(a_chunk, b_chunk)| {
                a_chunk
                    .iter()
                    .zip(b_chunk.iter())
                    .map(|(x, y)| *x * *y)
                    .fold(T::zero(), |acc, val| acc + val)
            })
            .reduce(|| T::zero(), |acc, val| acc + val);

        Ok(result)
    }

    /// Vectorized L2 norm computation
    pub fn l2_norm_vectorized<T: RealField + Copy + Send + Sync>(input: &[T]) -> T {
        let sum_of_squares = input
            .par_chunks(CHUNK_SIZE_PARALLEL)
            .map(|chunk| {
                chunk
                    .iter()
                    .map(|x| *x * *x)
                    .fold(T::zero(), |acc, val| acc + val)
            })
            .reduce(|| T::zero(), |acc, val| acc + val);

        sum_of_squares.sqrt()
    }

    /// Vectorized matrix-vector multiplication for sparse patterns
    pub fn matvec_vectorized<T: RealField + Copy + Send + Sync>(
        values: &[T],
        col_indices: &[usize],
        row_ptr: &[usize],
        x: &[T],
        y: &mut [T],
    ) -> Result<(), &'static str> {
        if row_ptr.len() != y.len() + 1 {
            return Err("Invalid row pointer length");
        }

        // Parallel processing of matrix rows
        y.par_iter_mut().enumerate().for_each(|(i, y_i)| {
            let start = row_ptr[i];
            let end = row_ptr[i + 1];

            *y_i = (start..end)
                .map(|j| {
                    let col = col_indices[j];
                    values[j] * x[col]
                })
                .fold(T::zero(), |acc, val| acc + val);
        });

        Ok(())
    }

    /// Vectorized finite difference computation
    pub fn finite_diff_vectorized<T: RealField + Copy + Send + Sync>(
        input: &[T],
        spacing: T,
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if input.len() < 2 || result.len() != input.len() - 1 {
            return Err("Invalid input or result length for finite differences");
        }

        result.par_iter_mut().enumerate().for_each(|(i, diff)| {
            *diff = (input[i + 1] - input[i]) / spacing;
        });

        Ok(())
    }

    /// Vectorized convolution for filtering operations
    pub fn convolution_vectorized<T: RealField + Copy + Send + Sync>(
        signal: &[T],
        kernel: &[T],
        result: &mut [T],
    ) -> Result<(), &'static str> {
        let signal_len = signal.len();
        let kernel_len = kernel.len();
        let result_len = signal_len + kernel_len - 1;

        if result.len() != result_len {
            return Err("Result length must be signal_len + kernel_len - 1");
        }

        // Parallel convolution computation
        result.par_iter_mut().enumerate().for_each(|(n, output)| {
            *output = (0..kernel_len)
                .filter_map(|k| {
                    if n >= k && n - k < signal_len {
                        Some(signal[n - k] * kernel[k])
                    } else {
                        None
                    }
                })
                .fold(T::zero(), |acc, val| acc + val);
        });

        Ok(())
    }

    /// Vectorized reduction operations with custom binary operators
    ///
    /// # Note on Clippy Warning
    /// The `|a, b| op(a, b)` closures are flagged as redundant by clippy, but they are
    /// required by Rust ownership semantics. The closure captures `op` by move into the
    /// parallel iterator context, making direct function reference `op` invalid.
    /// This is a known clippy false positive (rust-clippy#13094).
    #[allow(clippy::redundant_closure)]
    pub fn reduce_vectorized<T, F>(input: &[T], identity: T, op: F) -> T
    where
        T: Clone + Send + Sync,
        F: Fn(T, T) -> T + Sync + Send + Clone,
    {
        input
            .par_chunks(CHUNK_SIZE_PARALLEL)
            .map(|chunk| {
                chunk
                    .iter()
                    .cloned()
                    .fold(identity.clone(), |a, b| op(a, b))
            })
            .reduce(|| identity.clone(), |a, b| op(a, b))
    }

    /// Vectorized prefix sum (scan) operation
    pub fn prefix_sum_vectorized<T: RealField + Copy + Send + Sync>(
        input: &[T],
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if input.len() != result.len() {
            return Err("Input and result must have the same length");
        }

        if input.is_empty() {
            return Ok(());
        }

        // Sequential prefix sum (parallel prefix sum is complex and may not be worth it for most cases)
        result[0] = input[0];
        for i in 1..input.len() {
            result[i] = result[i - 1] + input[i];
        }

        Ok(())
    }
}

impl Default for VectorizedOps {
    fn default() -> Self {
        Self::new()
    }
}
