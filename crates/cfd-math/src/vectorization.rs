//! Vectorization utilities for SIMD operations in CFD computations.
//!
//! This module provides vectorized operations that can take advantage of SIMD
//! instructions for improved performance in numerical computations.

use nalgebra::RealField;
use rayon::prelude::*;

/// Vectorized operations for CFD computations
pub struct VectorizedOps;

impl VectorizedOps {
    /// Vectorized element-wise addition with potential SIMD optimization
    pub fn add_vectorized<T: RealField + Copy + Send + Sync>(
        a: &[T], 
        b: &[T], 
        result: &mut [T]
    ) -> Result<(), &'static str> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err("All slices must have the same length");
        }

        // Use parallel iterators for automatic vectorization
        result.par_iter_mut()
            .zip(a.par_iter().zip(b.par_iter()))
            .for_each(|(r, (a_val, b_val))| {
                *r = a_val + b_val;
            });

        Ok(())
    }

    /// Vectorized element-wise multiplication
    pub fn mul_vectorized<T: RealField + Copy + Send + Sync>(
        a: &[T], 
        b: &[T], 
        result: &mut [T]
    ) -> Result<(), &'static str> {
        if a.len() != b.len() || a.len() != result.len() {
            return Err("All slices must have the same length");
        }

        result.par_iter_mut()
            .zip(a.par_iter().zip(b.par_iter()))
            .for_each(|(r, (a_val, b_val))| {
                *r = a_val * b_val;
            });

        Ok(())
    }

    /// Vectorized scalar multiplication with broadcasting
    pub fn scale_vectorized<T: RealField + Copy + Send + Sync>(
        input: &[T], 
        scalar: T, 
        result: &mut [T]
    ) -> Result<(), &'static str> {
        if input.len() != result.len() {
            return Err("Input and result slices must have the same length");
        }

        result.par_iter_mut()
            .zip(input.par_iter())
            .for_each(|(r, val)| {
                *r = scalar * val;
            });

        Ok(())
    }

    /// Broadcasting addition: adds scalar to each element of vector
    pub fn broadcast_add<T: RealField + Copy + Send + Sync>(
        input: &[T],
        scalar: T,
        result: &mut [T]
    ) -> Result<(), &'static str> {
        if input.len() != result.len() {
            return Err("Input and result slices must have the same length");
        }

        // Use iterator chunks for cache-friendly access patterns
        const CHUNK_SIZE: usize = 64; // Optimize for cache line size

        input.par_chunks(CHUNK_SIZE)
            .zip(result.par_chunks_mut(CHUNK_SIZE))
            .for_each(|(input_chunk, result_chunk)| {
                input_chunk.iter()
                    .zip(result_chunk.iter_mut())
                    .for_each(|(input_val, result_val)| {
                        *result_val = input_val + scalar;
                    });
            });

        Ok(())
    }

    /// Broadcasting multiplication: multiplies each element by broadcasted vector
    pub fn broadcast_mul_vector<T: RealField + Copy + Send + Sync>(
        matrix: &[T], 
        matrix_rows: usize,
        matrix_cols: usize,
        vector: &[T],
        result: &mut [T]
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
        result.par_chunks_mut(matrix_cols)
            .zip(matrix.par_chunks(matrix_cols))
            .for_each(|(result_row, matrix_row)| {
                result_row.iter_mut()
                    .zip(matrix_row.iter())
                    .zip(vector.iter())
                    .for_each(|((r, m), v)| {
                        *r = m * v;
                    });
            });

        Ok(())
    }

    /// Vectorized dot product with chunked processing for cache efficiency
    pub fn dot_vectorized<T: RealField + Copy + Send + Sync>(a: &[T], b: &[T]) -> Result<T, &'static str> {
        if a.len() != b.len() {
            return Err("Vectors must have the same length");
        }

        // Use chunked parallel processing for better cache locality
        const CHUNK_SIZE: usize = 1024;
        
        let result = a.par_chunks(CHUNK_SIZE)
            .zip(b.par_chunks(CHUNK_SIZE))
            .map(|(a_chunk, b_chunk)| {
                a_chunk.iter()
                    .zip(b_chunk.iter())
                    .map(|(x, y)| x * y)
                    .fold(T::zero(), |acc, val| acc + val)
            })
            .reduce(|| T::zero(), |acc, val| acc + val);

        Ok(result)
    }

    /// Vectorized L2 norm computation
    pub fn l2_norm_vectorized<T: RealField + Copy + Send + Sync>(input: &[T]) -> T {
        const CHUNK_SIZE: usize = 1024;
        
        let sum_of_squares = input.par_chunks(CHUNK_SIZE)
            .map(|chunk| {
                chunk.iter()
                    .map(|x| x * x)
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
        y: &mut [T]
    ) -> Result<(), &'static str> {
        if row_ptr.len() != y.len() + 1 {
            return Err("Invalid row pointer length");
        }

        // Parallel processing of matrix rows
        y.par_iter_mut()
            .enumerate()
            .for_each(|(i, y_i)| {
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
        result: &mut [T]
    ) -> Result<(), &'static str> {
        if input.len() < 2 || result.len() != input.len() - 1 {
            return Err("Invalid input or result length for finite differences");
        }

        result.par_iter_mut()
            .enumerate()
            .for_each(|(i, diff)| {
                *diff = (input[i + 1] - input[i]) / spacing;
            });

        Ok(())
    }

    /// Vectorized convolution for filtering operations
    pub fn convolution_vectorized<T: RealField + Copy + Send + Sync>(
        signal: &[T],
        kernel: &[T],
        result: &mut [T]
    ) -> Result<(), &'static str> {
        let signal_len = signal.len();
        let kernel_len = kernel.len();
        let result_len = signal_len + kernel_len - 1;

        if result.len() != result_len {
            return Err("Result length must be signal_len + kernel_len - 1");
        }

        // Parallel convolution computation
        result.par_iter_mut()
            .enumerate()
            .for_each(|(n, output)| {
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
    pub fn reduce_vectorized<T, F>(input: &[T], identity: T, op: F) -> T
    where
        T: Clone + Send + Sync,
        F: Fn(T, T) -> T + Sync + Send + Clone,
    {
        const CHUNK_SIZE: usize = 1024;

        input.par_chunks(CHUNK_SIZE)
            .map(|chunk| {
                chunk.iter()
                    .cloned()
                    .fold(identity.clone(), |a, b| op(a, b))
            })
            .reduce(|| identity.clone(), |a, b| op(a, b))
    }

    /// Vectorized prefix sum (scan) operation
    pub fn prefix_sum_vectorized<T: RealField + Copy + Send + Sync>(
        input: &[T],
        result: &mut [T]
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

/// Vectorized operations specifically for CFD stencil computations
pub struct StencilOps;

impl StencilOps {
    /// 5-point stencil for 2D Laplacian (vectorized)
    pub fn laplacian_5point<T: RealField + Copy + Send + Sync>(
        field: &[T],
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        result: &mut [T]
    ) -> Result<(), &'static str> {
        if field.len() != nx * ny || result.len() != nx * ny {
            return Err("Field and result dimensions must match grid size");
        }

        let dx2 = dx * dx;
        let dy2 = dy * dy;

        // Process interior points (sequential for now due to mutable access patterns)
        for j in 1..ny-1 {
            for i in 1..nx-1 {
                let idx = j * nx + i;
                let center = field[idx];
                let left = field[idx - 1];
                let right = field[idx + 1];
                let bottom = field[idx - nx];
                let top = field[idx + nx];

                result[idx] = (left - center * T::from_f64(2.0).unwrap_or_else(|| T::zero()) + right) / dx2
                            + (bottom - center * T::from_f64(2.0).unwrap_or_else(|| T::zero()) + top) / dy2;
            }
        }

        Ok(())
    }

    /// Vectorized gradient computation using central differences
    pub fn gradient_central<T: RealField + Copy + Send + Sync>(
        field: &[T],
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        grad_x: &mut [T],
        grad_y: &mut [T]
    ) -> Result<(), &'static str> {
        if field.len() != nx * ny || grad_x.len() != nx * ny || grad_y.len() != nx * ny {
            return Err("All arrays must match grid size");
        }

        let dx_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dx);
        let dy_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dy);

        // Gradient computation (sequential for now due to mutable access patterns)
        for j in 1..ny-1 {
            for i in 1..nx-1 {
                let idx = j * nx + i;

                // X-gradient
                grad_x[idx] = (field[idx + 1] - field[idx - 1]) * dx_inv;

                // Y-gradient
                grad_y[idx] = (field[idx + nx] - field[idx - nx]) * dy_inv;
            }
        }

        Ok(())
    }

    /// Vectorized divergence computation for 3D vector fields on structured grids
    pub fn divergence_3d<T: RealField + Copy + Send + Sync>(
        u_field: &[T], v_field: &[T], w_field: &[T],
        nx: usize, ny: usize, nz: usize,
        dx: T, dy: T, dz: T,
        result: &mut [T]
    ) -> Result<(), &'static str> {
        if u_field.len() != nx * ny * nz ||
           v_field.len() != nx * ny * nz ||
           w_field.len() != nx * ny * nz ||
           result.len() != nx * ny * nz {
            return Err("All fields must match grid dimensions");
        }

        let dx_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dx);
        let dy_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dy);
        let dz_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dz);

        // Compute divergence: ∇·v = ∂u/∂x + ∂v/∂y + ∂w/∂z
        for k in 1..nz-1 {
            for j in 1..ny-1 {
                for i in 1..nx-1 {
                    let idx = k * nx * ny + j * nx + i;

                    // Central differences
                    let dudx = (u_field[idx + 1] - u_field[idx - 1]) * dx_inv;
                    let dvdy = (v_field[idx + nx] - v_field[idx - nx]) * dy_inv;
                    let dwdz = (w_field[idx + nx * ny] - w_field[idx - nx * ny]) * dz_inv;

                    result[idx] = dudx + dvdy + dwdz;
                }
            }
        }

        Ok(())
    }

    /// Vectorized curl computation for 3D vector fields
    pub fn curl_3d<T: RealField + Copy + Send + Sync>(
        u_field: &[T], v_field: &[T], w_field: &[T],
        nx: usize, ny: usize, nz: usize,
        dx: T, dy: T, dz: T,
        curl_x: &mut [T], curl_y: &mut [T], curl_z: &mut [T]
    ) -> Result<(), &'static str> {
        if u_field.len() != nx * ny * nz ||
           v_field.len() != nx * ny * nz ||
           w_field.len() != nx * ny * nz ||
           curl_x.len() != nx * ny * nz ||
           curl_y.len() != nx * ny * nz ||
           curl_z.len() != nx * ny * nz {
            return Err("All fields must match grid dimensions");
        }

        let dx_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dx);
        let dy_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dy);
        let dz_inv = T::one() / (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * dz);

        // Compute curl: ∇×v = (∂w/∂y - ∂v/∂z, ∂u/∂z - ∂w/∂x, ∂v/∂x - ∂u/∂y)
        for k in 1..nz-1 {
            for j in 1..ny-1 {
                for i in 1..nx-1 {
                    let idx = k * nx * ny + j * nx + i;

                    // Curl x-component: ∂w/∂y - ∂v/∂z
                    let dwdy = (w_field[idx + nx] - w_field[idx - nx]) * dy_inv;
                    let dvdz = (v_field[idx + nx * ny] - v_field[idx - nx * ny]) * dz_inv;
                    curl_x[idx] = dwdy - dvdz;

                    // Curl y-component: ∂u/∂z - ∂w/∂x
                    let dudz = (u_field[idx + nx * ny] - u_field[idx - nx * ny]) * dz_inv;
                    let dwdx = (w_field[idx + 1] - w_field[idx - 1]) * dx_inv;
                    curl_y[idx] = dudz - dwdx;

                    // Curl z-component: ∂v/∂x - ∂u/∂y
                    let dvdx = (v_field[idx + 1] - v_field[idx - 1]) * dx_inv;
                    let dudy = (u_field[idx + nx] - u_field[idx - nx]) * dy_inv;
                    curl_z[idx] = dvdx - dudy;
                }
            }
        }

        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_vectorized_add() {
        let a = vec![1.0, 2.0, 3.0, 4.0];
        let b = vec![5.0, 6.0, 7.0, 8.0];
        let mut result = vec![0.0; 4];

        VectorizedOps::add_vectorized(&a, &b, &mut result).expect("CRITICAL: Add proper error handling");
        assert_eq!(result, vec![6.0, 8.0, 10.0, 12.0]);
    }

    #[test]
    fn test_vectorized_dot() {
        let a = vec![1.0, 2.0, 3.0];
        let b = vec![4.0, 5.0, 6.0];

        let result = VectorizedOps::dot_vectorized(&a, &b).expect("CRITICAL: Add proper error handling");
        assert_eq!(result, 32.0); // 1*4 + 2*5 + 3*6 = 32
    }

    #[test]
    fn test_l2_norm_vectorized() {
        let input = vec![3.0, 4.0];
        let norm = VectorizedOps::l2_norm_vectorized(&input);
        assert!((norm - 5.0_f64).abs() < 1e-10);
    }
}
