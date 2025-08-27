//! Vectorized stencil operations for CFD computations

use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Constants for stencil operations
const STENCIL_CENTER_COEFFICIENT: f64 = 2.0;
const GRADIENT_DIVISOR: f64 = 2.0;

/// Vectorized operations specifically for CFD stencil computations
pub struct StencilOps;

impl StencilOps {
    /// 5-point stencil for 2D Laplacian (vectorized)
    pub fn laplacian_5point<T: RealField + Copy + Send + Sync + FromPrimitive>(
        field: &[T],
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if field.len() != nx * ny || result.len() != nx * ny {
            return Err("Field and result dimensions must match grid size");
        }

        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let center_coeff = T::from_f64(STENCIL_CENTER_COEFFICIENT).unwrap_or_else(|| T::zero());

        // Process interior points (sequential for now due to mutable access patterns)
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                let center = field[idx];
                let left = field[idx - 1];
                let right = field[idx + 1];
                let bottom = field[idx - nx];
                let top = field[idx + nx];

                result[idx] = (left - center * center_coeff + right) / dx2
                    + (bottom - center * center_coeff + top) / dy2;
            }
        }

        Ok(())
    }

    /// Vectorized gradient computation using central differences
    pub fn gradient_central<T: RealField + Copy + Send + Sync + FromPrimitive>(
        field: &[T],
        nx: usize,
        ny: usize,
        dx: T,
        dy: T,
        grad_x: &mut [T],
        grad_y: &mut [T],
    ) -> Result<(), &'static str> {
        if field.len() != nx * ny || grad_x.len() != nx * ny || grad_y.len() != nx * ny {
            return Err("All arrays must match grid size");
        }

        let gradient_div = T::from_f64(GRADIENT_DIVISOR).unwrap_or_else(|| T::zero());
        let dx_inv = T::one() / (gradient_div * dx);
        let dy_inv = T::one() / (gradient_div * dy);

        // Gradient computation (sequential for now due to mutable access patterns)
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
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
    pub fn divergence_3d<T: RealField + Copy + Send + Sync + FromPrimitive>(
        u_field: &[T],
        v_field: &[T],
        w_field: &[T],
        nx: usize,
        ny: usize,
        nz: usize,
        dx: T,
        dy: T,
        dz: T,
        result: &mut [T],
    ) -> Result<(), &'static str> {
        if u_field.len() != nx * ny * nz
            || v_field.len() != nx * ny * nz
            || w_field.len() != nx * ny * nz
            || result.len() != nx * ny * nz
        {
            return Err("All fields must match grid dimensions");
        }

        let gradient_div = T::from_f64(GRADIENT_DIVISOR).unwrap_or_else(|| T::zero());
        let dx_inv = T::one() / (gradient_div * dx);
        let dy_inv = T::one() / (gradient_div * dy);
        let dz_inv = T::one() / (gradient_div * dz);

        // Compute divergence: ∇·v = ∂u/∂x + ∂v/∂y + ∂w/∂z
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
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
    pub fn curl_3d<T: RealField + Copy + Send + Sync + FromPrimitive>(
        u_field: &[T],
        v_field: &[T],
        w_field: &[T],
        nx: usize,
        ny: usize,
        nz: usize,
        dx: T,
        dy: T,
        dz: T,
        curl_x: &mut [T],
        curl_y: &mut [T],
        curl_z: &mut [T],
    ) -> Result<(), &'static str> {
        if u_field.len() != nx * ny * nz
            || v_field.len() != nx * ny * nz
            || w_field.len() != nx * ny * nz
            || curl_x.len() != nx * ny * nz
            || curl_y.len() != nx * ny * nz
            || curl_z.len() != nx * ny * nz
        {
            return Err("All fields must match grid dimensions");
        }

        let gradient_div = T::from_f64(GRADIENT_DIVISOR).unwrap_or_else(|| T::zero());
        let dx_inv = T::one() / (gradient_div * dx);
        let dy_inv = T::one() / (gradient_div * dy);
        let dz_inv = T::one() / (gradient_div * dz);

        // Compute curl: ∇×v = (∂w/∂y - ∂v/∂z, ∂u/∂z - ∂w/∂x, ∂v/∂x - ∂u/∂y)
        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
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
    fn test_laplacian_5point() {
        // Create a field with known Laplacian
        let nx = 5;
        let ny = 5;
        let dx = 1.0;
        let dy = 1.0;

        let mut field = vec![0.0f64; nx * ny];
        let mut result = vec![0.0f64; nx * ny];

        // Set center point
        field[12] = 1.0; // Center of 5x5 grid

        StencilOps::laplacian_5point(&field, nx, ny, dx, dy, &mut result)
            .expect("Laplacian computation failed");

        // Check that Laplacian is computed correctly
        assert!(result[12].abs() > 0.0);
    }

    #[test]
    fn test_gradient_central() {
        let nx = 5;
        let ny = 5;
        let dx = 1.0;
        let dy = 1.0;

        // Create a field with linear gradient
        let mut field = vec![0.0f64; nx * ny];
        for j in 0..ny {
            for i in 0..nx {
                field[j * nx + i] = i as f64 + j as f64;
            }
        }

        let mut grad_x = vec![0.0f64; nx * ny];
        let mut grad_y = vec![0.0f64; nx * ny];

        StencilOps::gradient_central(&field, nx, ny, dx, dy, &mut grad_x, &mut grad_y)
            .expect("Gradient computation failed");

        // Check interior points have correct gradient
        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                assert!((grad_x[idx] - 1.0).abs() < 1e-10);
                assert!((grad_y[idx] - 1.0).abs() < 1e-10);
            }
        }
    }
}
