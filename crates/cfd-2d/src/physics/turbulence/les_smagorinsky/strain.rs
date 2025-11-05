//! Strain rate tensor computation for LES models
//!
//! Implements efficient computation of the strain rate magnitude
//! using central differences and proper boundary conditions.

use nalgebra::DMatrix;

/// Compute the magnitude of the strain rate tensor
///
/// Uses central differences for velocity gradients and applies
/// zero strain rate boundary conditions at domain boundaries.
pub fn compute_strain_rate_magnitude(
    velocity_u: &DMatrix<f64>,
    velocity_v: &DMatrix<f64>,
    dx: f64,
    dy: f64,
) -> DMatrix<f64> {
    let nx = velocity_u.nrows();
    let ny = velocity_u.ncols();
    let mut strain_magnitude = DMatrix::zeros(nx, ny);

    // Interior points - central differences
    for i in 1..nx-1 {
        for j in 1..ny-1 {
            // Velocity gradients (central differences)
            let du_dx = (velocity_u[(i+1, j)] - velocity_u[(i-1, j)]) / (2.0 * dx);
            let du_dy = (velocity_u[(i, j+1)] - velocity_u[(i, j-1)]) / (2.0 * dy);
            let dv_dx = (velocity_v[(i+1, j)] - velocity_v[(i-1, j)]) / (2.0 * dx);
            let dv_dy = (velocity_v[(i, j+1)] - velocity_v[(i, j-1)]) / (2.0 * dy);

            // Strain rate tensor components
            let s11 = du_dx;
            let s22 = dv_dy;
            let s12 = 0.5 * (du_dy + dv_dx);

            // Magnitude of strain rate tensor - optimized computation
            let s11_sq = s11 * s11;
            let s22_sq = s22 * s22;
            let s12_sq = 4.0 * s12 * s12; // 2*s12*s12 * 2

            strain_magnitude[(i, j)] = (2.0 * (s11_sq + s22_sq) + s12_sq).sqrt();
        }
    }

    // Boundary conditions (zero strain at boundaries)
    for i in 0..nx {
        strain_magnitude[(i, 0)] = 0.0;
        strain_magnitude[(i, ny-1)] = 0.0;
    }
    for j in 0..ny {
        strain_magnitude[(0, j)] = 0.0;
        strain_magnitude[(nx-1, j)] = 0.0;
    }

    strain_magnitude
}

/// Compute individual strain rate components
///
/// Returns (s11, s22, s12) components of the strain rate tensor.
/// Useful for more advanced SGS models that need full tensor information.
pub fn compute_strain_rate_components(
    velocity_u: &DMatrix<f64>,
    velocity_v: &DMatrix<f64>,
    dx: f64,
    dy: f64,
) -> (DMatrix<f64>, DMatrix<f64>, DMatrix<f64>) {
    let nx = velocity_u.nrows();
    let ny = velocity_u.ncols();
    let mut s11 = DMatrix::zeros(nx, ny);
    let mut s22 = DMatrix::zeros(nx, ny);
    let mut s12 = DMatrix::zeros(nx, ny);

    // Interior points - central differences
    for i in 1..nx-1 {
        for j in 1..ny-1 {
            // Velocity gradients (central differences)
            let du_dx = (velocity_u[(i+1, j)] - velocity_u[(i-1, j)]) / (2.0 * dx);
            let du_dy = (velocity_u[(i, j+1)] - velocity_u[(i, j-1)]) / (2.0 * dy);
            let dv_dx = (velocity_v[(i+1, j)] - velocity_v[(i-1, j)]) / (2.0 * dx);
            let dv_dy = (velocity_v[(i, j+1)] - velocity_v[(i, j-1)]) / (2.0 * dy);

            s11[(i, j)] = du_dx;
            s22[(i, j)] = dv_dy;
            s12[(i, j)] = 0.5 * (du_dy + dv_dx);
        }
    }

    (s11, s22, s12)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn create_test_velocity_fields(nx: usize, ny: usize) -> (DMatrix<f64>, DMatrix<f64>) {
        let mut velocity_u = DMatrix::zeros(nx, ny);
        let mut velocity_v = DMatrix::zeros(nx, ny);

        // Simple shear flow: u = y, v = 0
        for i in 0..nx {
            for j in 0..ny {
                velocity_u[(i, j)] = j as f64 * 0.1; // Linear shear
                velocity_v[(i, j)] = 0.0;
            }
        }

        (velocity_u, velocity_v)
    }

    #[test]
    fn test_strain_rate_computation() {
        let (velocity_u, velocity_v) = create_test_velocity_fields(10, 10);
        let strain = compute_strain_rate_magnitude(&velocity_u, &velocity_v, 0.1, 0.1);

        // Check that strain rate is computed (non-zero for shear flow)
        assert!(strain.iter().any(|&s| s > 0.0));

        // Check boundary conditions (should be zero)
        assert_eq!(strain[(0, 0)], 0.0);
        assert_eq!(strain[(9, 0)], 0.0);
        assert_eq!(strain[(0, 9)], 0.0);
        assert_eq!(strain[(9, 9)], 0.0);

        // Check dimensions
        assert_eq!(strain.nrows(), 10);
        assert_eq!(strain.ncols(), 10);
    }

    #[test]
    fn test_strain_rate_components() {
        let (velocity_u, velocity_v) = create_test_velocity_fields(10, 10);
        let (s11, s22, s12) = compute_strain_rate_components(&velocity_u, &velocity_v, 0.1, 0.1);

        // For shear flow u = 0.1*j (j is column index, representing y), v = 0:
        // With dy = 0.1, the discrete derivative gives:
        // du/dy = 1.0 (from finite difference), dv/dx = 0, dv/dy = 0
        // s11 = du/dx = 0, s22 = dv/dy = 0, s12 = 0.5*(du/dy + dv/dx) = 0.5

        // Check interior points
        for i in 1..9 {
            for j in 1..9 {
                assert_relative_eq!(s11[(i, j)], 0.0, epsilon = 1e-10);
                assert_relative_eq!(s22[(i, j)], 0.0, epsilon = 1e-10);
                assert_relative_eq!(s12[(i, j)], 0.5, epsilon = 1e-6);
            }
        }

        // Check dimensions
        assert_eq!(s11.nrows(), 10);
        assert_eq!(s11.ncols(), 10);
        assert_eq!(s22.nrows(), 10);
        assert_eq!(s22.ncols(), 10);
        assert_eq!(s12.nrows(), 10);
        assert_eq!(s12.ncols(), 10);
    }

    #[test]
    fn test_zero_velocity_field() {
        let velocity_u = DMatrix::zeros(10, 10);
        let velocity_v = DMatrix::zeros(10, 10);
        let strain = compute_strain_rate_magnitude(&velocity_u, &velocity_v, 0.1, 0.1);

        // Zero velocity field should give zero strain everywhere
        for &s in strain.iter() {
            assert_relative_eq!(s, 0.0, epsilon = 1e-15);
        }
    }
}
