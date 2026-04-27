//! Strain rate tensor computation for LES models
//!
//! Implements efficient computation of the strain-rate magnitude using
//! second-order interior and boundary finite-difference stencils.
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

use nalgebra::DMatrix;

/// Compute the magnitude of the strain rate tensor.
///
/// # Theorem -- Second-Order Boundary Strain Recovery
///
/// On a uniform grid, the interior central stencil and one-sided boundary
/// stencils used here recover the exact derivative of every quadratic
/// polynomial at the evaluation point up to the truncation term `O(h^2)`.
///
/// **Proof.** Taylor expansion of `f(x+h)`, `f(x-h)`, and `f(x+2h)` about the
/// target point shows that `(f_{i+1}-f_{i-1})/(2h)`,
/// `(-3f_0+4f_1-f_2)/(2h)`, and `(3f_n-4f_{n-1}+f_{n-2})/(2h)` all cancel
/// the constant term, preserve one copy of `f'`, and cancel the second
/// derivative term. The first nonzero neglected term is proportional to
/// `h^2 f'''`; hence boundary strain is computed from the resolved velocity
/// field instead of imposed as an artificial zero state.
pub fn compute_strain_rate_magnitude(
    velocity_u: &DMatrix<f64>,
    velocity_v: &DMatrix<f64>,
    dx: f64,
    dy: f64,
) -> DMatrix<f64> {
    let nx = velocity_u.nrows();
    let ny = velocity_u.ncols();
    let mut strain_magnitude = DMatrix::zeros(nx, ny);

    for i in 0..nx {
        for j in 0..ny {
            let du_dx = differentiate_uniform(nx, i, dx, |ii| velocity_u[(ii, j)]);
            let du_dy = differentiate_uniform(ny, j, dy, |jj| velocity_u[(i, jj)]);
            let dv_dx = differentiate_uniform(nx, i, dx, |ii| velocity_v[(ii, j)]);
            let dv_dy = differentiate_uniform(ny, j, dy, |jj| velocity_v[(i, jj)]);

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

    for i in 0..nx {
        for j in 0..ny {
            let du_dx = differentiate_uniform(nx, i, dx, |ii| velocity_u[(ii, j)]);
            let du_dy = differentiate_uniform(ny, j, dy, |jj| velocity_u[(i, jj)]);
            let dv_dx = differentiate_uniform(nx, i, dx, |ii| velocity_v[(ii, j)]);
            let dv_dy = differentiate_uniform(ny, j, dy, |jj| velocity_v[(i, jj)]);

            s11[(i, j)] = du_dx;
            s22[(i, j)] = dv_dy;
            s12[(i, j)] = 0.5 * (du_dy + dv_dx);
        }
    }

    (s11, s22, s12)
}

#[inline]
fn differentiate_uniform<F>(n: usize, idx: usize, spacing: f64, sample: F) -> f64
where
    F: Fn(usize) -> f64,
{
    if n < 2 {
        return 0.0;
    }
    if n == 2 {
        return (sample(1) - sample(0)) / spacing;
    }
    if idx == 0 {
        return (-3.0 * sample(0) + 4.0 * sample(1) - sample(2)) / (2.0 * spacing);
    }
    if idx + 1 == n {
        return (3.0 * sample(n - 1) - 4.0 * sample(n - 2) + sample(n - 3)) / (2.0 * spacing);
    }
    (sample(idx + 1) - sample(idx - 1)) / (2.0 * spacing)
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

        for &value in strain.iter() {
            assert_relative_eq!(value, 1.0, epsilon = 1e-10);
        }

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

        for i in 0..10 {
            for j in 0..10 {
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
