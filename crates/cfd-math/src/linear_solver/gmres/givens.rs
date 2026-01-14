//! Givens rotations for GMRES
//!
//! Efficient incremental solution of least-squares problem via Givens rotations

use cfd_core::error::{ConvergenceErrorKind, Error, Result};
use nalgebra::{DMatrix, DVector, RealField};

/// Apply previous Givens rotations to the new column of H
pub fn apply_previous_rotations<T: RealField + Copy>(
    h: &mut DMatrix<T>,
    c: &DVector<T>,
    s: &DVector<T>,
    k: usize,
) {
    for i in 0..k {
        let temp = c[i] * h[(i, k)] + s[i] * h[(i + 1, k)];
        h[(i + 1, k)] = -s[i] * h[(i, k)] + c[i] * h[(i + 1, k)];
        h[(i, k)] = temp;
    }
}

/// Compute new Givens rotation coefficients (c, s) to eliminate h[k+1, k]
pub fn compute_rotation<T: RealField + Copy>(h_kk: T, h_kp1_k: T) -> (T, T) {
    if h_kp1_k == T::zero() {
        (T::one(), T::zero())
    } else if h_kk.abs() > h_kp1_k.abs() {
        let t = h_kp1_k / h_kk;
        let cs = T::one() / (T::one() + t * t).sqrt();
        let sn = cs * t;
        (cs, sn)
    } else {
        let t = h_kk / h_kp1_k;
        let sn = T::one() / (T::one() + t * t).sqrt();
        let cs = sn * t;
        (cs, sn)
    }
}

/// Apply new Givens rotation to H and g
pub fn apply_new_rotation<T: RealField + Copy>(
    h: &mut DMatrix<T>,
    g: &mut DVector<T>,
    cs: T,
    sn: T,
    k: usize,
) {
    // Apply to Hessenberg matrix
    h[(k, k)] = cs * h[(k, k)] + sn * h[(k + 1, k)];
    h[(k + 1, k)] = T::zero();

    // Apply to right-hand side vector g
    let temp = cs * g[k] + sn * g[k + 1];
    g[k + 1] = -sn * g[k] + cs * g[k + 1];
    g[k] = temp;
}

/// Solve upper triangular system H[0..k, 0..k] * y = g[0..k]
pub fn solve_upper_triangular<T: RealField + Copy>(
    h: &DMatrix<T>,
    g: &DVector<T>,
    k: usize,
) -> Result<DVector<T>> {
    let mut y = DVector::zeros(k);
    for i in (0..k).rev() {
        let diag = h[(i, i)];
        if diag.abs() <= T::default_epsilon() {
            return Err(Error::Convergence(ConvergenceErrorKind::Breakdown));
        }
        let mut sum = g[i];
        for j in (i + 1)..k {
            sum -= h[(i, j)] * y[j];
        }
        y[i] = sum / diag;
    }
    Ok(y)
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_solve_upper_triangular() {
        // Test solving upper triangular system
        // [ 2  1 ] [ y0 ]   [ 5 ]
        // [ 0  3 ] [ y1 ] = [ 6 ]
        // Solution: y1 = 2, y0 = 1.5

        let mut h = DMatrix::zeros(2, 2);
        h[(0, 0)] = 2.0;
        h[(0, 1)] = 1.0;
        h[(1, 1)] = 3.0;

        let g = DVector::from_vec(vec![5.0, 6.0]);

        let y = solve_upper_triangular(&h, &g, 2).unwrap();
        assert_relative_eq!(y[0], 1.5);
        assert_relative_eq!(y[1], 2.0);
    }

    #[test]
    fn test_givens_rotation() {
        // Test Givens rotation on simple 2x2 Hessenberg
        let mut h = DMatrix::zeros(3, 2);
        h[(0, 0)] = 3.0;
        h[(1, 0)] = 4.0;

        let mut g = DVector::from_vec(vec![5.0, 0.0, 0.0]);

        // Compute rotation to eliminate H[1,0]
        let (cs, sn) = compute_rotation(h[(0, 0)], h[(1, 0)]);

        // Apply rotation to H and g
        apply_new_rotation(&mut h, &mut g, cs, sn, 0);

        // After rotation, H[1,0] should be zero
        assert_relative_eq!(h[(1, 0)], 0.0, epsilon = 1e-10);

        // H[0,0] should be hypotenuse = 5
        assert_relative_eq!(h[(0, 0)], 5.0, epsilon = 1e-10);

        // Check rotation coefficients
        assert_relative_eq!(cs, 0.6, epsilon = 1e-10); // 3/5
        assert_relative_eq!(sn, 0.8, epsilon = 1e-10); // 4/5
    }
}
