//! Givens rotations for GMRES
//!
//! Efficient incremental solution of least-squares problem via Givens rotations

use cfd_core::error::Result;
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
        let mut sum = g[i];
        for j in (i + 1)..k {
            sum -= h[(i, j)] * y[j];
        }
        y[i] = sum / h[(i, i)];
    }
    Ok(y)
}

/// Apply Givens rotation to eliminate H[i+1, i]
///
/// Transforms the Hessenberg matrix H to upper triangular form using
/// Givens rotations, enabling efficient solution of the least-squares problem.
///
/// # Arguments
///
/// * `h` - Upper Hessenberg matrix H
/// * `cs` - Cosine values of previous Givens rotations
/// * `sn` - Sine values of previous Givens rotations
/// * `g` - Right-hand side vector for least-squares problem
/// * `i` - Current column index
///
/// # Algorithm
///
/// 1. Apply previous rotations to column i of H
/// 2. Compute new rotation (cs[i], sn[i]) to zero H[i+1, i]
/// 3. Apply new rotation to H and g
///
/// Updates H to upper triangular form and propagates rotation to g
pub fn apply_givens_rotation<T: RealField + Copy>(
    h: &mut DMatrix<T>,
    cs: &mut [T],
    sn: &mut [T],
    g: &mut DVector<T>,
    i: usize,
) {
    // Apply previous Givens rotations to column i
    for k in 0..i {
        let temp = cs[k] * h[(k, i)] + sn[k] * h[(k + 1, i)];
        h[(k + 1, i)] = -sn[k] * h[(k, i)] + cs[k] * h[(k + 1, i)];
        h[(k, i)] = temp;
    }

    // Compute new Givens rotation to eliminate H[i+1, i]
    let h_ii = h[(i, i)];
    let h_ip1_i = h[(i + 1, i)];
    let hypotenuse = (h_ii * h_ii + h_ip1_i * h_ip1_i).sqrt();

    if hypotenuse > T::zero() {
        cs[i] = h_ii / hypotenuse;
        sn[i] = h_ip1_i / hypotenuse;
    } else {
        // Degenerate case: both entries are zero
        cs[i] = T::one();
        sn[i] = T::zero();
    }

    // Apply rotation to H
    h[(i, i)] = cs[i] * h_ii + sn[i] * h_ip1_i;
    h[(i + 1, i)] = T::zero();

    // Apply rotation to RHS vector
    let temp = cs[i] * g[i] + sn[i] * g[i + 1];
    g[i + 1] = -sn[i] * g[i] + cs[i] * g[i + 1];
    g[i] = temp;
}

/// Back-substitution to solve upper triangular system H*y = g
///
/// Solves R*y = g where R = H[0:k, 0:k] is upper triangular after Givens rotations
///
/// # Arguments
///
/// * `h` - Upper triangular Hessenberg matrix (after Givens rotations)
/// * `g` - Right-hand side vector
/// * `y` - Solution vector (output)
/// * `k` - Dimension of the system
pub fn back_substitution<T: RealField + Copy>(
    h: &DMatrix<T>,
    g: &DVector<T>,
    y: &mut DVector<T>,
    k: usize,
) {
    // Solve R*y = g where R = H[0:k, 0:k] is upper triangular
    for i in (0..k).rev() {
        let mut sum = g[i];
        for j in (i + 1)..k {
            sum -= h[(i, j)] * y[j];
        }
        y[i] = sum / h[(i, i)];
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_back_substitution() {
        // Test solving upper triangular system
        // [ 2  1 ] [ y0 ]   [ 5 ]
        // [ 0  3 ] [ y1 ] = [ 6 ]
        // Solution: y1 = 2, y0 = 1.5

        let mut h = DMatrix::zeros(2, 2);
        h[(0, 0)] = 2.0;
        h[(0, 1)] = 1.0;
        h[(1, 1)] = 3.0;

        let g = DVector::from_vec(vec![5.0, 6.0]);
        let mut y = DVector::zeros(2);

        back_substitution(&h, &g, &mut y, 2);

        assert_relative_eq!(y[1], 2.0, epsilon = 1e-10);
        assert_relative_eq!(y[0], 1.5, epsilon = 1e-10);
    }

    #[test]
    fn test_givens_rotation() {
        // Test Givens rotation on simple 2x2 Hessenberg
        let mut h = DMatrix::zeros(3, 2);
        h[(0, 0)] = 3.0;
        h[(1, 0)] = 4.0;

        let mut cs = vec![0.0; 2];
        let mut sn = vec![0.0; 2];
        let mut g = DVector::from_vec(vec![5.0, 0.0, 0.0]);

        apply_givens_rotation(&mut h, &mut cs, &mut sn, &mut g, 0);

        // After rotation, H[1,0] should be zero
        assert_relative_eq!(h[(1, 0)], 0.0, epsilon = 1e-10);

        // H[0,0] should be hypotenuse = 5
        assert_relative_eq!(h[(0, 0)], 5.0, epsilon = 1e-10);

        // Check rotation coefficients
        assert_relative_eq!(cs[0], 0.6, epsilon = 1e-10); // 3/5
        assert_relative_eq!(sn[0], 0.8, epsilon = 1e-10); // 4/5
    }
}
