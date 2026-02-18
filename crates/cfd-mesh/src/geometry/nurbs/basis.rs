//! # B-Spline Basis Function Evaluation
//!
//! Implements the **Cox–de Boor recursion** for computing B-spline basis
//! functions `N_{i,p}(ξ)` and their first derivatives `N'_{i,p}(ξ)`.
//!
//! ## Algorithm — de Boor's Algorithm
//!
//! Rather than computing all `N_{i,p}` recursively (which recomputes shared
//! sub-expressions), this module uses the *triangular table* approach:
//!
//! ```text
//! N[0] = [N_{i-p,0}, N_{i-p+1,0}, …, N_{i,0}]  (only the active span is 1)
//!
//! For d = 1 to p:
//!   For j = 0 to p-d:
//!     left  = ξ − ξ_{i-p+d+j}
//!     right = ξ_{i+1+j} − ξ
//!     saved = N[j] * left / (left + right)
//!     N[j]  = saved + N[j+1] * right / (left + right)
//! ```
//!
//! This yields the `p+1` non-zero basis functions `N_{i-p,p}(ξ)` through
//! `N_{i,p}(ξ)` in O(p²) time.
//!
//! ## Theorem — Partition of Unity
//!
//! For any valid knot span and any `ξ` in `[ξᵢ, ξᵢ₊₁)`:
//!
//! ```text
//! Σⱼ₌₀ᵖ N_{i-p+j, p}(ξ) = 1
//! ```
//!
//! All evaluations in this module maintain this invariant.

use crate::core::scalar::Real;
use super::knot::KnotVector;

/// Evaluate the `p+1` non-zero B-spline basis functions at parameter `t`.
///
/// Returns a fixed-size array — use [`eval_basis_alloc`] if `p` is not
/// known at compile time.
///
/// # Arguments
/// - `span`  — knot span index from [`KnotVector::find_span`].
/// - `t`     — parameter value.
/// - `p`     — degree.
/// - `knots` — the knot vector.
///
/// # Returns
/// `N[0..=p]` where `N[j] = N_{span-p+j, p}(t)`.
pub fn eval_basis(span: usize, t: Real, p: usize, knots: &KnotVector) -> Vec<Real> {
    let mut n = vec![0.0 as Real; p + 1];
    let mut left = vec![0.0 as Real; p + 1];
    let mut right = vec![0.0 as Real; p + 1];

    n[0] = 1.0;
    for j in 1..=p {
        left[j] = t - knots.get(span + 1 - j);
        right[j] = knots.get(span + j) - t;
        let mut saved = 0.0 as Real;
        for r in 0..j {
            let denom = right[r + 1] + left[j - r];
            let temp = if denom.abs() < 1e-15 {
                0.0
            } else {
                n[r] / denom
            };
            n[r] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        n[j] = saved;
    }
    n
}

/// Evaluate all non-zero basis functions **and their first derivatives**.
///
/// Returns `(N, dN)` where:
/// - `N[j] = N_{span-p+j, p}(t)`
/// - `dN[j] = N'_{span-p+j, p}(t)`
///
/// Uses the recurrence:
/// ```text
/// N'_{i,p}(t) = p * (N_{i,p-1}(t)     / (ξᵢ₊ₚ   − ξᵢ)
///                  − N_{i+1,p-1}(t)   / (ξᵢ₊ₚ₊₁ − ξᵢ₊₁))
/// ```
pub fn eval_basis_and_deriv(
    span: usize,
    t: Real,
    p: usize,
    knots: &KnotVector,
) -> (Vec<Real>, Vec<Real>) {
    // ndu[i][j]: N_{span-p+j, i}(t) triangular table
    let mut ndu = vec![vec![0.0 as Real; p + 1]; p + 1];
    let mut left = vec![0.0 as Real; p + 1];
    let mut right = vec![0.0 as Real; p + 1];

    ndu[0][0] = 1.0;
    for j in 1..=p {
        left[j] = t - knots.get(span + 1 - j);
        right[j] = knots.get(span + j) - t;
        let mut saved = 0.0 as Real;
        for r in 0..j {
            ndu[j][r] = right[r + 1] + left[j - r];
            let temp = if ndu[j][r].abs() < 1e-15 {
                0.0
            } else {
                ndu[r][j - 1] / ndu[j][r]
            };
            ndu[r][j] = saved + right[r + 1] * temp;
            saved = left[j - r] * temp;
        }
        ndu[j][j] = saved;
    }

    let n: Vec<Real> = ndu[p].to_vec();

    // Compute first-order derivatives
    let mut dn = vec![0.0 as Real; p + 1];
    let mut a = vec![vec![0.0 as Real; p + 1]; 2];
    a[0][0] = 1.0;
    for j in 0..=p {
        let mut s1 = 0usize;
        let mut s2 = 1usize;
        let r = (span as isize - p as isize + j as isize) as usize;
        for k in 1..=1usize {
            let mut d = 0.0 as Real;
            let rk = r as isize - k as isize;
            let pk = p as isize - k as isize;
            if r >= k {
                a[s2][0] = a[s1][0] / ndu[pk as usize + 1][r - k];
                d = a[s2][0] * ndu[r - k][pk as usize];
            }
            let j1 = if rk >= -1 { 1 } else { (-rk) as usize };
            let j2 = if (r as isize - 1) <= pk {
                k - 1
            } else {
                p - r
            };
            for i in j1..=j2 {
                a[s2][i] = (a[s1][i] - a[s1][i - 1]) / ndu[pk as usize + 1][r + i - k];
                d += a[s2][i] * ndu[r + i - k][pk as usize];
            }
            if r <= p - k {
                a[s2][k] = -a[s1][k - 1] / ndu[pk as usize + 1][r];
                d += a[s2][k] * ndu[r][pk as usize];
            }
            dn[j] = d;
            std::mem::swap(&mut s1, &mut s2);
        }
    }
    // Multiply by p
    let pp = p as Real;
    for d in dn.iter_mut() {
        *d *= pp;
    }

    (n, dn)
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::geometry::nurbs::knot::KnotVector;

    #[test]
    fn partition_of_unity() {
        // Cubic with 5 control points, clamped uniform
        let kv = KnotVector::clamped_uniform(4, 3);
        for &t in &[0.0, 0.1, 0.25, 0.5, 0.75, 0.9, 1.0] {
            let span = kv.find_span(t, 4);
            let n = eval_basis(span, t, 3, &kv);
            let sum: f64 = n.iter().sum();
            assert!(
                (sum - 1.0).abs() < 1e-12,
                "partition of unity violated at t={}: sum={}",
                t,
                sum
            );
        }
    }

    #[test]
    fn non_negativity() {
        let kv = KnotVector::clamped_uniform(5, 3);
        for i in 0..20 {
            let t = i as f64 / 20.0;
            let span = kv.find_span(t, 5);
            let n = eval_basis(span, t, 3, &kv);
            for (j, &v) in n.iter().enumerate() {
                assert!(v >= -1e-14, "basis N[{}] = {} < 0 at t={}", j, v, t);
            }
        }
    }

    #[test]
    fn linear_basis_midpoint() {
        // Linear (p=1), 2 control points: knots = [0,0,1,1]
        let kv = KnotVector::try_new(vec![0.0, 0.0, 1.0, 1.0]).unwrap();
        let span = kv.find_span(0.5, 1);
        let n = eval_basis(span, 0.5, 1, &kv);
        // Both basis functions should be 0.5 at t=0.5
        assert!((n[0] - 0.5).abs() < 1e-14);
        assert!((n[1] - 0.5).abs() < 1e-14);
    }
}
