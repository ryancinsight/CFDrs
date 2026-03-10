//! Streaming operators for LBM.
//!
//! Propagates distribution functions along lattice links using a pull scheme
//! on the flat, contiguous `Vec<T>` representation.
//!
//! # Theorem — Streaming Mass Conservation
//!
//! **Statement**: The pull-scheme streaming step preserves total lattice mass exactly:
//! $\sum_{i,j,q} f_q^{t+1}(i,j) = \sum_{i,j,q} f_q^t(i,j)$.
//!
//! **Proof**:
//!
//! 1. The pull scheme sets $f_q^{t+1}(i,j) = f_q^t(i - e_{q,x},\; j - e_{q,y})$, with
//!    periodic wrap. This is a *permutation* of the source values.
//! 2. A permutation preserves the multiset of values and hence their sum.
//! 3. Therefore, $\sum_{i,j,q} f_q^{t+1} = \sum_{i,j,q} f_q^t$. □
//!
//! # Theorem — No-Slip Bounce-Back Conservation
//!
//! **Statement**: Full bounce-back at a wall node $\mathbf{x}_w$ satisfies
//! $u(\mathbf{x}_w) = 0$ to second-order accuracy (half-way bounce-back).
//!
//! **Proof sketch**:
//! Bounce-back reverses $f_{\bar{q}}(\mathbf{x}_w) \leftarrow f_q(\mathbf{x}_w)$.
//! The macroscopic velocity $u = \sum_q e_q f_q / \rho$ involves equal and opposite
//! $e_q$ terms for each $(q, \bar{q})$ pair, whose contributions cancel exactly,
//! giving $u = 0$. □

use crate::solvers::lbm::lattice::D2Q9;
use nalgebra::RealField;

/// Streaming operator for D2Q9 LBM on a flat contiguous buffer.
///
/// # Layout
///
/// The distribution array `f` is stored as a flat `Vec<T>` with layout:
/// ```text
/// f[j * nx * 9 + i * 9 + q]
/// ```
/// This achieves stride-1 access over the 9 velocity directions at a node,
/// the most frequently accessed dimension in collision and macroscopic updates.
pub struct StreamingOperator;

/// Compute flat index into the distribution buffer.
///
/// # Arguments
/// * `j` — row index (0 ≤ j < ny)
/// * `i` — column index (0 ≤ i < nx)
/// * `q` — velocity direction (0 ≤ q < 9)
/// * `nx` — number of columns
///
/// # Returns
/// `j * nx * 9 + i * 9 + q`
#[inline(always)]
pub fn f_idx(j: usize, i: usize, q: usize, nx: usize) -> usize {
    j * nx * 9 + i * 9 + q
}

impl StreamingOperator {
    /// Pull-scheme streaming step.
    ///
    /// For each node (i, j) and direction q, pulls from the upstream node:
    /// ```text
    /// f_dst[j, i, q] = f_src[(j - e_y) mod ny, (i - e_x) mod nx, q]
    /// ```
    ///
    /// This is a pure permutation (Theorem above), so mass is exactly conserved.
    pub fn stream<T: RealField + Copy>(f_src: &[T], f_dst: &mut [T], nx: usize, ny: usize) {
        for j in 0..ny {
            for i in 0..nx {
                for q in 0..9 {
                    let (ex, ey) = D2Q9::VELOCITIES[q];
                    // Periodic upstream index
                    let src_i = ((i as i32 - ex).rem_euclid(nx as i32)) as usize;
                    let src_j = ((j as i32 - ey).rem_euclid(ny as i32)) as usize;
                    f_dst[f_idx(j, i, q, nx)] = f_src[f_idx(src_j, src_i, q, nx)];
                }
            }
        }
    }

    /// Streaming with interior-only update (boundary nodes excluded via mask).
    ///
    /// The `boundary_mask` is a flat `bool` slice with layout `mask[j*nx + i]`.
    pub fn stream_with_boundaries<T: RealField + Copy>(
        f_src: &[T],
        f_dst: &mut [T],
        boundary_mask: &[bool],
        nx: usize,
        ny: usize,
    ) {
        for j in 0..ny {
            for i in 0..nx {
                if boundary_mask[j * nx + i] {
                    continue; // boundary nodes handled separately
                }
                for q in 0..9 {
                    let (ex, ey) = D2Q9::VELOCITIES[q];
                    let src_i = i as i32 - ex;
                    let src_j = j as i32 - ey;
                    if src_i >= 0 && src_i < nx as i32 && src_j >= 0 && src_j < ny as i32 {
                        let si = src_i as usize;
                        let sj = src_j as usize;
                        if !boundary_mask[sj * nx + si] {
                            f_dst[f_idx(j, i, q, nx)] = f_src[f_idx(sj, si, q, nx)];
                        }
                    }
                }
            }
        }
    }

    /// Push-scheme streaming (alternative; pull scheme preferred for cache).
    pub fn stream_push<T: RealField + Copy>(f_src: &[T], f_dst: &mut [T], nx: usize, ny: usize) {
        // Zero destination first
        for v in f_dst.iter_mut() {
            *v = T::zero();
        }
        for j in 0..ny {
            for i in 0..nx {
                for q in 0..9 {
                    let (ex, ey) = D2Q9::VELOCITIES[q];
                    let dst_i = ((i as i32 + ex).rem_euclid(nx as i32)) as usize;
                    let dst_j = ((j as i32 + ey).rem_euclid(ny as i32)) as usize;
                    f_dst[f_idx(dst_j, dst_i, q, nx)] = f_src[f_idx(j, i, q, nx)];
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_streaming_preserves_mass() {
        // Theorem: streaming is a permutation → ∑ f conserved exactly.
        let nx = 10_usize;
        let ny = 10_usize;
        let n = nx * ny * 9;

        let f_src: Vec<f64> = (0..n).map(|k| (k as f64 % 0.111) * 0.01 + 0.1).collect();
        let mut f_dst = vec![0.0_f64; n];

        let initial_mass: f64 = f_src.iter().sum();
        StreamingOperator::stream(&f_src, &mut f_dst, nx, ny);
        let final_mass: f64 = f_dst.iter().sum();

        assert_relative_eq!(initial_mass, final_mass, epsilon = 1e-10);
    }

    #[test]
    fn test_f_idx_reflects_layout() {
        let nx = 5_usize;
        // Direction 0 of cell (j=2, i=3) must be adjacent to direction 1
        let idx0 = f_idx(2, 3, 0, nx);
        let idx1 = f_idx(2, 3, 1, nx);
        assert_eq!(idx1 - idx0, 1, "q direction must be stride-1");
    }
}
