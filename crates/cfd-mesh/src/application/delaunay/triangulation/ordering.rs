//! Hilbert-curve vertex insertion ordering for optimal Delaunay performance.
//!
//! # Theorem — Hilbert-Curve Insertion Optimality
//!
//! **Statement**: By sorting 2D points along a Hilbert space-filling curve
//! before incremental insertion, the expected Lawson walk distance per
//! insertion is reduced from $O(\sqrt{n})$ to near-$O(1)$, yielding an
//! overall expected construction time of $O(n \log n)$.
//!
//! **Proof sketch**: The Hilbert curve preserves 2D locality — nearby points
//! on the curve are nearby in 2D.  Since the Lawson walk starts from the
//! last-inserted triangle (the walk hint), consecutive spatially-adjacent
//! insertions traverse only $O(1)$ triangles on average.  The $\log n$
//! factor arises from the curve's recursive quadrant structure, which
//! ensures the maximum distance between consecutive insertions scales
//! logarithmically with the number of points at each recursion level.
//!
//! # Reference
//!
//! - Liu, Y. & Yan, D. (2013). "Fast Delaunay triangulation via
//!   kd-tree parent walk hints."
//! - Amenta, N., Choi, S., Rote, G. (2003). "Incremental constructions
//!   with biased randomized insertion order (BRIO)."

use crate::application::delaunay::pslg::vertex::PslgVertex;
use crate::domain::core::scalar::Real;

/// Compute a Hilbert-curve-based insertion order for a set of 2D vertices.
///
/// 1. Maps each vertex to integer coordinates in a $2^{16} \times 2^{16}$ grid.
/// 2. Computes the Hilbert index for each grid cell.
/// 3. Sorts vertices by Hilbert index; returns the sorted permutation.
pub(crate) fn hilbert_order(vertices: &[PslgVertex]) -> Vec<usize> {
    if vertices.is_empty() {
        return Vec::new();
    }

    // Compute bounding box.
    let (mut min_x, mut min_y) = (Real::INFINITY, Real::INFINITY);
    let (mut max_x, mut max_y) = (Real::NEG_INFINITY, Real::NEG_INFINITY);
    for v in vertices {
        if v.x < min_x {
            min_x = v.x;
        }
        if v.y < min_y {
            min_y = v.y;
        }
        if v.x > max_x {
            max_x = v.x;
        }
        if v.y > max_y {
            max_y = v.y;
        }
    }

    let range_x = (max_x - min_x).max(1e-30);
    let range_y = (max_y - min_y).max(1e-30);
    let scale = (1u32 << 16) - 1;

    // Compute Hilbert index for each vertex.
    let mut indexed: Vec<(u64, usize)> = vertices
        .iter()
        .enumerate()
        .map(|(i, v)| {
            let ix = (((v.x - min_x) / range_x) * f64::from(scale)).round() as u32;
            let iy = (((v.y - min_y) / range_y) * f64::from(scale)).round() as u32;
            let hilbert_idx = xy_to_hilbert(ix.min(scale), iy.min(scale), 16);
            (hilbert_idx, i)
        })
        .collect();

    indexed.sort_unstable_by_key(|&(hilbert_idx, _)| hilbert_idx);
    indexed.into_iter().map(|(_, i)| i).collect()
}

/// Convert 2D grid coordinates to a Hilbert curve index.
///
/// Uses the iterative algorithm for an `order`-bit Hilbert curve
/// (grid size $2^{\text{order}} \times 2^{\text{order}}$).
///
/// # Reference
///
/// Wikipedia: Hilbert curve mapping algorithms.
fn xy_to_hilbert(gx: u32, gy: u32, order: u32) -> u64 {
    let mut rot_x: u32;
    let mut rot_y: u32;
    let mut dist: u64 = 0;
    let grid = 1u32 << order;
    let mut tx = gx;
    let mut ty = gy;

    let mut half = grid >> 1;
    while half > 0 {
        rot_x = if (tx & half) > 0 { 1 } else { 0 };
        rot_y = if (ty & half) > 0 { 1 } else { 0 };
        dist += u64::from(half) * u64::from(half) * u64::from((3 * rot_x) ^ rot_y);
        // Rotate quadrant.
        if rot_y == 0 {
            if rot_x == 1 {
                tx = half.wrapping_mul(2).wrapping_sub(1).wrapping_sub(tx);
                ty = half.wrapping_mul(2).wrapping_sub(1).wrapping_sub(ty);
            }
            std::mem::swap(&mut tx, &mut ty);
        }
        half >>= 1;
    }
    dist
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn hilbert_deterministic() {
        let verts = vec![
            PslgVertex::new(0.0, 0.0),
            PslgVertex::new(1.0, 0.0),
            PslgVertex::new(0.5, 1.0),
            PslgVertex::new(0.5, 0.3),
        ];
        let order1 = hilbert_order(&verts);
        let order2 = hilbert_order(&verts);
        assert_eq!(order1, order2);
    }

    #[test]
    fn hilbert_empty() {
        assert!(hilbert_order(&[]).is_empty());
    }

    #[test]
    fn hilbert_single() {
        let order = hilbert_order(&[PslgVertex::new(42.0, 7.0)]);
        assert_eq!(order, vec![0]);
    }

    #[test]
    fn hilbert_preserves_all_indices() {
        let verts: Vec<PslgVertex> = (0..20)
            .map(|i| PslgVertex::new(f64::from(i) * 0.1, f64::from(i) * 0.05))
            .collect();
        let order = hilbert_order(&verts);
        assert_eq!(order.len(), 20);
        let mut sorted = order.clone();
        sorted.sort_unstable();
        let expected: Vec<usize> = (0..20).collect();
        assert_eq!(sorted, expected);
    }
}
