//! SAH-BVH construction: recursive build, bin-based SAH split, partition.
//!
//! # Build strategy
//!
//! 1. Recurse into plain `Vec<Aabb>` / `Vec<BvhNodeKind>` (no token required).
//! 2. Transfer connectivity into a `PermissionedArena` in one pass after build.
//!
//! This two-phase approach decouples the mutable build logic from the branded
//! permission system, keeping both sides clean.

use super::geo::{axis_extent, axis_min, axis_value, longest_axis, surface_area};
use super::node::{BvhNodeKind, MAX_LEAF_PRIMITIVES, SAH_TRAVERSAL_COST};
use crate::domain::geometry::aabb::Aabb;
use nalgebra::Point3;
use std::cmp::Ordering;

// ── Public build entry ────────────────────────────────────────────────────────

/// Compute the union AABB of `indices[start..end]`.
pub(super) fn range_aabb(aabbs: &[Aabb], indices: &[usize], start: usize, end: usize) -> Aabb {
    let mut a = Aabb::empty();
    for &idx in &indices[start..end] {
        a = a.union(&aabbs[idx]);
    }
    a
}

/// Recursively build SAH-BVH into parallel plain `Vec`s.
///
/// Appends to `out_aabbs` and `out_kinds` in sync and returns the index of the
/// root node just pushed.
///
/// # Panics (debug only)
/// Asserts `start < end`.
pub(super) fn build_recursive(
    aabbs: &[Aabb],
    centroids: &[Point3<f64>],
    indices: &mut [usize],
    start: usize,
    end: usize,
    out_aabbs: &mut Vec<Aabb>,
    out_kinds: &mut Vec<BvhNodeKind>,
) -> u32 {
    debug_assert!(start < end, "empty range in build_recursive");

    let node_aabb = range_aabb(aabbs, indices, start, end);
    let count = end - start;

    // Leaf threshold.
    if count <= MAX_LEAF_PRIMITIVES {
        let idx = out_aabbs.len() as u32;
        out_aabbs.push(node_aabb);
        out_kinds.push(BvhNodeKind::Leaf {
            start: start as u32,
            end: end as u32,
        });
        return idx;
    }

    let (split_axis, split_pos) = sah_split(aabbs, centroids, indices, start, end, &node_aabb);
    let mut mid = partition(centroids, indices, start, end, split_axis, split_pos);

    // Forced median split when SAH partition degenerates (all on one side).
    if mid == start || mid == end {
        mid = partition_median(centroids, indices, start, end, split_axis);
    }

    // Reserve parent slot; children are appended after this index.
    let parent_idx = out_aabbs.len() as u32;
    out_aabbs.push(node_aabb);
    out_kinds.push(BvhNodeKind::Inner { left: 0, right: 0 }); // patched below

    let left_idx = build_recursive(aabbs, centroids, indices, start, mid, out_aabbs, out_kinds);
    let right_idx = build_recursive(aabbs, centroids, indices, mid, end, out_aabbs, out_kinds);

    // Patch parent connectivity.  `parent_idx` is valid: it was pushed before
    // the recursive calls that only append beyond it.
    out_kinds[parent_idx as usize] = BvhNodeKind::Inner {
        left: left_idx,
        right: right_idx,
    };
    parent_idx
}

// ── SAH split ─────────────────────────────────────────────────────────────────

/// Choose the best splitting axis and position via 8-bin SAH.
///
/// Evaluates 7 split planes per axis (8 bins, 3 axes = 21 candidates) and
/// returns `(axis, centroid_split_value)`.
fn sah_split(
    aabbs: &[Aabb],
    centroids: &[Point3<f64>],
    indices: &[usize],
    start: usize,
    end: usize,
    parent_aabb: &Aabb,
) -> (usize, f64) {
    const N_BINS: usize = 8;

    let parent_sa = surface_area(parent_aabb);
    let count = (end - start) as f64;

    // Centroid AABB determines bin placement.
    let mut centroid_aabb = Aabb::empty();
    for &idx in &indices[start..end] {
        let c = centroids[idx];
        centroid_aabb.expand(&c);
    }

    let mut best_cost = f64::INFINITY;
    let mut best_axis = 0usize;
    let mut best_split = 0.0f64;

    for axis in 0..3usize {
        let extent = axis_extent(&centroid_aabb, axis);
        if extent < 1e-15 {
            continue;
        } // degenerate – all centroids coplanar

        let inv_extent = 1.0 / extent;
        let min_c = axis_min(&centroid_aabb, axis);

        let mut bin_aabb = [Aabb::empty(); N_BINS];
        let mut bin_count = [0u32; N_BINS];

        for &idx in &indices[start..end] {
            let c = axis_value(&centroids[idx], axis);
            let b = ((c - min_c) * inv_extent * N_BINS as f64) as usize;
            let b = b.min(N_BINS - 1);
            bin_aabb[b] = bin_aabb[b].union(&aabbs[idx]);
            bin_count[b] += 1;
        }

        // Left-prefix scan.
        let mut left_aabb = Aabb::empty();
        let mut left_count = 0u32;
        let mut prefix_sa = [0.0f64; N_BINS];
        let mut prefix_cnt = [0u32; N_BINS];

        for k in 0..(N_BINS - 1) {
            left_aabb = left_aabb.union(&bin_aabb[k]);
            left_count += bin_count[k];
            prefix_sa[k] = surface_area(&left_aabb);
            prefix_cnt[k] = left_count;
        }

        // Right-suffix scan.
        let mut right_aabb = Aabb::empty();
        let mut right_count = 0u32;

        for k in (1..N_BINS).rev() {
            right_aabb = right_aabb.union(&bin_aabb[k]);
            right_count += bin_count[k];

            let split_k = k - 1;
            let l_sa = prefix_sa[split_k];
            let l_cnt = f64::from(prefix_cnt[split_k]);
            let r_sa = surface_area(&right_aabb);
            let r_cnt = f64::from(right_count);

            if l_cnt == 0.0 || r_cnt == 0.0 {
                continue;
            }

            let cost = SAH_TRAVERSAL_COST + (l_sa / parent_sa) * l_cnt + (r_sa / parent_sa) * r_cnt;

            if cost < best_cost {
                best_cost = cost;
                best_axis = axis;
                best_split = min_c + (split_k as f64 + 1.0) / N_BINS as f64 * extent;
            }
        }
    }

    // Fallback: longest-axis centroid midpoint.
    if best_cost == f64::INFINITY || best_cost >= count {
        let (axis, _) = longest_axis(parent_aabb);
        return (axis, axis_value(&parent_aabb.center(), axis));
    }

    (best_axis, best_split)
}

// ── Partition ─────────────────────────────────────────────────────────────────

/// In-place partition of `indices[start..end]` around `split_value` on `axis`.
///
/// Returns `mid` such that centroids `< split_value` occupy `[start, mid)`.
/// Uses an in-place two-pointer walk with no auxiliary allocation.
///
/// # Theorem — Partition Soundness and Completeness
///
/// For the return value `mid`, every index in `[start, mid)` satisfies
/// `centroid_axis < split_value`, and every index in `[mid, end)` satisfies
/// `centroid_axis >= split_value`.
///
/// **Proof sketch**: The left pointer only advances past values that satisfy
/// the left predicate; the right pointer only retreats past values that satisfy
/// the right predicate. Whenever both pointers stop with `left < right`, the
/// two offending elements are swapped, restoring predicate validity on both
/// sides before pointers continue. Termination occurs when `left == right`,
/// at which point no violating pair remains. Since only swaps are used, the
/// output is a permutation of the input range. ∎
pub(super) fn partition(
    centroids: &[Point3<f64>],
    indices: &mut [usize],
    start: usize,
    end: usize,
    axis: usize,
    split_value: f64,
) -> usize {
    let mut left = start;
    let mut right = end;

    while left < right {
        while left < right && axis_value(&centroids[indices[left]], axis) < split_value {
            left += 1;
        }
        while left < right && axis_value(&centroids[indices[right - 1]], axis) >= split_value {
            right -= 1;
        }
        if left < right {
            indices.swap(left, right - 1);
            left += 1;
            right -= 1;
        }
    }

    left
}

/// Deterministic median partition fallback for degenerate SAH splits.
///
/// Uses `select_nth_unstable_by` on centroid coordinates along `axis`.
///
/// # Theorem — Balanced fallback partition
///
/// Let `n = end - start` and `k = n / 2`. After selecting nth on
/// `indices[start..end]`, exactly `k` items occupy `[start, start + k)`, and
/// `n-k` items occupy `[start + k, end)`. Therefore both recursive ranges are
/// non-empty for `n >= 2`, and their cardinalities differ by at most 1. ∎
pub(super) fn partition_median(
    centroids: &[Point3<f64>],
    indices: &mut [usize],
    start: usize,
    end: usize,
    axis: usize,
) -> usize {
    debug_assert!(start < end, "empty range in partition_median");
    let local_mid = (end - start) / 2;
    indices[start..end].select_nth_unstable_by(local_mid, |lhs, rhs| {
        axis_value(&centroids[*lhs], axis)
            .partial_cmp(&axis_value(&centroids[*rhs], axis))
            .unwrap_or(Ordering::Equal)
    });
    start + local_mid
}

/// Build a centroid cache aligned with `aabbs` by index.
///
/// # Theorem — Centroid cache equivalence
///
/// If `centroids[i] = aabbs[i].center()` for every `i`, then any decision that
/// depends only on centroid axis coordinates is identical whether computed from
/// `aabbs[i].center()` on demand or read from `centroids[i]`.
///
/// **Proof sketch**: each lookup returns the same real value by construction,
/// so all axis comparisons and bucket index computations are unchanged. ∎
#[must_use]
pub(super) fn build_centroids(aabbs: &[Aabb]) -> Vec<Point3<f64>> {
    aabbs.iter().map(Aabb::center).collect()
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;
    use proptest::prelude::*;
    use std::collections::BTreeSet;

    fn pt(x: f64, y: f64, z: f64) -> Point3<f64> {
        Point3::new(x, y, z)
    }

    fn point_aabb(x: f64, y: f64, z: f64) -> Aabb {
        Aabb::new(pt(x, y, z), pt(x, y, z))
    }

    proptest! {
        #[test]
        fn partition_predicate_and_permutation_hold(
            vals in prop::collection::vec(-50_i16..50_i16, 2..128),
            split in -50_i16..50_i16,
        ) {
            let aabbs: Vec<Aabb> = vals
                .iter()
                .map(|&x| point_aabb(f64::from(x), 0.0, 0.0))
                .collect();
            let centroids = build_centroids(&aabbs);
            let mut indices: Vec<usize> = (0..aabbs.len()).collect();
            let before: BTreeSet<usize> = indices.iter().copied().collect();
            let end = indices.len();
            let mid = partition(&centroids, &mut indices, 0, end, 0, f64::from(split));

            for &idx in &indices[..mid] {
                prop_assert!(axis_value(&aabbs[idx].center(), 0) < f64::from(split));
            }
            for &idx in &indices[mid..] {
                prop_assert!(axis_value(&aabbs[idx].center(), 0) >= f64::from(split));
            }

            let after: BTreeSet<usize> = indices.iter().copied().collect();
            prop_assert_eq!(before, after);
        }
    }

    #[test]
    fn median_partition_balances_degenerate_centroids() {
        let aabbs = vec![point_aabb(1.0, 0.0, 0.0); 11];
        let centroids = build_centroids(&aabbs);
        let mut indices: Vec<usize> = (0..aabbs.len()).collect();
        let mid = partition_median(&centroids, &mut indices, 0, aabbs.len(), 0);
        assert_eq!(mid, 5);
    }

    #[test]
    fn median_partition_keeps_all_indices() {
        let aabbs = vec![
            point_aabb(4.0, 0.0, 0.0),
            point_aabb(1.0, 0.0, 0.0),
            point_aabb(9.0, 0.0, 0.0),
            point_aabb(2.0, 0.0, 0.0),
            point_aabb(7.0, 0.0, 0.0),
            point_aabb(3.0, 0.0, 0.0),
        ];
        let centroids = build_centroids(&aabbs);
        let mut indices = vec![5, 2, 4, 0, 1, 3];
        let before: BTreeSet<usize> = indices.iter().copied().collect();
        let mid = partition_median(&centroids, &mut indices, 0, 6, 0);
        assert_eq!(mid, 3);
        let after: BTreeSet<usize> = indices.iter().copied().collect();
        assert_eq!(before, after);
    }

    #[test]
    fn centroid_cache_matches_direct_centers() {
        let aabbs = vec![
            point_aabb(1.0, 2.0, 3.0),
            point_aabb(-1.0, 0.0, 9.0),
            point_aabb(4.5, -7.0, 0.25),
        ];
        let centroids = build_centroids(&aabbs);
        assert_eq!(centroids.len(), aabbs.len());
        for (c, a) in centroids.iter().zip(aabbs.iter()) {
            assert_eq!(*c, a.center());
        }
    }
}
