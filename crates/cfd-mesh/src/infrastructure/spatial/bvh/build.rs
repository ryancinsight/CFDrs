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
    indices: &mut Vec<usize>,
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

    let (split_axis, split_pos) = sah_split(aabbs, indices, start, end, &node_aabb);
    let mid = partition(aabbs, indices, start, end, split_axis, split_pos);

    // Forced median split when SAH partition degenerates (all on one side).
    let mid = if mid == start || mid == end {
        start + count / 2
    } else {
        mid
    };

    // Reserve parent slot; children are appended after this index.
    let parent_idx = out_aabbs.len() as u32;
    out_aabbs.push(node_aabb);
    out_kinds.push(BvhNodeKind::Inner { left: 0, right: 0 }); // patched below

    let left_idx = build_recursive(aabbs, indices, start, mid, out_aabbs, out_kinds);
    let right_idx = build_recursive(aabbs, indices, mid, end, out_aabbs, out_kinds);

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
        let c = aabbs[idx].center();
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
            let c = axis_value(&aabbs[idx].center(), axis);
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
            let l_cnt = prefix_cnt[split_k] as f64;
            let r_sa = surface_area(&right_aabb);
            let r_cnt = right_count as f64;

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

/// Stable partition of `indices[start..end]` around `split_value` on `axis`.
///
/// Returns `mid` such that centroids `< split_value` occupy `[start, mid)`.
/// Uses scratch buffers to avoid `usize` underflow in an in-place two-pointer
/// approach.
pub(super) fn partition(
    aabbs: &[Aabb],
    indices: &mut Vec<usize>,
    start: usize,
    end: usize,
    axis: usize,
    split_value: f64,
) -> usize {
    let slice = &mut indices[start..end];
    let n = slice.len();
    let mut left_buf: Vec<usize> = Vec::with_capacity(n);
    let mut right_buf: Vec<usize> = Vec::new();

    for &idx in slice.iter() {
        if axis_value(&aabbs[idx].center(), axis) < split_value {
            left_buf.push(idx);
        } else {
            right_buf.push(idx);
        }
    }
    let mid = start + left_buf.len();
    left_buf.extend(right_buf);
    slice.copy_from_slice(&left_buf);
    mid
}
