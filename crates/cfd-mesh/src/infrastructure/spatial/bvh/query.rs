//! Iterative BVH traversal using a fixed-size stack.
//!
//! Avoids recursion and heap allocation during queries.  Node AABBs are read
//! without a token (pure spatial culling data); connectivity is read through
//! the [`GhostToken`] (branded structural data).

use super::node::{BvhNodeKind, MAX_STACK_DEPTH};
use crate::domain::geometry::aabb::Aabb;
use crate::infrastructure::permission::{GhostToken, PermissionedArena};

/// Iterative AABB-overlap query over a flat-arena BVH.
///
/// Traverses with a `[u32; MAX_STACK_DEPTH]` stack.  Exact per-primitive
/// checks at leaves eliminate false positives from conservative union AABBs.
///
/// # Arguments
///
/// - `node_aabbs`  — node bounding boxes (token-free, pure geometry)
/// - `node_kinds`  — node connectivity (requires `token`)
/// - `indices`     — permuted primitive index table
/// - `prim_aabbs`  — source primitive AABBs for exact leaf-level checks
/// - `query`       — the AABB to test against
/// - `token`       — branded permission token matching `node_kinds`
/// - `out`         — accumulator; not cleared before appending
pub(super) fn query_overlapping<'brand>(
    node_aabbs: &[Aabb],
    node_kinds: &PermissionedArena<'brand, BvhNodeKind>,
    indices: &[usize],
    prim_aabbs: &[Aabb],
    query: &Aabb,
    token: &GhostToken<'brand>,
    out: &mut Vec<usize>,
) {
    if node_aabbs.is_empty() {
        return;
    }

    let mut stack = [0u32; MAX_STACK_DEPTH];
    let mut top: usize = 0;
    stack[top] = 0;
    top += 1;

    while top > 0 {
        top -= 1;
        let node_idx = stack[top] as usize;

        // Token-free AABB cull — ends with no borrow held.
        if !node_aabbs[node_idx].intersects(query) {
            continue;
        }

        // `Copy` clone is a 12-byte register copy; the arena borrow ends
        // before we access `indices` or `prim_aabbs` below.
        let kind = *node_kinds.get(node_idx, token);

        match kind {
            BvhNodeKind::Leaf { start, end } => {
                for i in (start as usize)..(end as usize) {
                    let prim_idx = indices[i];
                    // Exact per-primitive check; eliminates false positives
                    // from the conservative union AABB of the leaf node.
                    if prim_aabbs[prim_idx].intersects(query) {
                        out.push(prim_idx);
                    }
                }
            }
            BvhNodeKind::Inner { left, right } => {
                // Push right first so left is popped first (DFS order).
                stack[top] = right;
                top += 1;
                stack[top] = left;
                top += 1;
            }
        }
    }
}
