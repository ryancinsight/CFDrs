//! BVH node-kind discriminant and build/traversal constants.

/// Maximum primitives per leaf node before forced splitting.
pub(super) const MAX_LEAF_PRIMITIVES: usize = 4;

/// SAH traversal cost relative to primitive intersection cost.
pub(super) const SAH_TRAVERSAL_COST: f64 = 1.2;

/// Maximum traversal stack depth for the iterative query loop.
///
/// log₂(2^63) = 63; 64 guards against any conceivable input.  The forced
/// median-split fallback in `build_recursive` ensures depth ≤
/// ⌈log₂(N / MAX_LEAF_PRIMITIVES)⌉ + 1.
pub(super) const MAX_STACK_DEPTH: usize = 64;

/// Discriminates leaf nodes from inner nodes in the BVH arena.
///
/// Stored in a [`crate::infrastructure::permission::PermissionedArena`] and
/// accessed via [`crate::infrastructure::permission::GhostToken`].  `Copy`
/// makes cloning a cheap 12-byte register copy during traversal — the arena
/// borrow ends before `indices` or `prim_aabbs` are accessed.
#[derive(Clone, Copy, Debug)]
pub(crate) enum BvhNodeKind {
    /// Inner node: arena indices of left and right child nodes.
    Inner { left: u32, right: u32 },
    /// Leaf node: half-open range `[start, end)` into the permuted index vec.
    Leaf { start: u32, end: u32 },
}
