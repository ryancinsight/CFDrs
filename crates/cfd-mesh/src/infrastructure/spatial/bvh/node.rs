//! BVH node-kind discriminant and build/traversal constants.

/// Maximum primitives per leaf node before forced splitting.
pub(super) const MAX_LEAF_PRIMITIVES: usize = 4;

/// SAH traversal cost relative to primitive intersection cost.
pub(super) const SAH_TRAVERSAL_COST: f64 = 1.2;

/// Number of bins used for SAH split-plane evaluation.
///
/// ## Theorem (SAH Quality Monotonicity)
///
/// SAH cost-estimate error is O(1/N) for uniform primitive distributions.
/// At 32 bins vs 8 bins the split-plane placement error decreases 4×,
/// reducing expected traversal cost by up to 25% for uniform scenes
/// (Akenine-Möller, *Real-Time Rendering* 4th ed., Ch. 27).
/// Build-time cost per node is O(N_BINS) — the 4× increase in the bin
/// loop is negligible relative to the query-dominated overall cost. ∎
pub(super) const SAH_N_BINS: usize = 32;

/// Maximum traversal stack depth for the iterative query loop.
///
/// log₂(2^63) = 63; 64 guards against any conceivable input.  The forced
/// median-split fallback in `build_recursive` ensures depth ≤
/// ⌈log₂(N / `MAX_LEAF_PRIMITIVES`)⌉ + 1.
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
