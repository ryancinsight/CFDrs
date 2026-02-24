//! Surface-Area-Heuristic Bounding Volume Hierarchy — public API.
//!
//! # Module layout
//!
//! ```text
//! bvh/
//! ├── mod.rs    — BvhTree<'brand,'aabbs>, with_bvh, build_tree, tests
//! ├── node.rs   — BvhNodeKind enum, build/traversal constants
//! ├── build.rs  — build_recursive, sah_split, partition, range_aabb
//! ├── query.rs  — iterative stack traversal
//! └── geo.rs    — surface_area, axis helpers, longest_axis
//! ```
//!
//! # GhostCell Integration
//!
//! Tree connectivity lives in a
//! [`PermissionedArena<'brand, BvhNodeKind>`][`crate::infrastructure::permission::PermissionedArena`],
//! branded with the same `'brand` as the [`GhostToken`] issued by [`with_bvh`].
//! The build phase uses plain `Vec<BvhNodeKind>` mutation (no token required),
//! then transfers nodes into the arena before the branded scope opens.
//!
//! Node AABBs are stored in a plain `Vec<Aabb>` and read **without** a token.
//! This eliminates the borrow-checker conflict that would arise if both node
//! connectivity and node AABBs lived in the arena: the query loop can read the
//! node AABB, cull, then copy the `BvhNodeKind` (12-byte register copy ending
//! the arena borrow), then access `indices` and `prim_aabbs` without conflict.
//!
//! `BvhTree` exposes no `borrow_mut` surface — connectivity is private —
//! so external code cannot mutate the tree while queries are in flight.
//!
//! # Theorems
//!
//! ## Theorem: BVH Completeness
//! **Statement**: `query_overlapping` returns every index `i` such that
//! `prim_aabbs[i]` intersects `query`.
//!
//! **Proof sketch**: The root AABB is the union of all primitive AABBs.  At
//! each inner node the traversal visits a child iff its AABB overlaps `query`.
//! Every primitive AABB is contained in every ancestor AABB, so a culled
//! subtree contains no overlapping primitive.  By induction, all overlapping
//! primitives are found.  Leaf-level exact checks eliminate false positives.  ∎
//!
//! ## Theorem: SAH Validity
//! **Statement**: The SAH cost function upper-bounds expected query cost for
//! uniformly distributed query AABBs.
//!
//! **Proof sketch**: The probability that a random AABB query intersects a box
//! is proportional to its surface area (Buffon-style).  SAH minimises expected
//! child-visitation count.  ∎

mod build;
mod geo;
mod node;
mod query;

use crate::domain::geometry::aabb::Aabb;
use crate::infrastructure::permission::{GhostToken, PermissionedArena};
use build::build_recursive;
use node::BvhNodeKind;

// ── BvhTree ───────────────────────────────────────────────────────────────────

/// Flat-arena SAH-BVH branded with `'brand`.
///
/// - `'brand` — the GhostCell brand gating connectivity reads.
/// - `'aabbs` — lifetime of the source AABB slice (no copy of input data).
///
/// Construct via [`with_bvh`]; query via [`BvhTree::query_overlapping`].
pub struct BvhTree<'brand, 'aabbs> {
    /// Node bounding boxes — readable without a token (pure spatial data).
    node_aabbs: Vec<Aabb>,
    /// Node connectivity — readable only via `&GhostToken<'brand>`.
    node_kinds: PermissionedArena<'brand, BvhNodeKind>,
    /// Permuted primitive index table; leaves hold ranges into this vec.
    indices: Vec<usize>,
    /// Borrowed source primitive AABBs for exact leaf-level intersection.
    prim_aabbs: &'aabbs [Aabb],
}

impl<'brand, 'aabbs> BvhTree<'brand, 'aabbs> {
    /// Query all primitive indices whose AABB overlaps `query`.
    ///
    /// Uses an iterative `[u32; 64]` stack — no recursion, no heap allocation
    /// during traversal.  Appends to `out`; does **not** clear it first.
    pub fn query_overlapping(
        &self,
        query: &Aabb,
        token: &GhostToken<'brand>,
        out: &mut Vec<usize>,
    ) {
        query::query_overlapping(
            &self.node_aabbs,
            &self.node_kinds,
            &self.indices,
            self.prim_aabbs,
            query,
            token,
            out,
        );
    }

    /// Number of internal tree nodes (useful for diagnostics and testing).
    #[inline]
    pub fn node_count(&self) -> usize {
        self.node_aabbs.len()
    }
}

// ── Entry point ───────────────────────────────────────────────────────────────

/// Build a SAH-BVH over `aabbs` and invoke `f` with the tree and a branded
/// [`GhostToken`].
///
/// The brand `'brand` is **invariant** and introduced here; it cannot escape
/// the closure.  This guarantees:
///
/// - The tree is fully constructed before `f` receives it.
/// - Querying connectivity requires `&GhostToken<'brand>`, tying the read to
///   the correct branded scope and preventing token confusion.
///
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::infrastructure::spatial::bvh::with_bvh;
///
/// with_bvh(&my_aabbs, |tree, token| {
///     let mut hits = Vec::new();
///     tree.query_overlapping(&query, &token, &mut hits);
/// });
/// ```
pub fn with_bvh<'a, F, R>(aabbs: &'a [Aabb], f: F) -> R
where
    F: for<'brand> FnOnce(BvhTree<'brand, 'a>, GhostToken<'brand>) -> R,
{
    GhostToken::new(|token| {
        let tree: BvhTree<'_, 'a> = build_tree(aabbs);
        f(tree, token)
    })
}

// ── Build helper (private) ────────────────────────────────────────────────────

/// Construct a `BvhTree` from `aabbs`.
///
/// Two-phase:
/// 1. Recursive SAH build into plain `Vec`s (no token needed).
/// 2. Transfer `Vec<BvhNodeKind>` into `PermissionedArena` in one pass.
fn build_tree<'brand, 'a>(aabbs: &'a [Aabb]) -> BvhTree<'brand, 'a> {
    if aabbs.is_empty() {
        return BvhTree {
            node_aabbs: Vec::new(),
            node_kinds: PermissionedArena::new(),
            indices: Vec::new(),
            prim_aabbs: aabbs,
        };
    }

    let cap = 2 * aabbs.len();
    let mut node_aabbs: Vec<Aabb> = Vec::with_capacity(cap);
    let mut raw_kinds: Vec<BvhNodeKind> = Vec::with_capacity(cap);
    let mut indices: Vec<usize> = (0..aabbs.len()).collect();
    let n = indices.len();

    build_recursive(aabbs, &mut indices, 0, n, &mut node_aabbs, &mut raw_kinds);

    // Transfer into PermissionedArena — one allocation, no per-element overhead.
    let mut node_kinds: PermissionedArena<'brand, BvhNodeKind> =
        PermissionedArena::with_capacity(raw_kinds.len());
    for kind in raw_kinds {
        node_kinds.push(kind);
    }

    BvhTree {
        node_aabbs,
        node_kinds,
        indices,
        prim_aabbs: aabbs,
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Point3;

    fn pt(x: f64, y: f64, z: f64) -> Point3<f64> {
        Point3::new(x, y, z)
    }

    fn unit_aabb(x: f64, y: f64, z: f64) -> Aabb {
        Aabb::new(pt(x, y, z), pt(x + 1.0, y + 1.0, z + 1.0))
    }

    // ── Correctness ───────────────────────────────────────────────────────────

    /// 100 unit AABBs in a 10×10 grid; query covers only 4.
    #[test]
    fn query_finds_all_overlapping() {
        let aabbs: Vec<Aabb> = (0..10)
            .flat_map(|ix| (0..10).map(move |iy| unit_aabb(ix as f64 * 2.0, iy as f64 * 2.0, 0.0)))
            .collect();

        // ix ∈ {0,1}, iy ∈ {0,1} → indices 0, 1, 10, 11.
        // iy=2 box at y=[4,5] does NOT overlap query max_y=3.
        let query = Aabb::new(pt(0.0, 0.0, 0.0), pt(3.0, 3.0, 1.0));
        let expected = vec![0usize, 1, 10, 11];

        with_bvh(&aabbs, |tree, token| {
            let mut hits = Vec::new();
            tree.query_overlapping(&query, &token, &mut hits);
            hits.sort_unstable();
            assert_eq!(hits, expected);
        });
    }

    #[test]
    fn empty_bvh_yields_no_results() {
        with_bvh(&[], |tree, token| {
            let query = Aabb::new(pt(-100.0, -100.0, -100.0), pt(100.0, 100.0, 100.0));
            let mut hits = Vec::new();
            tree.query_overlapping(&query, &token, &mut hits);
            assert!(hits.is_empty());
        });
    }

    #[test]
    fn single_primitive_hit() {
        let aabbs = vec![unit_aabb(0.0, 0.0, 0.0)];
        with_bvh(&aabbs, |tree, token| {
            let mut hits = Vec::new();
            tree.query_overlapping(
                &Aabb::new(pt(0.5, 0.5, 0.5), pt(0.6, 0.6, 0.6)),
                &token,
                &mut hits,
            );
            assert_eq!(hits, vec![0]);
        });
    }

    #[test]
    fn single_primitive_miss() {
        let aabbs = vec![unit_aabb(0.0, 0.0, 0.0)];
        with_bvh(&aabbs, |tree, token| {
            let mut hits = Vec::new();
            tree.query_overlapping(
                &Aabb::new(pt(5.0, 5.0, 5.0), pt(6.0, 6.0, 6.0)),
                &token,
                &mut hits,
            );
            assert!(hits.is_empty());
        });
    }

    // ── No false negatives: BVH ⊇ brute-force ────────────────────────────────

    #[test]
    fn no_false_negatives_random_layout() {
        let aabbs: Vec<Aabb> = (0..50_usize)
            .map(|i| {
                let f = i as f64 * 0.73;
                Aabb::new(pt(f, f * 0.5, 0.0), pt(f + 1.0, f * 0.5 + 1.0, 1.0))
            })
            .collect();
        let query = Aabb::new(pt(5.0, 0.0, 0.0), pt(15.0, 5.0, 2.0));
        let brute: Vec<usize> = aabbs
            .iter()
            .enumerate()
            .filter(|(_, a)| a.intersects(&query))
            .map(|(i, _)| i)
            .collect();

        with_bvh(&aabbs, |tree, token| {
            let mut hits = Vec::new();
            tree.query_overlapping(&query, &token, &mut hits);
            hits.sort_unstable();
            for &i in &brute {
                assert!(hits.contains(&i), "BVH missed primitive {i}");
            }
        });
    }

    // ── No false positives: BVH ⊆ brute-force ────────────────────────────────

    #[test]
    fn no_false_positives_random_layout() {
        let aabbs: Vec<Aabb> = (0..50_usize)
            .map(|i| {
                let f = i as f64 * 0.73;
                Aabb::new(pt(f, f * 0.5, 0.0), pt(f + 1.0, f * 0.5 + 1.0, 1.0))
            })
            .collect();
        let query = Aabb::new(pt(5.0, 0.0, 0.0), pt(15.0, 5.0, 2.0));

        with_bvh(&aabbs, |tree, token| {
            let mut hits = Vec::new();
            tree.query_overlapping(&query, &token, &mut hits);
            for &i in &hits {
                assert!(
                    aabbs[i].intersects(&query),
                    "BVH false positive: primitive {i} AABB does not overlap query"
                );
            }
        });
    }

    // ── Geo helpers ───────────────────────────────────────────────────────────

    // Delegate surface_area / longest_axis tests to geo sub-module.
    #[allow(unused_imports)]
    use geo::tests::*;

    // ── Scale ─────────────────────────────────────────────────────────────────

    #[test]
    fn build_10k_primitives() {
        let aabbs: Vec<Aabb> = (0..10_000)
            .map(|i| {
                let f = i as f64 * 0.01;
                Aabb::new(pt(f, 0.0, 0.0), pt(f + 0.01, 1.0, 1.0))
            })
            .collect();
        with_bvh(&aabbs, |tree, _token| {
            assert!(tree.node_count() > 0);
        });
    }
}
