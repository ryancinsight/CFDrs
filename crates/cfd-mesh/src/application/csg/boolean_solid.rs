//! Dependency-inversion trait for CSG Boolean operations on closed solid meshes.
//!
//! ## Design
//!
//! [`BooleanSolid`] decouples call sites from the concrete [`csg_boolean_indexed`]
//! implementation via the Dependency Inversion Principle (DIP).  Any future
//! alternative CSG backend (e.g., CGAL bindings, GPU-accelerated intersection)
//! can be substituted by implementing this trait without changing call sites.
//!
//! ## Implemented For
//!
//! - [`IndexedMesh`] — delegates to [`csg_boolean_indexed`].
//!
//! ## Backward Compatibility
//!
//! The free functions [`csg_boolean`] and [`csg_boolean_indexed`] are **not removed**;
//! existing call sites continue to compile unchanged.

use crate::application::csg::boolean::{csg_boolean, BooleanOp};
use crate::domain::core::error::MeshResult;
use crate::domain::mesh::IndexedMesh;

// ── Trait ────────────────────────────────────────────────────────────────────

/// DIP-compliant interface for CSG Boolean operations on closed solid meshes.
///
/// Implementors must produce watertight output for watertight inputs (the
/// "watertightness guarantee").
///
/// # Example
///
/// ```ignore
/// use cfd_mesh::application::csg::BooleanSolid;
///
/// let union = cube.union(&sphere)?;
/// let diff  = cube.difference(&sphere)?;
/// let inter = cube.intersection(&sphere)?;
/// ```
pub trait BooleanSolid: Sized {
    /// Compute `self ∪ other`.
    fn union(&self, other: &Self) -> MeshResult<Self>;

    /// Compute `self ∩ other`.
    fn intersection(&self, other: &Self) -> MeshResult<Self>;

    /// Compute `self \ other` (self minus other).
    fn difference(&self, other: &Self) -> MeshResult<Self>;
}

// ── Implementation for IndexedMesh ────────────────────────────────────────────

impl BooleanSolid for IndexedMesh {
    fn union(&self, other: &Self) -> MeshResult<Self> {
        csg_boolean(BooleanOp::Union, self, other)
    }

    fn intersection(&self, other: &Self) -> MeshResult<Self> {
        csg_boolean(BooleanOp::Intersection, self, other)
    }

    fn difference(&self, other: &Self) -> MeshResult<Self> {
        csg_boolean(BooleanOp::Difference, self, other)
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, PrimitiveMesh};

    /// Axis-aligned unit cube centred at the origin.
    fn cube_a() -> IndexedMesh {
        Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_a build")
    }

    /// Offset cube — overlaps cube_a by half its width along X.
    fn cube_b() -> IndexedMesh {
        Cube {
            origin: Point3r::new(-0.5, -0.5, -0.5),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube_b build")
    }

    /// `BooleanSolid::union` produces a valid result for overlapping cubes.
    #[test]
    fn boolean_solid_union_succeeds() {
        assert!(cube_a().union(&cube_b()).is_ok(), "union should succeed");
    }

    /// `BooleanSolid::intersection` produces a valid result for overlapping cubes.
    #[test]
    fn boolean_solid_intersection_succeeds() {
        assert!(
            cube_a().intersection(&cube_b()).is_ok(),
            "intersection should succeed"
        );
    }

    /// `BooleanSolid::difference` produces a valid result for overlapping cubes.
    #[test]
    fn boolean_solid_difference_succeeds() {
        assert!(
            cube_a().difference(&cube_b()).is_ok(),
            "difference should succeed"
        );
    }

    /// Trait and free function agree on result type for union.
    #[test]
    fn boolean_solid_matches_free_fn_union() {
        let a = cube_a();
        let b = cube_b();
        let via_trait = a.union(&b);
        let via_fn = csg_boolean(BooleanOp::Union, &a, &b);
        assert_eq!(via_trait.is_ok(), via_fn.is_ok());
    }
}
