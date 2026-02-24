//! # Half-Edge Topology Kernel
//!
//! The half-edge data structure is the topological backbone of [`crate::domain::mesh::Mesh`].
//! Every undirected edge in the mesh is represented as a pair of directed
//! *half-edges*, one per adjacent face. This gives O(1) local traversal for
//! all CFD stencil operations: face-to-face adjacency, vertex one-ring
//! neighborhoods, and boundary loop iteration.
//!
//! ## Structure
//!
//! ```text
//!          A
//!         / \
//!        /   \
//!      he0   he2      he0: A→B, face=△ABС, twin=he3, next=he1
//!      /  △1  \       he1: B→C, face=△ABC, twin=he_ext, next=he2
//!     B──he1──►C      he2: C→A, face=△ABC, twin=he_ext, next=he0
//!     ◄─he4───C       he3: B→A, face=△BDA, twin=he0, next=he5
//!      \  △2  /       he4: C→B, face=△CBD, twin=he1, next=he3 (approx)
//!      he5   he3      he5: D→C  ...
//!        \   /
//!          D
//! ```
//!
//! ## Invariants
//!
//! All of the following hold for every half-edge `he` in a valid mesh:
//!
//! 1. **Twin involution**: `twin(twin(he)) == he`
//! 2. **Next-prev consistency**: `next(prev(he)) == he`  and  `prev(next(he)) == he`
//! 3. **Face-loop closure**: following `next` from any half-edge returns to `he`
//!    in a finite number of steps
//! 4. **Vertex-ring closure**: following `twin → next` from any half-edge
//!    returns to `he` in a finite number of steps (for interior vertices)
//! 5. **Face consistency**: `face(next(he)) == face(he)` (all half-edges in a
//!    face loop reference the same face)
//!
//! ## Theorem: Twin Involution
//!
//! **Statement**: The `twin` mapping is an involution on the set of half-edges:
//! `twin ∘ twin = id`.
//!
//! **Proof**: Each undirected edge `{u, v}` stores exactly two half-edges:
//! `he(u→v)` and `he(v→u)`. By construction, `twin(he(u→v)) = he(v→u)` and
//! `twin(he(v→u)) = he(u→v)`, so composing twice returns the original. ∎
//!
//! ## Theorem: 2-Manifold Condition
//!
//! **Statement**: A mesh is a combinatorial 2-manifold if and only if:
//! - Every interior edge has exactly 2 incident faces.
//! - Every boundary edge has exactly 1 incident face.
//! - The faces incident to each vertex form a single connected cycle
//!   (topological disk or half-disk for boundary vertices).
//!
//! **Proof sketch**: The local homeomorphism to ℝ² or ℝ²₊ at each point is
//! exactly encoded by these combinatorial conditions. Interior vertices have a
//! disk link; boundary vertices have a half-disk link. ∎

use crate::domain::core::index::{FaceKey, HalfEdgeKey, PatchKey, VertexKey};
use crate::domain::core::scalar::Real;
use nalgebra::{Point3, UnitVector3, Vector3};

// ── Boundary patch types ──────────────────────────────────────────────────────

/// Physical role of a CFD boundary patch.
///
/// Each named boundary patch groups faces that share the same physical boundary
/// condition for the CFD solver (inlet velocity, outlet pressure, wall, etc.).
///
/// # CFD Semantics
///
/// | Variant | OpenFOAM type | Typical BC |
/// |---------|--------------|------------|
/// | `Inlet` | `patch` (physicalType: inlet) | Fixed velocity |
/// | `Outlet` | `patch` (physicalType: outlet) | Fixed pressure |
/// | `Wall` | `wall` | No-slip or slip |
/// | `Symmetry` | `symmetry` | Zero-gradient normal |
/// | `Periodic` | `cyclicAMI` | Matched pair |
/// | `Channel` | `patch` | Millifluidic channel boundary |
/// | `Custom` | user-defined | Any |
#[derive(Debug, Clone, PartialEq, Eq, Hash, serde::Serialize, serde::Deserialize)]
pub enum PatchType {
    /// Velocity or pressure inlet.
    Inlet,
    /// Pressure or zero-gradient outlet.
    Outlet,
    /// No-slip or slip wall boundary.
    Wall,
    /// Symmetry plane (zero-gradient normal flux).
    Symmetry,
    /// Periodic boundary — must come in matched pairs.
    Periodic,
    /// Millifluidic channel boundary (legacy label from builders).
    Channel,
    /// Custom named boundary type.
    Custom(String),
}

impl PatchType {
    /// Returns the OpenFOAM boundary type string for this patch.
    pub fn openfoam_type(&self) -> &str {
        match self {
            PatchType::Wall => "wall",
            PatchType::Symmetry => "symmetry",
            PatchType::Periodic => "cyclicAMI",
            _ => "patch",
        }
    }
}

/// A named CFD boundary patch.
///
/// Groups a set of boundary faces with a physical role ([`PatchType`]) and a
/// human-readable name (e.g. `"inlet"`, `"top_wall"`).
///
/// # Invariants
/// - `name` is non-empty.
/// - All faces assigned to this patch must be boundary faces (i.e. their
///   bounding half-edges have `face == None` on the exterior side).
#[derive(Debug, Clone, serde::Serialize, serde::Deserialize)]
pub struct BoundaryPatch {
    /// Human-readable name (e.g. `"inlet"`, `"wall_top"`).
    pub name: String,
    /// Physical role of this boundary.
    pub patch_type: PatchType,
}

impl BoundaryPatch {
    /// Create a new boundary patch.
    pub fn new(name: impl Into<String>, patch_type: PatchType) -> Self {
        Self {
            name: name.into(),
            patch_type,
        }
    }
}

// ── Half-edge data ────────────────────────────────────────────────────────────

/// Data stored per directed half-edge.
///
/// The half-edge is the atom of the topology. Every undirected edge `{u,v}` is
/// stored as two half-edges: `he(u→v)` and `he(v→u)`, which are each other's
/// `twin`.
///
/// # Invariants
///
/// - `twin(twin(self)) == self`
/// - `next(prev(self)) == self`
/// - `face(next(self)) == self.face`
///
/// # Diagram
///
/// ```text
///  ... ──► prev(he) ──► he ──► next(he) ──► ...  (face loop, CCW)
///                        │
///                       twin
///                        │
///  ... ◄── prev(twin) ◄── twin(he) ◄── next(twin) ◄── ...  (adjacent face, CW)
/// ```
#[derive(Debug, Clone, Copy)]
pub struct HalfEdgeData {
    /// The **tip** (head) vertex of this directed half-edge.
    ///
    /// For half-edge `u→v`, `vertex` holds the key of `v`.
    pub vertex: VertexKey,

    /// The face to the **left** of this directed half-edge.
    ///
    /// `None` for boundary sentinel half-edges (the "exterior" half of a
    /// boundary edge). Boundary half-edges form closed boundary loops via
    /// `next` links even though they have no face.
    pub face: Option<FaceKey>,

    /// The opposing half-edge for the same undirected edge.
    ///
    /// # Invariant
    /// `twin(twin(he)) == he` always.
    pub twin: HalfEdgeKey,

    /// The next half-edge around the **same face** (counter-clockwise).
    ///
    /// Following `next` repeatedly from any half-edge traces the entire
    /// boundary loop of its face and returns to the starting half-edge.
    pub next: HalfEdgeKey,

    /// The previous half-edge around the **same face** (clockwise).
    ///
    /// `prev` is the inverse of `next`: `next(prev(he)) == he`.
    pub prev: HalfEdgeKey,
}

// ── Vertex data ───────────────────────────────────────────────────────────────

/// Data stored per vertex.
///
/// The `half_edge` field holds *any* outgoing half-edge from this vertex.
/// That is sufficient to traverse the complete one-ring neighborhood of the
/// vertex via `twin → next` links.
///
/// # One-ring traversal
///
/// ```text
/// start = vertex.half_edge          // any outgoing he from v
/// loop:
///     visit start.twin.vertex        // a neighbour of v
///     start = start.twin.next        // next spoke in the fan
/// until start == vertex.half_edge   // ring closed
/// ```
#[derive(Debug, Clone)]
pub struct VertexData {
    /// Position in 3D space.
    pub position: Point3<Real>,

    /// Any outgoing half-edge from this vertex.
    ///
    /// Used as the entry point for one-ring traversal. The specific choice
    /// of which outgoing half-edge does not affect correctness.
    pub half_edge: HalfEdgeKey,
}

impl VertexData {
    /// Create vertex data at `position` with a placeholder half-edge key.
    ///
    /// The `half_edge` field must be updated before any topology traversal.
    pub fn new(position: Point3<Real>) -> Self {
        // Use a dummy key; the mesh builder overwrites this immediately.
        Self {
            position,
            half_edge: HalfEdgeKey::default(),
        }
    }
}

// ── Face data ─────────────────────────────────────────────────────────────────

/// Data stored per face (polygon).
///
/// Typically a triangle in tessellated meshes, but can be an n-gon when
/// importing quadrilateral or polygonal surface descriptions.
///
/// # Traversal
///
/// ```text
/// start = face.half_edge
/// loop:
///     visit start.vertex    // collect vertex of each edge
///     start = start.next
/// until start == face.half_edge
/// ```
#[derive(Debug, Clone)]
pub struct FaceData {
    /// Any half-edge bounding this face.
    ///
    /// The entry point for face-loop traversal.
    pub half_edge: HalfEdgeKey,

    /// The CFD boundary patch this face belongs to, if any.
    ///
    /// `None` for interior faces (faces shared by two cells in a volume mesh,
    /// or non-boundary surface faces). `Some(key)` for boundary faces that
    /// represent inlet/outlet/wall/etc. boundaries.
    pub patch: Option<PatchKey>,

    /// Outward-pointing unit normal (recomputed on geometry changes).
    ///
    /// Cached for performance. Recompute via cross product of face edges
    /// whenever vertex positions change.
    pub normal: UnitVector3<Real>,
}

impl FaceData {
    /// Create face data with a given entry half-edge and precomputed normal.
    pub fn new(half_edge: HalfEdgeKey, normal: UnitVector3<Real>) -> Self {
        Self {
            half_edge,
            patch: None,
            normal,
        }
    }

    /// Create face data with a zero normal (placeholder; call `recompute_normal`).
    pub fn with_half_edge(half_edge: HalfEdgeKey) -> Self {
        Self {
            half_edge,
            patch: None,
            normal: UnitVector3::new_unchecked(Vector3::z()),
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn patch_type_openfoam_strings() {
        assert_eq!(PatchType::Wall.openfoam_type(), "wall");
        assert_eq!(PatchType::Symmetry.openfoam_type(), "symmetry");
        assert_eq!(PatchType::Periodic.openfoam_type(), "cyclicAMI");
        assert_eq!(PatchType::Inlet.openfoam_type(), "patch");
        assert_eq!(PatchType::Outlet.openfoam_type(), "patch");
    }

    #[test]
    fn boundary_patch_name() {
        let p = BoundaryPatch::new("inlet", PatchType::Inlet);
        assert_eq!(p.name, "inlet");
        assert_eq!(p.patch_type, PatchType::Inlet);
    }

    #[test]
    fn vertex_data_position() {
        let v = VertexData::new(Point3::new(1.0, 2.0, 3.0));
        assert_eq!(v.position.x, 1.0);
        assert_eq!(v.position.y, 2.0);
        assert_eq!(v.position.z, 3.0);
    }
}
