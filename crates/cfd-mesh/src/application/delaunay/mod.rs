//! # Delaunay Triangulation — PSLG / CDT / Ruppert Refinement
//!
//! A complete, self-contained 2-D Constrained Delaunay Triangulation (CDT)
//! engine with Ruppert's mesh refinement, built on the `cfd-mesh` GhostCell /
//! `IndexedMesh` architectural paradigm.  Replaces the external `spade` crate
//! for all triangulation needs.
//!
//! ## Module Hierarchy
//!
//! ```text
//! delaunay/
//! ├── pslg/              Planar Straight-Line Graph (input representation)
//! │   ├── vertex          Vertex type with exact-arithmetic position
//! │   ├── segment         Constraint segment between two vertices
//! │   └── graph           PSLG builder — vertices + segments + holes
//! ├── triangulation/     Core Delaunay kernel
//! │   ├── triangle        Triangle data structure with adjacency slots
//! │   ├── adjacency       Triangle-triangle adjacency bookkeeping
//! │   ├── locate          Point-location via Lawson walk
//! │   └── bowyer_watson   Incremental Bowyer-Watson insertion
//! ├── constraint/        Constrained Delaunay extension
//! │   ├── enforce         Edge-flip constraint enforcement
//! │   └── cavity          Cavity polygon re-triangulation
//! ├── refinement/        Ruppert's algorithm with quality guarantees
//! │   ├── quality         Triangle quality metrics (radius-edge ratio)
//! │   ├── circumcenter    Circumcenter and off-center computation
//! │   ├── encroachment    Diametral lens segment encroachment tests
//! │   └── ruppert         Main refinement loop with termination proof
//! └── convert/           Mesh format conversions
//!     ├── indexed_mesh    Emit an IndexedMesh from the triangulation
//!     └── halfedge_mesh   Emit a HalfEdgeMesh via with_mesh + GhostToken
//! ```
//!
//! ## Mathematical Foundations
//!
//! ### Theorem — Delaunay Empty-Circumcircle Property
//!
//! **Statement**: A triangulation is Delaunay if and only if the open
//! circumdisk of every triangle contains no vertex of the triangulation.
//!
//! **Proof sketch** (Lawson): An edge flip that replaces a non-Delaunay edge
//! `(a, c)` shared by triangles `(a, b, c)` and `(a, c, d)` with `(b, d)`
//! strictly increases the minimum angle of the quad `(a, b, c, d)`.  Since
//! the set of possible triangulations is finite and the minimum angle is
//! bounded, the flip sequence terminates at a Delaunay triangulation.
//! Uniqueness follows from the empty-circumcircle characterisation when no
//! four points are co-circular.
//!
//! ### Theorem — CDT Existence
//!
//! **Statement**: Given a PSLG with no two constraint segments that cross in
//! their interiors, a constrained Delaunay triangulation (CDT) exists and is
//! unique (up to co-circular degeneracies).
//!
//! **Proof sketch**: After inserting all PSLG vertices into a Delaunay
//! triangulation, each missing constraint edge can be recovered by a
//! sequence of edge flips across the cavity of triangles that intersect the
//! constraint.  The resulting triangulation is "Delaunay among all
//! triangulations that contain the constraint edges" — i.e., no vertex
//! lies inside the circumcircle of a triangle unless blocked by a
//! constraint segment.
//!
//! ### Theorem — Ruppert's Refinement Termination & Quality
//!
//! **Statement** (Ruppert 1995, Shewchuk 1997): For any PSLG with minimum
//! input angle ≥ 60° and a user-specified radius-edge ratio bound
//! $B \geq \sqrt{2}$, Ruppert's algorithm terminates in $O(n)$ Steiner
//! point insertions (where $n$ is the output size) and produces a mesh
//! where:
//! 1. Every triangle has radius-edge ratio $\leq B$.
//! 2. No segment is encroached (no vertex lies in the diametral circle).
//! 3. The triangulation is a CDT of the refined PSLG.
//!
//! **Proof sketch**: Each circumcenter insertion either (a) removes a skinny
//! triangle or (b) triggers a segment split that shortens the local feature
//! size.  The local feature size function `lfs(x)` is 1-Lipschitz, so the
//! spacing between Steiner points is bounded below by $lfs(x) / C$ for a
//! constant $C$ depending on $B$.  A packing argument then bounds the total
//! number of insertions.
//!
//! ## Usage
//!
//! ```rust,ignore
//! use cfd_mesh::application::delaunay::{Pslg, Cdt, RuppertRefiner};
//!
//! // Build the PSLG
//! let mut pslg = Pslg::new();
//! let a = pslg.add_vertex(0.0, 0.0);
//! let b = pslg.add_vertex(1.0, 0.0);
//! let c = pslg.add_vertex(0.5, 1.0);
//! pslg.add_segment(a, b);
//! pslg.add_segment(b, c);
//! pslg.add_segment(c, a);
//!
//! // Compute CDT
//! let cdt = Cdt::from_pslg(&pslg);
//!
//! // Refine with Ruppert's algorithm
//! let refined = RuppertRefiner::new(cdt)
//!     .max_radius_edge_ratio(1.5)
//!     .max_area(0.01)
//!     .refine();
//!
//! // Convert to IndexedMesh
//! let mesh = refined.to_indexed_mesh();
//! ```
//!
//! ## Exact Arithmetic
//!
//! All geometric predicates (`orient_2d`, `incircle`) delegate to Shewchuk
//! adaptive-precision arithmetic via [`crate::domain::geometry::predicates`].  This
//! guarantees correct topology even for nearly-degenerate inputs — no
//! epsilon-based fallbacks.

pub mod constraint;
pub mod convert;
pub mod pslg;
pub mod refinement;
pub mod triangulation;

#[cfg(test)]
mod tests;

// ── Primary re-exports ────────────────────────────────────────────────────────

/// The Planar Straight-Line Graph input representation.
pub use pslg::graph::Pslg;
/// PSLG validation error.
pub use pslg::graph::PslgValidationError;

/// A single vertex handle in the PSLG / triangulation.
pub use pslg::vertex::PslgVertexId;

/// A constraint segment handle.
pub use pslg::segment::PslgSegmentId;

/// The core Delaunay triangulation data structure.
pub use triangulation::triangle::TriangleId;

/// Constrained Delaunay Triangulation.
pub use constraint::enforce::Cdt;

/// Ruppert's mesh refinement engine.
pub use refinement::ruppert::RuppertRefiner;

/// Triangle quality metric.
pub use refinement::quality::TriangleQuality;
