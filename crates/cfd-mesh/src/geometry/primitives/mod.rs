//! Analytic mesh primitives.
//!
//! Each primitive builder produces a **watertight, outward-oriented
//! [`IndexedMesh`]** that can be immediately used with the CSG pipeline,
//! STL export, or watertight validation.
//!
//! ## Included primitives
//!
//! | Builder | Shape | Analytical volume |
//! |---------|-------|-------------------|
//! | [`Tetrahedron`] | Regular tetrahedron | `(8√3)/3 · r³` |
//! | [`Cube`] | Axis-aligned box | `a · b · c` |
//! | [`UvSphere`] | UV-parametric sphere | `4π r³ / 3` |
//! | [`Cylinder`] | Closed right cylinder | `π r² h` |
//! | [`Cone`] | Right circular cone | `π r² h / 3` |
//! | [`Torus`] | Ring torus | `2π² R r²` |
//! | [`LinearSweep`] | Prismatic solid (polygon × extrusion) | `A_profile × h` |
//! | [`RevolutionSweep`] | Solid of revolution (Pappus) | `2π R̄ A angle/(2π)` |
//!
//! ## Winding convention
//!
//! All primitives use **outward CCW** face winding: when viewed from outside
//! the solid, vertices go counter-clockwise. The right-hand rule gives the
//! outward surface normal. `signed_volume > 0` for all primitives.
//!
//! ## Example
//!
//! ```rust,ignore
//! use cfd_mesh::geometry::primitives::UvSphere;
//!
//! let mesh = UvSphere { radius: 1.0, segments: 32, stacks: 16 }
//!     .build()
//!     .expect("sphere");
//! assert!(mesh.signed_volume() > 0.0);
//! ```

pub mod tetrahedron;
pub mod cube;
pub mod sphere;
pub mod cylinder;
pub mod cone;
pub mod torus;
pub mod linear_sweep;
pub mod revolution_sweep;

pub use tetrahedron::Tetrahedron;
pub use cube::Cube;
pub use sphere::UvSphere;
pub use cylinder::Cylinder;
pub use cone::Cone;
pub use torus::Torus;
pub use linear_sweep::LinearSweep;
pub use revolution_sweep::RevolutionSweep;

use crate::mesh::IndexedMesh;

/// Common error type for primitive mesh builders.
#[derive(Debug, thiserror::Error)]
pub enum PrimitiveError {
    /// A dimension or resolution parameter was zero or negative.
    #[error("invalid parameter: {0}")]
    InvalidParam(String),
    /// Fewer than 3 angular segments were requested.
    #[error("segments must be ≥ 3, got {0}")]
    TooFewSegments(usize),
    /// The resulting mesh is not watertight (internal consistency failure).
    #[error("internal mesh error: {0}")]
    Mesh(#[from] crate::core::error::MeshError),
}

/// Trait for all primitive mesh builders.
pub trait PrimitiveMesh {
    /// Build the mesh, returning a watertight, outward-oriented [`IndexedMesh`].
    fn build(&self) -> Result<IndexedMesh, PrimitiveError>;
}
