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
//! | [`Octahedron`] | Regular octahedron | `(4/3) R³` |
//! | [`Icosahedron`] | Regular icosahedron | `(5(3+√5)/12) a³` |
//! | [`Ellipsoid`] | Triaxial ellipsoid | `(4/3)π a b c` |
//! | [`Frustum`] | Truncated cone | `(π h/3)(r₀²+r₀r₁+r₁²)` |
//! | [`Capsule`] | Cylinder + hemisphere caps | `π r²(h + 4r/3)` |
//! | [`Pipe`] | Hollow cylinder (χ=0) | `π(r_o²−r_i²)h` |
//! | [`Elbow`] | Circular-arc pipe bend | `π r² R θ` |
//! | [`BiconcaveDisk`] | Evans-Fung RBC shape | ~94 fL at 8 µm |
//! | [`SphericalShell`] | Hollow sphere (χ=0) | `(4π/3)(r_o³−r_i³)` |
//! | [`StadiumPrism`] | Rounded-rectangle prism | `(π r² + 2r·flat)·h` |
//! | [`Dodecahedron`] | Regular dodecahedron | `(1/4)(15+7√5) a³` |
//! | [`GeodesicSphere`] | Subdivided icosahedron | `≈ (4π/3) R³` |
//! | [`HelixSweep`] | Helical tube sweep | `π r² · arc_length` |
//! | [`RoundedCube`] | Filleted box | `≈ w·h·d` |
//! | [`Cuboctahedron`] | Archimedean solid | `(5√2/3) a³` |
//! | [`Pyramid`] | Right n-gon pyramid | `A_base h / 3` |
//! | [`Antiprism`] | n-gon antiprism | varies |
//! | [`TruncatedIcosahedron`] | Soccer ball / C₆₀ | varies |
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
//! use cfd_mesh::domain::geometry::primitives::UvSphere;
//!
//! let mesh = UvSphere { radius: 1.0, segments: 32, stacks: 16 }
//!     .build()
//!     .expect("sphere");
//! assert!(mesh.signed_volume() > 0.0);
//! ```

pub mod antiprism;
pub mod biconcave_disk;
pub mod capsule;
pub mod cone;
pub mod cube;
pub mod cuboctahedron;
pub mod cylinder;
pub mod disk;
pub mod dodecahedron;
pub mod elbow;
pub mod ellipsoid;
pub mod frustum;
pub mod geodesic_sphere;
pub mod helix_sweep;
pub mod icosahedron;
pub mod linear_sweep;
pub mod octahedron;
pub mod pipe;
pub mod pyramid;
pub mod revolution_sweep;
pub mod rounded_cube;
pub mod sphere;
pub mod spherical_shell;
pub mod stadium_prism;
pub mod tetrahedron;
pub mod torus;
pub mod truncated_icosahedron;

pub use antiprism::Antiprism;
pub use biconcave_disk::BiconcaveDisk;
pub use capsule::Capsule;
pub use cone::Cone;
pub use cube::Cube;
pub use cuboctahedron::Cuboctahedron;
pub use cylinder::Cylinder;
pub use disk::Disk;
pub use dodecahedron::Dodecahedron;
pub use elbow::Elbow;
pub use ellipsoid::Ellipsoid;
pub use frustum::Frustum;
pub use geodesic_sphere::GeodesicSphere;
pub use helix_sweep::HelixSweep;
pub use icosahedron::Icosahedron;
pub use linear_sweep::LinearSweep;
pub use octahedron::Octahedron;
pub use pipe::Pipe;
pub use pyramid::Pyramid;
pub use revolution_sweep::RevolutionSweep;
pub use rounded_cube::RoundedCube;
pub use sphere::UvSphere;
pub use spherical_shell::SphericalShell;
pub use stadium_prism::StadiumPrism;
pub use tetrahedron::Tetrahedron;
pub use torus::Torus;
pub use truncated_icosahedron::TruncatedIcosahedron;

use crate::domain::mesh::IndexedMesh;

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
    Mesh(#[from] crate::domain::core::error::MeshError),
}

/// Trait for all primitive mesh builders.
pub trait PrimitiveMesh {
    /// Build the mesh, returning a watertight, outward-oriented [`IndexedMesh`].
    fn build(&self) -> Result<IndexedMesh, PrimitiveError>;
}
