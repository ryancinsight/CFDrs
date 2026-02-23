//! Geometric primitives, exact predicates, and mesh builders.

pub mod plane;
pub mod aabb;
pub mod normal;
pub mod measure;
pub mod distance;
pub mod predicates;

pub use plane::Plane;
pub use aabb::Aabb;
pub use predicates::Orientation;

pub mod venturi;
pub mod serpentine;
pub mod branching;

pub mod nurbs;
pub use nurbs::{NurbsCurve, NurbsSurface, TessellationOptions};

pub mod primitives;
pub use primitives::{
    Tetrahedron, Cube, UvSphere, Cylinder, Cone, Torus, LinearSweep, RevolutionSweep,
    Octahedron, Icosahedron, Ellipsoid, Frustum, Capsule, Pipe, Elbow,
    BiconcaveDisk, SphericalShell, StadiumPrism, Dodecahedron, GeodesicSphere,
    HelixSweep, RoundedCube, Cuboctahedron, Pyramid, Antiprism, TruncatedIcosahedron,
    Disk,
};
