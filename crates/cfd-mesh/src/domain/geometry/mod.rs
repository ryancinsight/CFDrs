//! Geometric primitives, exact predicates, and mesh builders.

pub mod aabb;
pub mod measure;
pub mod normal;
pub mod plane;
pub mod predicates;

pub use aabb::Aabb;
pub use plane::Plane;
pub use predicates::Orientation;


pub mod nurbs;
pub use nurbs::{NurbsCurve, NurbsSurface, TessellationOptions};

pub mod primitives;
pub use primitives::{
    Antiprism, BiconcaveDisk, Capsule, Cone, Cube, Cuboctahedron, Cylinder, Disk, Dodecahedron,
    Elbow, Ellipsoid, Frustum, GeodesicSphere, HelixSweep, Icosahedron, LinearSweep, Octahedron,
    Pipe, Pyramid, RevolutionSweep, RoundedCube, SerpentineTube, SphericalShell, StadiumPrism,
    Tetrahedron, Torus, TruncatedIcosahedron, UvSphere,
};
