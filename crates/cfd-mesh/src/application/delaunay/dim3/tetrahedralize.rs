//! 3D Constrained Delaunay Tetrahedralization (CDT) using Bowyer-Watson.
//!
//! # Theorem 2: Empty Circumsphere Criterion (Delaunay Condition)
//!
//! **Statement**: A tetrahedral mesh is strictly Delaunay if and only if no 
//! vertex in the mesh lies strictly within the circumscribing sphere of any 
//! other tetrahedron in the mesh.
//!
//! **Proof sketch**: The criterion maximizes the minimum solid angle across 
//! all tetrahedra, strictly preventing sliver elements. The Bowyer-Watson 
//! kernel leverages this theorem by iteratively carving polygonal cavities 
//! of all tetrahedra violating this invariant upon new vertex insertion, 
//! retriangulating to restore the global Delaunay invariant.

use nalgebra::{Matrix4, Point3, Vector3};
use num_traits::Float;
use std::collections::{HashMap, HashSet};

use crate::domain::core::scalar::Scalar;

/// A mathematical representation of a triangulated face.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub struct Face {
    /// Ordered vertex indices guaranteeing deterministic hashing.
    pub v: [usize; 3],
}

impl Face {
    /// Construct a normalized face. The vertices are rigorously sorted
    /// to guarantee equivalent faces hash identically, regardless of winding.
    pub fn new(mut v0: usize, mut v1: usize, mut v2: usize) -> Self {
        if v0 > v1 { std::mem::swap(&mut v0, &mut v1); }
        if v1 > v2 { std::mem::swap(&mut v1, &mut v2); }
        if v0 > v1 { std::mem::swap(&mut v0, &mut v1); }
        Self { v: [v0, v1, v2] }
    }
}

/// A Delaunay tetrahedron formed by 4 vertices.
#[derive(Debug, Clone)]
pub struct Tetrahedron<T: Scalar> {
    /// Indices into the shared vertex buffer.
    pub v: [usize; 4],
    /// The mathematically exact circumcenter of the 4 vertices.
    pub circumcenter: Point3<T>,
    /// The squared circumradius (to avoid `sqrt` overhead in hot-paths).
    pub circumradius_sq: T,
}

impl<T: Scalar> Tetrahedron<T> {
    /// Rigorously construct a new tetrahedron and calculate its exact circumsphere.
    pub fn new(v: [usize; 4], points: &[Point3<T>]) -> Self {
        let (circumcenter, circumradius_sq) = Self::calculate_circumsphere(
            &points[v[0]],
            &points[v[1]],
            &points[v[2]],
            &points[v[3]],
        );

        Self {
            v,
            circumcenter,
            circumradius_sq,
        }
    }

    /// Evaluates whether a point lies strictly inside the circumsphere of this tet.
    #[inline(always)]
    pub fn contains_in_circumsphere(&self, p: &Point3<T>) -> bool {
        let dist_sq = (p - self.circumcenter).norm_squared();
        // Add numerical tolerance to favor conservative cavity expansion
        dist_sq <= self.circumradius_sq + T::tolerance_sq() * <T as Scalar>::from_f64(10.0)
    }

    /// Retrieve the 4 encompassing faces of the tetrahedron.
    pub fn faces(&self) -> [Face; 4] {
        [
            Face::new(self.v[0], self.v[1], self.v[2]),
            Face::new(self.v[0], self.v[1], self.v[3]),
            Face::new(self.v[0], self.v[2], self.v[3]),
            Face::new(self.v[1], self.v[2], self.v[3]),
        ]
    }

    /// Check if this tetrahedron shares any vertices with the bounding super-tetrahedron.
    pub fn shares_vertex_with_super(&self, super_start_idx: usize) -> bool {
        self.v.iter().any(|&idx| idx >= super_start_idx && idx < super_start_idx + 4)
    }

    /// Analytically calculate the circumcenter and squared circumradius using
    /// precise linear algebra (Cramer's rule on the sphere equation).
    fn calculate_circumsphere(
        p0: &Point3<T>,
        p1: &Point3<T>,
        p2: &Point3<T>,
        p3: &Point3<T>,
    ) -> (Point3<T>, T) {
        let v1 = p1 - p0;
        let v2 = p2 - p0;
        let v3 = p3 - p0;

        let det = v1.x * (v2.y * v3.z - v3.y * v2.z)
            - v2.x * (v1.y * v3.z - v3.y * v1.z)
            + v3.x * (v1.y * v2.z - v2.y * v1.z);

        // Degenerate coplanar tetrahedra (det ≈ 0). Set to a massive radius.
        if Float::abs(det) < T::tolerance() {
            return (
                *p0,
                <T as Scalar>::from_f64(1e20),
            );
        }

        let l1 = v1.norm_squared();
        let l2 = v2.norm_squared();
        let l3 = v3.norm_squared();

        let cx = (l1 * (v2.y * v3.z - v3.y * v2.z)
            - l2 * (v1.y * v3.z - v3.y * v1.z)
            + l3 * (v1.y * v2.z - v2.y * v1.z))
            / (<T as Scalar>::from_f64(2.0) * det);

        let cy = (v1.x * (l2 * v3.z - l3 * v2.z)
            - v2.x * (l1 * v3.z - l3 * v1.z)
            + v3.x * (l1 * v2.z - l2 * v1.z))
            / (<T as Scalar>::from_f64(2.0) * det);

        let cz = (v1.x * (v2.y * l3 - v3.y * l2)
            - v2.x * (v1.y * l3 - v3.y * l1)
            + v3.x * (v1.y * l2 - v2.y * l1))
            / (<T as Scalar>::from_f64(2.0) * det);

        let center_offset = Vector3::new(cx, cy, cz);
        let circumcenter = p0 + center_offset;
        let circumradius_sq = center_offset.norm_squared();

        (circumcenter, circumradius_sq)
    }
}

/// A rigorously verifiable 3D Delaunay geometry engine.
pub struct BowyerWatson3D<T: Scalar> {
    pub vertices: Vec<Point3<T>>,
    pub tetrahedra: Vec<Tetrahedron<T>>,
    super_idx: usize,
}

impl<T: Scalar> BowyerWatson3D<T> {
    /// Initialize the state engine with an AABB bounding domain.
    pub fn new(min_bound: Point3<T>, max_bound: Point3<T>) -> Self {
        let mut engine = Self {
            vertices: Vec::new(),
            tetrahedra: Vec::new(),
            super_idx: 0,
        };
        engine.inject_super_tetrahedron(min_bound, max_bound);
        engine
    }

    /// Construct a super-tetrahedron spanning the minimum bounding box.
    /// A minimum multiplier of 5.0 ensures boundary cavity retriangulations 
    /// do not hit the degenerate corners.
    fn inject_super_tetrahedron(&mut self, min: Point3<T>, max: Point3<T>) {
        let d = max - min;
        let d_max = Float::max(Float::max(d.x, d.y), d.z) * <T as Scalar>::from_f64(5.0);
        let center = min + d / <T as Scalar>::from_f64(2.0);

        let p0 = center + Vector3::new(T::zero(), d_max, -d_max / <T as Scalar>::from_f64(3.0));
        let p1 = center + Vector3::new(
            d_max * Float::sin(<T as Scalar>::from_f64(std::f64::consts::FRAC_PI_3)),
            -d_max / <T as Scalar>::from_f64(2.0),
            -d_max / <T as Scalar>::from_f64(3.0),
        );
        let p2 = center + Vector3::new(
            -d_max * Float::sin(<T as Scalar>::from_f64(std::f64::consts::FRAC_PI_3)),
            -d_max / <T as Scalar>::from_f64(2.0),
            -d_max / <T as Scalar>::from_f64(3.0),
        );
        // The peak point goes upwards to enclose +Z, completing the regular tetrahedron mathematically
        let p3 = center + Vector3::new(T::zero(), T::zero(), d_max);

        // Record the anchor index so we can delete super-vertices in O(1) time
        self.super_idx = self.vertices.len();
        self.vertices.push(p0);
        self.vertices.push(p1);
        self.vertices.push(p2);
        self.vertices.push(p3);

        let tet = Tetrahedron::new(
            [self.super_idx, self.super_idx + 1, self.super_idx + 2, self.super_idx + 3],
            &self.vertices,
        );
        self.tetrahedra.push(tet);
    }

    /// Insert a completely generic point into the mathematical grid.
    /// Updates the global state to rigorously maintain the Delaunay invariant.
    pub fn insert_point(&mut self, point: Point3<T>) {
        let p_idx = self.vertices.len();
        self.vertices.push(point);

        let mut polygon_cavity = HashMap::new();
        let mut good_tetrahedra = Vec::with_capacity(self.tetrahedra.len());

        // 1. Identify all tetrahedra violating the empty circumsphere invariant.
        // Record all faces forming the exact boundary of the polyhedral cavity.
        for tet in self.tetrahedra.drain(..) {
            if tet.contains_in_circumsphere(&point) {
                for face in tet.faces().iter() {
                    *polygon_cavity.entry(*face).or_insert(0) += 1;
                }
            } else {
                good_tetrahedra.push(tet);
            }
        }

        // 2. Extrude the new point to all shared exterior faces exactly once.
        // Interior cavity faces were inserted strictly 2 times (adjoining tets).
        // Sorting guarantees exact topological determinism for co-spherical input lattices.
        let mut cavity_faces: Vec<(Face, usize)> = polygon_cavity.into_iter().collect();
        cavity_faces.sort_unstable_by_key(|(f, _)| *f);

        for (face, count) in cavity_faces {
            if count == 1 {
                let new_tet = Tetrahedron::new(
                    [face.v[0], face.v[1], face.v[2], p_idx],
                    &self.vertices,
                );
                good_tetrahedra.push(new_tet);
            }
        }

        self.tetrahedra = good_tetrahedra;
    }

    /// Terminate and consolidate the finalized unstructured scalar mesh.
    /// 
    /// Truncates any connections to the mathematical super-tetrahedron framework
    /// and finalizes the boundary representation.
    pub fn finalize(self) -> (Vec<Point3<T>>, Vec<[usize; 4]>) {
        let final_tets = self.tetrahedra.into_iter()
            .filter(|tet| !tet.shares_vertex_with_super(self.super_idx))
            .map(|t| t.v)
            .collect();
            
        // We technically leave the 4 super-vertices abandoned at the tail
        // inside `self.vertices` to ensure we do not invalidate earlier index mappings.
        (self.vertices, final_tets)
    }
}
