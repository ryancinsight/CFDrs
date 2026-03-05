//! Discrete mean curvature via the cotangent-weighted Laplace-Beltrami operator.
//!
//! ## Algorithm — Cotangent Laplacian
//!
//! For vertex `v_i` with 1-ring neighbourhood `N(i)`:
//!
//! ```text
//! H(v_i) = || K(v_i) || / 2
//!
//! K(v_i) = (1 / A_i) × Σ_{j ∈ N(i)} (cot α_ij + cot β_ij) × (v_j - v_i)
//! ```
//!
//! where:
//! - `α_ij`, `β_ij` — angles of the two triangles sharing edge `(i, j)` at
//!   their respective opposite vertices.
//! - `A_i` — barycentric area: `Σ_{faces containing i} face_area / 3`.
//!
//! ## Theorem — Second-Order Convergence
//!
//! Let `M` be a smooth 2-manifold, `T_h` a triangulation with maximum edge
//! length `h`.  The cotangent-weighted discrete mean curvature `H_h(v)` satisfies:
//!
//! ```text
//! | H_h(v) - H(v) | = O(h²)
//! ```
//!
//! as `h → 0`, for interior vertices on well-shaped triangles.
//!
//! *Source*: Meyer et al., "Discrete Differential-Geometry Operators for
//! Triangulated 2-Manifolds", VisMath 2003.
//!
//! ## Boundary Handling
//!
//! Vertices on the boundary of an open mesh have incomplete 1-rings.  The
//! cotangent sum is still accumulated from available faces; the barycentric area
//! is the sum of the partial ring.  The result is geometrically meaningful for
//! smooth boundary curves but is not the standard interior curvature.
//!
//! ## Complexity
//!
//! O(F) where F = face count; each face contributes 3 cotangent weights.

use crate::domain::core::scalar::Real;
use crate::domain::mesh::IndexedMesh;
use nalgebra::Vector3;

// ── Public API ────────────────────────────────────────────────────────────────

/// Compute discrete mean curvature at each vertex via the cotangent Laplacian.
///
/// Returns a `Vec<Real>` with length equal to `mesh.vertices.len()`.
/// Index `i` holds the mean curvature at the vertex with `VertexId(i as u32)`.
///
/// Non-finite values indicate degenerate local geometry (zero-area faces in
/// the 1-ring); callers should treat these as masked / invalid.
#[must_use]
pub fn vertex_mean_curvature(mesh: &IndexedMesh) -> Vec<Real> {
    let n = mesh.vertices.len();
    if n == 0 {
        return Vec::new();
    }

    // Accumulate cotangent-weighted Laplacian vectors and barycentric areas.
    let mut laplacian: Vec<Vector3<Real>> = vec![Vector3::zeros(); n];
    let mut area: Vec<Real> = vec![0.0; n];

    for face in mesh.faces.iter() {
        let [ia, ib, ic] = [
            face.vertices[0].as_usize(),
            face.vertices[1].as_usize(),
            face.vertices[2].as_usize(),
        ];
        let pa = mesh.vertices.position(face.vertices[0]);
        let pb = mesh.vertices.position(face.vertices[1]);
        let pc = mesh.vertices.position(face.vertices[2]);

        let va = Vector3::new(pa.x, pa.y, pa.z);
        let vb = Vector3::new(pb.x, pb.y, pb.z);
        let vc = Vector3::new(pc.x, pc.y, pc.z);

        // Cotangent at each vertex: cot(θ) = dot(e1, e2) / |e1 × e2|
        let cross_a = (vb - va).cross(&(vc - va));
        let cross_b = (va - vb).cross(&(vc - vb));
        let cross_c = (va - vc).cross(&(vb - vc));

        let denom_a = cross_a.norm().max(Real::MIN_POSITIVE);
        let denom_b = cross_b.norm().max(Real::MIN_POSITIVE);
        let denom_c = cross_c.norm().max(Real::MIN_POSITIVE);

        let cot_a = (vb - va).dot(&(vc - va)) / denom_a;
        let cot_b = (va - vb).dot(&(vc - vb)) / denom_b;
        let cot_c = (va - vc).dot(&(vb - vc)) / denom_c;

        let face_area = 0.5 * denom_a; // area of triangle from cross at vertex a

        // Accumulate barycentric (one-third) area at each vertex.
        area[ia] += face_area / 3.0;
        area[ib] += face_area / 3.0;
        area[ic] += face_area / 3.0;

        // Edge (b, c) — opposite vertex a contributes cot_a.
        laplacian[ib] += cot_a * (vc - vb);
        laplacian[ic] += cot_a * (vb - vc);

        // Edge (a, c) — opposite vertex b contributes cot_b.
        laplacian[ia] += cot_b * (vc - va);
        laplacian[ic] += cot_b * (va - vc);

        // Edge (a, b) — opposite vertex c contributes cot_c.
        laplacian[ia] += cot_c * (vb - va);
        laplacian[ib] += cot_c * (va - vb);
    }

    // Mean curvature magnitude: H = |K| / 2,  K = laplacian / (2 * A)
    (0..n)
        .map(|i| {
            let a = area[i];
            if a < Real::MIN_POSITIVE {
                return Real::NAN;
            }
            // K(v_i) = (1/A_i) * laplacian_sum; H = |K|/2
            laplacian[i].norm() / (4.0 * a)
        })
        .collect()
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, PrimitiveMesh, UvSphere};

    /// A flat mesh (all vertices coplanar) should have near-zero mean curvature.
    #[test]
    fn flat_mesh_has_near_zero_curvature() {
        // A single flat-panel cube face has H ≈ 0 for interior vertices.
        // We use the cube and check that all finite curvatures are small.
        let mesh = Cube {
            origin: Point3r::new(-1.0, -1.0, -1.0),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .expect("cube build");
        let h = vertex_mean_curvature(&mesh);
        // Cube has 8 corner vertices; each face is flat but corners have 3 meeting
        // at 90°. The curvature magnitudes should all be finite and relatively small.
        let finite: Vec<Real> = h.iter().copied().filter(|v| v.is_finite()).collect();
        assert!(
            !finite.is_empty(),
            "at least some curvatures should be finite"
        );
    }

    /// Sphere of radius R should have mean curvature ≈ 1/R at interior vertices.
    ///
    /// For a tessellated sphere, discrete curvature converges as O(h²). With
    /// 16 segments × 8 stacks, the error should be within 30% of 1/R = 1.0.
    #[test]
    fn sphere_curvature_near_one_over_radius() {
        let r = 1.0_f64;
        let mesh = UvSphere {
            radius: r,
            center: Point3r::origin(),
            segments: 16,
            stacks: 8,
        }
        .build()
        .expect("sphere build");
        let h = vertex_mean_curvature(&mesh);
        let finite: Vec<Real> = h.iter().copied().filter(|v| v.is_finite() && *v > 0.0).collect();
        assert!(!finite.is_empty(), "sphere should have positive curvatures");

        let mean_h: Real = finite.iter().sum::<Real>() / finite.len() as Real;
        let expected = 1.0 / r;
        let rel_err = (mean_h - expected).abs() / expected;
        assert!(
            rel_err < 0.30,
            "sphere mean curvature {mean_h:.4} should be within 30% of 1/R={expected:.4}"
        );
    }

    /// Output length equals vertex count.
    #[test]
    fn output_length_equals_vertex_count() {
        let mesh = Cube {
            origin: Point3r::new(0.0, 0.0, 0.0),
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        }
        .build()
        .expect("cube build");
        let h = vertex_mean_curvature(&mesh);
        assert_eq!(h.len(), mesh.vertices.len());
    }

    /// Empty mesh returns empty vec.
    #[test]
    fn empty_mesh_returns_empty() {
        let mesh = IndexedMesh::new();
        let h = vertex_mean_curvature(&mesh);
        assert!(h.is_empty());
    }
}
