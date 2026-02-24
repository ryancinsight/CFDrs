//! Normal-orientation analysis for `IndexedMesh` surfaces.
//!
//! Provides [`NormalAnalysis`] and [`analyze_normals`] — routinely used by CSG
//! examples and validation tools to report face-winding consistency and
//! vertex-normal alignment across a mesh.
//!
//! ## Algorithm
//!
//! 1. Compute the mesh **centroid** (mean vertex position).
//! 2. For each non-degenerate triangle:
//!    - Compute the **face normal** via the cross product of two edges.
//!    - Dot the face normal with the centroid-to-face-centre vector.
//!      Positive → outward-facing; negative → inward-facing.
//!    - Compute the **average stored vertex normal** for the face and dot it
//!      with the face normal to obtain the face–vertex alignment score.
//! 3. Aggregate counts and alignment statistics.
//!
//! ## Interpretation
//!
//! | `inward_faces / total_faces` | Likely cause                              |
//! |------------------------------|-------------------------------------------|
//! | 0%                           | All faces outward — ideal convex mesh     |
//! | < 5%                         | Acceptable; seam artefacts common in CSG  |
//! | > 10%                        | Winding problem; check Boolean op result  |
//! | ≈ 50%                        | Mixed winding; mesh likely non-manifold   |
//!
//! `face_vertex_alignment_mean` near 1.0 means stored vertex normals agree with
//! computed face normals (good for smooth-shaded rendering and CFD post-processing).

use crate::domain::core::scalar::{Point3r, Real, Vector3r};
use crate::domain::geometry::normal::triangle_normal;
use crate::domain::mesh::IndexedMesh;

// ── Public types ──────────────────────────────────────────────────────────────

/// Per-mesh normal-orientation statistics.
///
/// Returned by [`analyze_normals`].
///
/// # Invariants
///
/// - `outward_faces + inward_faces + degenerate_faces == total triangles checked`
/// - `face_vertex_alignment_mean` ∈ [−1, 1]; 1.0 = perfect agreement
/// - `face_vertex_alignment_min`  ∈ [−1, 1]; < 0 indicates at least one
///   face whose stored vertex normals point opposite to the winding normal
#[derive(Debug, Clone, PartialEq)]
pub struct NormalAnalysis {
    /// Number of faces whose computed normal points away from the mesh centroid.
    pub outward_faces: usize,
    /// Number of faces whose computed normal points toward the mesh centroid.
    pub inward_faces: usize,
    /// Number of degenerate (zero-area) faces skipped during analysis.
    pub degenerate_faces: usize,
    /// Mean dot product of the computed face normal vs the averaged stored
    /// vertex normals across all non-degenerate faces.
    pub face_vertex_alignment_mean: Real,
    /// Minimum dot product across all non-degenerate faces.
    pub face_vertex_alignment_min: Real,
}

impl NormalAnalysis {
    /// Total number of faces inspected (degenerate faces included).
    #[inline]
    pub fn total_faces(&self) -> usize {
        self.outward_faces + self.inward_faces + self.degenerate_faces
    }

    /// Fraction of non-degenerate faces that are inward-facing (0.0 – 1.0).
    ///
    /// Returns `0.0` when the mesh is empty.
    #[inline]
    pub fn inward_fraction(&self) -> Real {
        let n = (self.outward_faces + self.inward_faces) as Real;
        if n > 0.0 {
            self.inward_faces as Real / n
        } else {
            0.0
        }
    }

    /// Returns `true` when every non-degenerate face is outward-facing.
    #[inline]
    pub fn all_outward(&self) -> bool {
        self.inward_faces == 0
    }
}

// ── Public function ───────────────────────────────────────────────────────────

/// Analyse the normal orientation of every face in `mesh`.
///
/// Uses the mesh centroid as a reference point to classify faces as outward or
/// inward.  This heuristic is reliable for convex and near-convex meshes; for
/// highly non-convex shapes (e.g. a torus interior or a CSG `Difference` solid
/// with re-entrant cavities), a small number of geometrically-correct faces
/// near concavities may be falsely classified as inward.
///
/// # Arguments
///
/// * `mesh` — The surface mesh to analyse.  Modified only if the caller is
///   passing `&mut IndexedMesh`; this function takes a shared reference.
///
/// # Returns
///
/// A [`NormalAnalysis`] struct with per-category counts and alignment
/// statistics.
///
/// # Examples
///
/// ```rust,ignore
/// use cfd_mesh::{UvSphere, geometry::primitives::PrimitiveMesh};
/// use cfd_mesh::application::quality::normals::analyze_normals;
///
/// let sphere = UvSphere { radius: 1.0, segments: 32, stacks: 16, ..Default::default() }
///     .build().unwrap();
/// let report = analyze_normals(&sphere);
/// assert_eq!(report.inward_faces, 0, "sphere should be all-outward");
/// ```
pub fn analyze_normals(mesh: &IndexedMesh) -> NormalAnalysis {
    // ── Step 1: compute mesh centroid ────────────────────────────────────────
    let mut centroid_sum = Vector3r::zeros();
    let mut cnt = 0usize;
    for (_, v) in mesh.vertices.iter() {
        centroid_sum += v.position.coords;
        cnt += 1;
    }
    let center = if cnt > 0 {
        Point3r::from(centroid_sum / cnt as Real)
    } else {
        Point3r::origin()
    };

    // ── Step 2: per-face statistics ──────────────────────────────────────────
    let mut outward = 0usize;
    let mut inward = 0usize;
    let mut degen = 0usize;
    let mut asum: Real = 0.0;
    let mut acnt = 0usize;
    let mut amin: Real = 1.0;

    for face in mesh.faces.iter() {
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);

        let Some(face_n) = triangle_normal(a, b, c) else {
            degen += 1;
            continue;
        };

        // Classify by direction from centroid to face centre.
        let fc = Point3r::new(
            (a.x + b.x + c.x) / 3.0,
            (a.y + b.y + c.y) / 3.0,
            (a.z + b.z + c.z) / 3.0,
        );
        let dir = fc - center;
        if dir.norm() > 1e-12 {
            if face_n.dot(&dir.normalize()) >= 0.0 {
                outward += 1;
            } else {
                inward += 1;
            }
        }

        // Face ↔ vertex-normal alignment.
        let avg_n = (*mesh.vertices.normal(face.vertices[0])
            + *mesh.vertices.normal(face.vertices[1])
            + *mesh.vertices.normal(face.vertices[2]))
            / 3.0;
        let l = avg_n.norm();
        if l > 1e-12 {
            let al = face_n.dot(&(avg_n / l));
            asum += al;
            acnt += 1;
            amin = amin.min(al);
        }
    }

    NormalAnalysis {
        outward_faces: outward,
        inward_faces: inward,
        degenerate_faces: degen,
        face_vertex_alignment_mean: if acnt > 0 { asum / acnt as Real } else { 0.0 },
        face_vertex_alignment_min: if acnt > 0 { amin } else { 0.0 },
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::primitives::{Cube, PrimitiveMesh, UvSphere};

    #[test]
    fn sphere_all_outward() {
        let mesh = UvSphere {
            radius: 1.0,
            segments: 32,
            stacks: 16,
            ..Default::default()
        }
        .build()
        .unwrap();
        let r = analyze_normals(&mesh);
        assert_eq!(r.inward_faces, 0, "UV sphere should have zero inward faces");
        assert_eq!(
            r.degenerate_faces, 0,
            "UV sphere should have no degenerate faces"
        );
        assert!(
            r.face_vertex_alignment_mean > 0.9,
            "face-vertex alignment mean should be > 0.9, got {}",
            r.face_vertex_alignment_mean
        );
    }

    #[test]
    fn cube_all_outward() {
        let mesh = Cube {
            origin: Point3r::origin(),
            width: 2.0,
            height: 2.0,
            depth: 2.0,
        }
        .build()
        .unwrap();
        let r = analyze_normals(&mesh);
        assert_eq!(r.inward_faces, 0, "cube should have zero inward faces");
    }

    #[test]
    fn empty_mesh_returns_zeros() {
        let mesh = IndexedMesh::new();
        let r = analyze_normals(&mesh);
        assert_eq!(r.outward_faces, 0);
        assert_eq!(r.inward_faces, 0);
        assert_eq!(r.degenerate_faces, 0);
        assert_eq!(r.face_vertex_alignment_mean, 0.0);
        assert_eq!(r.face_vertex_alignment_min, 0.0);
    }

    #[test]
    fn inward_fraction_zero_on_clean_mesh() {
        let mesh = UvSphere {
            radius: 1.0,
            segments: 16,
            stacks: 8,
            ..Default::default()
        }
        .build()
        .unwrap();
        let r = analyze_normals(&mesh);
        assert_eq!(r.inward_fraction(), 0.0);
        assert!(r.all_outward());
    }

    #[test]
    fn total_faces_matches_mesh() {
        let mesh = UvSphere {
            radius: 1.0,
            segments: 16,
            stacks: 8,
            ..Default::default()
        }
        .build()
        .unwrap();
        let r = analyze_normals(&mesh);
        assert_eq!(r.total_faces(), mesh.face_count());
    }
}
