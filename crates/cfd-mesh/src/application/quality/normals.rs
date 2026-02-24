//! Normal-orientation analysis for `IndexedMesh` surfaces.
//!
//! Provides [`NormalAnalysis`] and [`analyze_normals`] — routinely used by CSG
//! examples and validation tools to report face-winding consistency and
//! vertex-normal alignment across a mesh.
//!
//! ## Algorithm
//!
//! 1. Build a half-edge adjacency map `(v_i, v_j) → face_idx`.
//! 2. Identify the face with the most extreme vertex (maximum X coordinate).
//!    For this face the outward direction is unambiguous: the face normal must
//!    have a **positive X component** to point away from the solid.
//! 3. Assign that seed face an `Outward` orientation based on sign of `n.x`.
//! 4. **Manifold BFS flood**: propagate orientation to all adjacent faces via
//!    the shared half-edge graph.  Two faces sharing edge (A→B) and (B→A) have
//!    consistent winding (both outward or both inward); sharing (A→B) and (A→B)
//!    (same direction) indicates a winding flip between the two faces.
//! 5. Faces unreachable from the seed (disconnected patches) are re-seeded
//!    from an unvisited extremal face.
//! 6. Count outward / inward from BFS labels; compute vertex-normal alignment.
//!
//! ## Properties
//!
//! - **Correct for non-convex meshes**: CSG difference, tori, concave shapes
//!   all produce `inward_faces = 0` when winding is globally consistent.
//! - **O(F + E)** time and space.
//! - **No centroid**: eliminates the false-positive bias of the old heuristic.
//!
//! ## Interpretation
//!
//! | `inward_faces / total_faces` | Likely cause                              |
//! |------------------------------|-------------------------------------------|
//! | 0%                           | Globally consistent outward winding       |
//! | < 5%                         | Acceptable; isolated CDT seam artefacts   |
//! | > 10%                        | Winding problem; check Boolean op result  |
//! | ≈ 50%                        | Mixed winding; mesh likely non-manifold   |
//!
//! `face_vertex_alignment_mean` near 1.0 means stored vertex normals agree with
//! computed face normals (good for smooth-shaded rendering and CFD post-processing).

use std::collections::{HashMap, VecDeque};

use crate::domain::core::index::VertexId;
use crate::domain::core::scalar::{Real, Vector3r};
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
    /// Number of faces whose computed normal is consistent with the outward
    /// manifold orientation (determined by BFS flood from the extremal seed face).
    pub outward_faces: usize,
    /// Number of faces whose computed normal is inconsistent with the outward
    /// manifold orientation (flipped winding relative to neighbours).
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
/// Uses a **manifold BFS flood** seeded from the face with the most extreme
/// vertex to determine globally consistent outward orientation.  This is
/// correct for any closed orientable 2-manifold, including non-convex CSG
/// difference solids, tori, and re-entrant geometries.
///
/// # Arguments
///
/// * `mesh` — The surface mesh to analyse.  Takes a shared reference.
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
    // ── Step 1: collect face normals and detect degenerates ──────────────────
    let face_list: Vec<_> = mesh.faces.iter().collect();
    let n_faces = face_list.len();

    if n_faces == 0 {
        return NormalAnalysis {
            outward_faces: 0,
            inward_faces: 0,
            degenerate_faces: 0,
            face_vertex_alignment_mean: 0.0,
            face_vertex_alignment_min: 0.0,
        };
    }

    // Per-face computed normals (None = degenerate).
    let mut face_normals: Vec<Option<Vector3r>> = Vec::with_capacity(n_faces);
    for face in &face_list {
        let a = mesh.vertices.position(face.vertices[0]);
        let b = mesh.vertices.position(face.vertices[1]);
        let c = mesh.vertices.position(face.vertices[2]);
        face_normals.push(triangle_normal(a, b, c));
    }

    // ── Step 2: build half-edge adjacency (directed edge → face index) ───────
    //
    // half_edge[(v_i, v_j)] = face_idx of the face that has directed edge i→j.
    // For a manifold mesh every directed edge appears in exactly one face.
    let mut half_edge: HashMap<(VertexId, VertexId), usize> =
        HashMap::with_capacity(n_faces * 3);
    for (fi, face) in face_list.iter().enumerate() {
        let v = face.vertices;
        for k in 0..3 {
            let j = (k + 1) % 3;
            half_edge.insert((v[k], v[j]), fi);
        }
    }

    // ── Step 3: BFS / flood orientation from extremal seed ───────────────────
    //
    // Invariant: `orientation[fi]` = true  → face fi is outward-consistent
    //                               = false → face fi is inward-consistent
    // None = unvisited.
    let mut orientation: Vec<Option<bool>> = vec![None; n_faces];

    // Seed-selection helper: find the non-degenerate face with the vertex
    // carrying the highest X coordinate.  Any consistent axis works; +X is
    // canonical.  The outward normal of such a face MUST have nx > 0.
    let find_seed = |orientation: &[Option<bool>]| -> Option<usize> {
        let mut best_x = f64::NEG_INFINITY;
        let mut best_fi: Option<usize> = None;
        for (fi, face) in face_list.iter().enumerate() {
            if orientation[fi].is_some() || face_normals[fi].is_none() {
                continue; // already visited or degenerate
            }
            for &vid in &face.vertices {
                let px = mesh.vertices.position(vid).x;
                if px > best_x {
                    best_x = px;
                    best_fi = Some(fi);
                }
            }
        }
        best_fi
    };

    // Outer loop handles disconnected patches (multiple connected components).
    loop {
        let Some(seed_fi) = find_seed(&orientation) else {
            break;
        };

        // Orient seed: outward if face normal has positive X component.
        let seed_normal = face_normals[seed_fi].unwrap(); // guaranteed non-None by find_seed
        let seed_is_outward = seed_normal.x >= 0.0;
        orientation[seed_fi] = Some(seed_is_outward);

        let mut queue: VecDeque<usize> = VecDeque::new();
        queue.push_back(seed_fi);

        while let Some(fi) = queue.pop_front() {
            let is_outward = orientation[fi].unwrap();
            let v = face_list[fi].vertices;

            // Inspect all three directed edges of this face.
            for k in 0..3 {
                let j = (k + 1) % 3;
                let va = v[k];
                let vb = v[j];

                // The manifold-adjacent face shares the REVERSE edge (vb→va).
                // Consistent adjacency → neighbour inherits same orientation.
                if let Some(&nfi) = half_edge.get(&(vb, va)) {
                    if orientation[nfi].is_none() && face_normals[nfi].is_some() {
                        orientation[nfi] = Some(is_outward);
                        queue.push_back(nfi);
                    }
                // Parallel edge (va→vb) in another face → winding flip.
                } else if let Some(&nfi) = half_edge.get(&(va, vb)) {
                    if orientation[nfi].is_none() && face_normals[nfi].is_some() {
                        orientation[nfi] = Some(!is_outward);
                        queue.push_back(nfi);
                    }
                }
            }
        }
    }

    // ── Step 4: count outward / inward / degenerate ──────────────────────────
    let mut outward = 0usize;
    let mut inward = 0usize;
    let mut degen = 0usize;

    for fi in 0..n_faces {
        if face_normals[fi].is_none() {
            degen += 1;
        } else {
            match orientation[fi] {
                Some(true) => outward += 1,
                Some(false) => inward += 1,
                None => degen += 1, // unreachable non-manifold fragment
            }
        }
    }

    // ── Step 5: face ↔ vertex-normal alignment statistics ───────────────────
    let mut asum: Real = 0.0;
    let mut acnt = 0usize;
    let mut amin: Real = 1.0;

    for (fi, face) in face_list.iter().enumerate() {
        let Some(face_n) = face_normals[fi] else {
            continue;
        };
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
