//! Mesh refinement module

use nalgebra::RealField;
use std::collections::HashMap;
use crate::topology::{Cell, Face, Vertex, ElementType};

pub mod criteria;

// Re-export main types
pub use criteria::*;

/// Refinement strategy trait
pub trait RefinementStrategy<T: RealField + Copy>: Send + Sync {
    /// Apply refinement to mesh
    fn refine(&self, mesh: &mut crate::mesh::Mesh<T>) -> Result<(), crate::error::MeshError>;

    /// Get strategy name
    fn name(&self) -> &'static str;
}

/// Uniform refinement strategy
pub struct UniformRefinement;

impl<T: RealField + Copy> RefinementStrategy<T> for UniformRefinement {
    fn refine(&self, mesh: &mut crate::mesh::Mesh<T>) -> Result<(), crate::error::MeshError> {
        if mesh.cell_count() == 0 {
            return Ok(());
        }

        // Check if we have only tetrahedra (simplification for now)
        for cell in mesh.cells() {
            if cell.element_type != ElementType::Tetrahedron {
                return Err(crate::error::MeshError::NotImplemented(
                    format!("Uniform refinement only implemented for Tetrahedra, found {:?}", cell.element_type)
                ));
            }
        }

        let mut new_mesh = crate::mesh::Mesh::new();
        
        // 1. Copy existing vertices and keep track of their new indices (identity map for now)
        // We assume new_mesh.vertices will be [old_vertices ..., new_midpoints ...]
        for v in mesh.vertices() {
            new_mesh.add_vertex(*v);
        }

        // Map from sorted edge vertex pair (v_min, v_max) to new midpoint vertex index
        let mut midpoint_cache: HashMap<(usize, usize), usize> = HashMap::new();
        
        // Map from sorted face vertex list to new face index (to deduplicate faces)
        let mut face_cache: HashMap<Vec<usize>, usize> = HashMap::new();

        // Helper to get or create midpoint
        let mut get_midpoint = |v1_idx: usize, v2_idx: usize, new_mesh: &mut crate::mesh::Mesh<T>, mesh: &crate::mesh::Mesh<T>| -> usize {
            let key = if v1_idx < v2_idx { (v1_idx, v2_idx) } else { (v2_idx, v1_idx) };
            if let Some(&idx) = midpoint_cache.get(&key) {
                idx
            } else {
                let v1 = mesh.vertex(v1_idx).unwrap();
                let v2 = mesh.vertex(v2_idx).unwrap();
                let center = nalgebra::center(&v1.position, &v2.position);
                // Create new vertex, inheriting partition info if both parents share it?
                // For now, just simple new vertex.
                let new_v = Vertex::new(center);
                let idx = new_mesh.add_vertex(new_v);
                midpoint_cache.insert(key, idx);
                idx
            }
        };

        // Helper to add face (deduplicated)
        let mut add_face_dedup = |v_indices: Vec<usize>, new_mesh: &mut crate::mesh::Mesh<T>| -> usize {
            let mut sorted_indices = v_indices.clone();
            sorted_indices.sort_unstable();
            
            if let Some(&idx) = face_cache.get(&sorted_indices) {
                idx
            } else {
                let face = Face::from_vertices(v_indices);
                let idx = new_mesh.add_face(face);
                face_cache.insert(sorted_indices, idx);
                idx
            }
        };

        // 2. Iterate over cells and refine
        for cell in mesh.cells() {
            // Assume Tetrahedron: vertices v0, v1, v2, v3
            // Faces: (0,1,2), (0,1,3), (0,2,3), (1,2,3)
            // Edges: (0,1), (0,2), (0,3), (1,2), (1,3), (2,3)
            
            // Get vertex indices
            // We need to extract them from the faces.
            // Mesh::element_vertices returns references, we need indices.
            // Let's assume standard ordering from the first face + opposite vertex?
            // Or just collect unique vertices.
            let mut v_indices = Vec::new();
            for &f_idx in &cell.faces {
                if let Some(face) = mesh.face(f_idx) {
                    for &v_idx in &face.vertices {
                        if !v_indices.contains(&v_idx) {
                            v_indices.push(v_idx);
                        }
                    }
                }
            }
            
            if v_indices.len() != 4 {
                 return Err(crate::error::MeshError::InvalidMesh(
                    format!("Tetrahedron must have 4 vertices, found {}", v_indices.len())
                ));
            }
            
            // Sort to ensure deterministic processing if needed, but we need to know topology.
            // Actually, for a tet, any permutation is a valid tet, but we need to know which edges exist.
            // In a full mesh, all pairs of vertices in a tet are edges.
            let v0 = v_indices[0];
            let v1 = v_indices[1];
            let v2 = v_indices[2];
            let v3 = v_indices[3];

            // Get midpoints
            let m01 = get_midpoint(v0, v1, &mut new_mesh, mesh);
            let m02 = get_midpoint(v0, v2, &mut new_mesh, mesh);
            let m03 = get_midpoint(v0, v3, &mut new_mesh, mesh);
            let m12 = get_midpoint(v1, v2, &mut new_mesh, mesh);
            let m13 = get_midpoint(v1, v3, &mut new_mesh, mesh);
            let m23 = get_midpoint(v2, v3, &mut new_mesh, mesh);

            // Define the 8 new tets (vertices)
            // 4 Corners
            let t0 = vec![v0, m01, m02, m03];
            let t1 = vec![v1, m01, m12, m13];
            let t2 = vec![v2, m02, m12, m23];
            let t3 = vec![v3, m03, m13, m23];

            // 4 Inner (Octahedron split)
            // Diagonal choice: (m01, m23) -> indices (m01, m23)
            // Split octahedron (m01, m02, m03, m12, m13, m23) along diagonal m01-m23
            // The 4 tets around this diagonal:
            // T4: m01, m23, m02, m12
            let t4 = vec![m01, m23, m02, m12];
            // T5: m01, m23, m12, m13
            let t5 = vec![m01, m23, m12, m13];
            // T6: m01, m23, m13, m03
            let t6 = vec![m01, m23, m13, m03];
            // T7: m01, m23, m03, m02
            let t7 = vec![m01, m23, m03, m02];

            let new_tets = vec![t0, t1, t2, t3, t4, t5, t6, t7];

            for tet_verts in new_tets {
                // Create faces for new tet
                let f0 = add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[2]], &mut new_mesh);
                let f1 = add_face_dedup(vec![tet_verts[0], tet_verts[1], tet_verts[3]], &mut new_mesh);
                let f2 = add_face_dedup(vec![tet_verts[0], tet_verts[2], tet_verts[3]], &mut new_mesh);
                let f3 = add_face_dedup(vec![tet_verts[1], tet_verts[2], tet_verts[3]], &mut new_mesh);

                new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
            }
        }

        // Replace mesh
        *mesh = new_mesh;

        Ok(())
    }

    fn name(&self) -> &'static str {
        "Uniform"
    }
}

/// Adaptive refinement strategy based on a per-cell error indicator.
///
/// # Algorithm (Dörfler marking)
///
/// Given a per-cell error indicator η_K:
/// 1. Sort cells by decreasing η_K.
/// 2. Mark cells for refinement until the cumulative squared error
///    exceeds `theta * Σ η_K²` (Dörfler's criterion, θ ∈ (0,1]).
///    With θ = 1.0 this degenerates to refining every cell above `threshold`.
/// 3. Apply the `UniformRefinement` split (1:8 tet subdivision) to marked cells.
///
/// # References
///
/// - Dörfler, W. (1996) "A convergent adaptive algorithm for Poisson's equation",
///   SIAM J. Numer. Anal.
/// - Verfürth, R. (1996) "A Review of A Posteriori Error Estimation and Adaptive
///   Mesh-Refinement Techniques", Wiley-Teubner.
pub struct AdaptiveRefinement<T: RealField + Copy> {
    /// Per-cell error indicator field.
    /// Length must equal the number of cells in the mesh when `refine` is called.
    /// Commonly this is the gradient jump estimator or recovery-based estimator.
    pub error_field: Vec<T>,
    /// Absolute refinement threshold η_min.
    /// Cells with η_K < threshold are never refined regardless of Dörfler marking.
    pub threshold: T,
    /// Dörfler marking parameter θ ∈ (0, 1].
    /// θ = 0.5 refines about 50 % of the total error; θ = 1.0 refines all cells
    /// above `threshold`.
    pub theta: T,
}

impl<T: RealField + Copy> AdaptiveRefinement<T> {
    /// Convenience constructor.
    pub fn new(error_field: Vec<T>, threshold: T) -> Self {
        Self {
            error_field,
            threshold,
            theta: T::from_f64(0.5).unwrap_or_else(T::one),
        }
    }

    /// Construct with explicit Dörfler parameter.
    pub fn with_theta(error_field: Vec<T>, threshold: T, theta: T) -> Self {
        Self { error_field, threshold, theta }
    }
}

impl<T: RealField + Copy> RefinementStrategy<T> for AdaptiveRefinement<T> {
    fn refine(&self, mesh: &mut crate::mesh::Mesh<T>) -> Result<(), crate::error::MeshError> {
        let n_cells = mesh.cell_count();
        if n_cells == 0 {
            return Ok(());
        }

        if self.error_field.len() != n_cells {
            return Err(crate::error::MeshError::InvalidMesh(format!(
                "Error field length ({}) does not match cell count ({})",
                self.error_field.len(),
                n_cells,
            )));
        }

        // Only tetrahedral meshes are supported (same constraint as UniformRefinement).
        for cell in mesh.cells() {
            if cell.element_type != ElementType::Tetrahedron {
                return Err(crate::error::MeshError::NotImplemented(format!(
                    "Adaptive refinement only implemented for Tetrahedra, found {:?}",
                    cell.element_type
                )));
            }
        }

        // ------------------------------------------------------------------
        // Step 1: Dörfler marking
        // ------------------------------------------------------------------
        // Compute total squared error
        let total_eta_sq: T = self
            .error_field
            .iter()
            .fold(T::zero(), |acc, &e| acc + e * e);

        if total_eta_sq <= T::zero() {
            return Ok(()); // nothing to refine
        }

        // Sort cell indices by descending error
        let mut sorted_indices: Vec<usize> = (0..n_cells).collect();
        sorted_indices.sort_by(|&a, &b| {
            self.error_field[b]
                .partial_cmp(&self.error_field[a])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        let target = self.theta * total_eta_sq;
        let mut cumulative = T::zero();
        let mut marked = vec![false; n_cells];

        for &idx in &sorted_indices {
            let eta = self.error_field[idx];
            if eta < self.threshold {
                break; // remaining cells are below absolute threshold
            }
            marked[idx] = true;
            cumulative += eta * eta;
            if cumulative >= target {
                break;
            }
        }

        // ------------------------------------------------------------------
        // Step 2: Build refined mesh — marked cells get 1:8 split,
        //         unmarked cells are copied verbatim.
        // ------------------------------------------------------------------
        let mut new_mesh = crate::mesh::Mesh::new();

        // Copy all existing vertices
        for v in mesh.vertices() {
            new_mesh.add_vertex(*v);
        }

        let mut midpoint_cache: HashMap<(usize, usize), usize> = HashMap::new();
        let mut face_cache: HashMap<Vec<usize>, usize> = HashMap::new();

        let mesh_ref = &*mesh; // immutable borrow for helpers

        let mut get_midpoint = |v1: usize, v2: usize, nm: &mut crate::mesh::Mesh<T>| -> usize {
            let key = if v1 < v2 { (v1, v2) } else { (v2, v1) };
            if let Some(&idx) = midpoint_cache.get(&key) {
                idx
            } else {
                let p1 = mesh_ref.vertex(v1).unwrap().position;
                let p2 = mesh_ref.vertex(v2).unwrap().position;
                let center = nalgebra::center(&p1, &p2);
                let idx = nm.add_vertex(Vertex::new(center));
                midpoint_cache.insert(key, idx);
                idx
            }
        };

        let mut add_face_dedup = |vs: Vec<usize>, nm: &mut crate::mesh::Mesh<T>| -> usize {
            let mut sorted = vs.clone();
            sorted.sort_unstable();
            if let Some(&idx) = face_cache.get(&sorted) {
                idx
            } else {
                let f = Face::from_vertices(vs);
                let idx = nm.add_face(f);
                face_cache.insert(sorted, idx);
                idx
            }
        };

        for (cell_idx, cell) in mesh.cells().iter().enumerate() {
            // Collect vertex indices for this tet
            let mut v_indices = Vec::new();
            for &f_idx in &cell.faces {
                if let Some(face) = mesh.face(f_idx) {
                    for &v in &face.vertices {
                        if !v_indices.contains(&v) {
                            v_indices.push(v);
                        }
                    }
                }
            }
            if v_indices.len() != 4 {
                return Err(crate::error::MeshError::InvalidMesh(format!(
                    "Tet must have 4 vertices, found {}",
                    v_indices.len()
                )));
            }

            if !marked[cell_idx] {
                // Copy cell unchanged
                let f0 = add_face_dedup(vec![v_indices[0], v_indices[1], v_indices[2]], &mut new_mesh);
                let f1 = add_face_dedup(vec![v_indices[0], v_indices[1], v_indices[3]], &mut new_mesh);
                let f2 = add_face_dedup(vec![v_indices[0], v_indices[2], v_indices[3]], &mut new_mesh);
                let f3 = add_face_dedup(vec![v_indices[1], v_indices[2], v_indices[3]], &mut new_mesh);
                new_mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
            } else {
                // 1:8 uniform refinement of this tet
                let (v0, v1, v2, v3) = (v_indices[0], v_indices[1], v_indices[2], v_indices[3]);
                let m01 = get_midpoint(v0, v1, &mut new_mesh);
                let m02 = get_midpoint(v0, v2, &mut new_mesh);
                let m03 = get_midpoint(v0, v3, &mut new_mesh);
                let m12 = get_midpoint(v1, v2, &mut new_mesh);
                let m13 = get_midpoint(v1, v3, &mut new_mesh);
                let m23 = get_midpoint(v2, v3, &mut new_mesh);

                // 4 corner tets + 4 interior tets (octahedron via m01-m23 diagonal)
                let new_tets: Vec<[usize; 4]> = vec![
                    [v0, m01, m02, m03],
                    [v1, m01, m12, m13],
                    [v2, m02, m12, m23],
                    [v3, m03, m13, m23],
                    [m01, m23, m02, m12],
                    [m01, m23, m12, m13],
                    [m01, m23, m13, m03],
                    [m01, m23, m03, m02],
                ];

                for t in new_tets {
                    let fa = add_face_dedup(vec![t[0], t[1], t[2]], &mut new_mesh);
                    let fb = add_face_dedup(vec![t[0], t[1], t[3]], &mut new_mesh);
                    let fc = add_face_dedup(vec![t[0], t[2], t[3]], &mut new_mesh);
                    let fd = add_face_dedup(vec![t[1], t[2], t[3]], &mut new_mesh);
                    new_mesh.add_cell(Cell::tetrahedron(fa, fb, fc, fd));
                }
            }
        }

        *mesh = new_mesh;
        Ok(())
    }

    fn name(&self) -> &'static str {
        "Adaptive"
    }
}
