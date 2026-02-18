//! Generic mesh type for CFD solvers.
//!
//! This module provides two mesh types:
//!
//! - **`Mesh<T>`** — legacy FEM/FVM mesh with typed vertices, faces, and cells.
//! - **`IndexedMesh`** — the new indexed-triangle surface mesh combining
//!   `VertexPool` (spatial-hash dedup), `FaceStore`, `EdgeStore`, and
//!   `AttributeStore` for watertight CFD geometry.

use nalgebra::{Point3, RealField};
use std::collections::HashMap;

use crate::topology::{Cell, ElementType, Face, Vertex};

/// Statistics about a mesh.
#[derive(Clone, Debug, Default)]
pub struct MeshStatistics {
    /// Number of vertices.
    pub vertex_count: usize,
    /// Number of cells.
    pub cell_count: usize,
    /// Number of faces on the boundary.
    pub boundary_face_count: usize,
}

/// A generic CFD mesh with vertices, faces, and cells.
///
/// `T` is the floating-point scalar type (typically `f64`).
#[derive(Clone, Debug)]
pub struct Mesh<T: Copy + RealField> {
    // ---- high-level typed arrays -----------------------------------------
    vertices: Vec<Vertex<T>>,
    faces: Vec<Face>,
    cells: Vec<Cell>,
    boundary_labels: HashMap<usize, String>,

    // ---- low-level raw arrays (for tests and legacy access) ---------------
    /// Raw node positions (mirrors `vertices[i].position`).
    pub nodes: Vec<Point3<T>>,
    /// Raw element connectivity in terms of vertex indices.
    pub elements: Vec<Vec<usize>>,
}

impl<T: Copy + RealField> Default for Mesh<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: Copy + RealField> Mesh<T> {
    /// Create an empty mesh.
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            faces: Vec::new(),
            cells: Vec::new(),
            boundary_labels: HashMap::new(),
            nodes: Vec::new(),
            elements: Vec::new(),
        }
    }

    // ---- vertex access ---------------------------------------------------

    /// Append a vertex; returns its index.
    pub fn add_vertex(&mut self, v: Vertex<T>) -> usize {
        let idx = self.vertices.len();
        self.nodes.push(v.position);
        self.vertices.push(v);
        idx
    }

    /// Number of vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of nodes (alias for `vertex_count`).
    pub fn num_nodes(&self) -> usize {
        self.vertices.len()
    }

    /// Immutable slice of all vertices.
    pub fn vertices(&self) -> &[Vertex<T>] {
        &self.vertices
    }

    /// Get a vertex by index.
    pub fn vertex(&self, idx: usize) -> Option<&Vertex<T>> {
        self.vertices.get(idx)
    }

    // ---- face access -----------------------------------------------------

    /// Append a face; returns its index (used for building cells).
    pub fn add_face(&mut self, f: Face) -> usize {
        let idx = self.faces.len();
        self.faces.push(f);
        idx
    }

    /// Number of faces.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Immutable slice of all faces.
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    /// Get a face by index.
    pub fn face(&self, idx: usize) -> Option<&Face> {
        self.faces.get(idx)
    }

    // ---- cell access -----------------------------------------------------

    /// Append a cell.
    pub fn add_cell(&mut self, c: Cell) {
        // also mirror into raw elements
        let verts: Vec<usize> = c.faces.iter()
            .flat_map(|&fi| self.faces.get(fi).map(|f| f.vertices.clone()).unwrap_or_default())
            .collect::<std::collections::BTreeSet<_>>()
            .into_iter()
            .collect();
        self.elements.push(verts);
        self.cells.push(c);
    }

    /// Number of cells.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Number of elements (alias for `cell_count`).
    pub fn num_elements(&self) -> usize {
        self.cells.len()
    }

    /// Immutable slice of all cells.
    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    /// Get a cell by index.
    pub fn cell(&self, idx: usize) -> Option<&Cell> {
        self.cells.get(idx)
    }

    /// Return the number of vertices per element (4 for tetrahedra).
    pub fn nodes_per_element(&self) -> usize {
        self.cells.first().map(|c| match c.element_type {
            ElementType::Tetrahedron => 4,
            ElementType::Hexahedron => 8,
            ElementType::Triangle => 3,
            ElementType::Quadrilateral => 4,
            ElementType::Wedge => 6,
            ElementType::Pyramid => 5,
        }).unwrap_or(4)
    }

    // ---- boundary management --------------------------------------------

    /// Label a face as a boundary with the given name.
    pub fn mark_boundary(&mut self, face_idx: usize, label: String) {
        self.boundary_labels.insert(face_idx, label);
    }

    /// Return the boundary label of a face, if any.
    pub fn boundary_label(&self, face_idx: usize) -> Option<&str> {
        self.boundary_labels.get(&face_idx).map(String::as_str)
    }

    /// Return all face indices that have been labeled as boundaries.
    pub fn marked_boundary_faces(&self) -> Vec<usize> {
        self.boundary_labels.keys().copied().collect()
    }

    /// Return face indices on the geometric boundary (faces belonging to
    /// exactly one cell). This is computed from face-cell adjacency.
    pub fn boundary_faces(&self) -> Vec<usize> {
        let mut face_cell_count: HashMap<usize, usize> = HashMap::new();
        for cell in &self.cells {
            for &fi in &cell.faces {
                *face_cell_count.entry(fi).or_insert(0) += 1;
            }
        }
        // Faces referenced by exactly one cell are on the boundary.
        // If no cells exist, all faces are boundary.
        if self.cells.is_empty() {
            return (0..self.faces.len()).collect();
        }
        face_cell_count
            .into_iter()
            .filter(|&(_, count)| count == 1)
            .map(|(fi, _)| fi)
            .collect()
    }

    // ---- statistics ------------------------------------------------------

    /// Compute mesh statistics.
    pub fn statistics(&self) -> MeshStatistics {
        let bf = self.boundary_faces();
        MeshStatistics {
            vertex_count: self.vertices.len(),
            cell_count: self.cells.len(),
            boundary_face_count: bf.len(),
        }
    }
}

// =========================================================================
// IndexedMesh — the new watertight-first surface mesh
// =========================================================================

use crate::core::index::{VertexId, FaceId, RegionId};
use crate::core::scalar::{Real, Point3r, Vector3r, TOLERANCE};
use crate::storage::vertex_pool::VertexPool;
use crate::storage::face_store::{FaceStore, FaceData};
use crate::storage::edge_store::EdgeStore;
use crate::storage::attribute::AttributeStore;
use crate::geometry::aabb::Aabb;

/// A deduplicated, indexed triangle surface mesh.
///
/// Combines:
/// - `VertexPool` — spatial-hash-deduplicated vertex storage
/// - `FaceStore` — indexed triangles with region tags
/// - `EdgeStore` — persistent adjacency (rebuilt on demand)
/// - `AttributeStore` — named per-face scalar channels
///
/// Use `MeshBuilder` for ergonomic construction.
pub struct IndexedMesh {
    /// Shared, deduplicated vertex positions + normals.
    pub vertices: VertexPool,
    /// Indexed triangular faces.
    pub faces: FaceStore,
    /// Edge adjacency (lazily built from faces).
    edges: Option<EdgeStore>,
    /// Per-face scalar attributes.
    pub attributes: AttributeStore<FaceId>,
}

impl IndexedMesh {
    /// Create an empty `IndexedMesh` with default millifluidic tolerances.
    pub fn new() -> Self {
        Self {
            vertices: VertexPool::default_millifluidic(),
            faces: FaceStore::new(),
            edges: None,
            attributes: AttributeStore::new(),
        }
    }

    /// Create with a custom welding tolerance.
    pub fn with_tolerance(cell_size: Real, tolerance: Real) -> Self {
        Self {
            vertices: VertexPool::new(cell_size, tolerance),
            faces: FaceStore::new(),
            edges: None,
            attributes: AttributeStore::new(),
        }
    }

    // ── Vertex operations ─────────────────────────────────────

    /// Insert a vertex (deduplicated via spatial hash).
    pub fn add_vertex(&mut self, position: Point3r, normal: Vector3r) -> VertexId {
        self.edges = None; // invalidate edge cache
        self.vertices.insert_or_weld(position, normal)
    }

    /// Insert a vertex by position only (zero normal).
    pub fn add_vertex_pos(&mut self, position: Point3r) -> VertexId {
        self.edges = None;
        self.vertices.insert_or_weld(position, Vector3r::zeros())
    }

    /// Number of unique vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    // ── Face operations ───────────────────────────────────────

    /// Add a triangle face from three vertex IDs.
    pub fn add_face(&mut self, v0: VertexId, v1: VertexId, v2: VertexId) -> FaceId {
        self.edges = None;
        self.faces.add_triangle(v0, v1, v2)
    }

    /// Add a triangle face with a region tag.
    pub fn add_face_with_region(
        &mut self,
        v0: VertexId,
        v1: VertexId,
        v2: VertexId,
        region: RegionId,
    ) -> FaceId {
        self.edges = None;
        self.faces.add_triangle_with_region(v0, v1, v2, region)
    }

    /// Number of faces.
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    // ── Edge / adjacency access ───────────────────────────────

    /// Get (or lazily build) the edge store.
    pub fn edges(&mut self) -> &EdgeStore {
        if self.edges.is_none() {
        self.edges = Some(EdgeStore::from_face_store(&self.faces));
        }
        self.edges.as_ref().unwrap()
    }

    /// Force rebuild of edge adjacency.
    pub fn rebuild_edges(&mut self) {
        self.edges = Some(EdgeStore::from_face_store(&self.faces));
    }

    /// Get existing edges (immutable).
    pub fn edges_ref(&self) -> Option<&EdgeStore> {
        self.edges.as_ref()
    }

    // ── Geometric queries ─────────────────────────────────────

    /// Axis-aligned bounding box.
    pub fn bounding_box(&self) -> Aabb {
        Aabb::from_points(self.vertices.positions())
    }

    /// Total surface area.
    pub fn surface_area(&self) -> Real {
        use crate::geometry::measure;
        measure::total_surface_area(
            self.faces.iter_enumerated().map(|(_, f)| {
                (
                    self.vertices.position(f.vertices[0]),
                    self.vertices.position(f.vertices[1]),
                    self.vertices.position(f.vertices[2]),
                )
            }),
        )
    }

    /// Signed volume (positive for outward-oriented closed mesh).
    pub fn signed_volume(&self) -> Real {
        use crate::geometry::measure;
        measure::total_signed_volume(
            self.faces.iter_enumerated().map(|(_, f)| {
                (
                    self.vertices.position(f.vertices[0]),
                    self.vertices.position(f.vertices[1]),
                    self.vertices.position(f.vertices[2]),
                )
            }),
        )
    }

    // ── Validation ────────────────────────────────────────────

    /// Check watertightness (rebuilds edges if needed).
    pub fn is_watertight(&mut self) -> bool {
        self.rebuild_edges();
        let edges = self.edges.as_ref().unwrap();
        let report = crate::watertight::check::check_watertight(
            &self.vertices,
            &self.faces,
            edges,
        );
        report.is_watertight
    }

    /// Run quality validation against default thresholds.
    pub fn quality_report(&self) -> crate::quality::validation::QualityReport {
        let validator = crate::quality::validation::MeshValidator::default();
        validator.validate(&self.faces, &self.vertices)
    }

    // ── Normal Recomputation ──────────────────────────────────

    /// Recompute all vertex normals from face geometry.
    ///
    /// After CSG operations, face winding may have changed, but the vertex
    /// normals stored in the pool are the original normals. This method
    /// recalculates normals based on the current face winding, averaging
    /// contributions from all faces that share each vertex.
    ///
    /// For each face, the normal is computed from the cross product:
    /// `n = normalize((v1 - v0) × (v2 - v0))`
    ///
    /// Each vertex's normal is the average of all face normals that use it.
    pub fn recompute_normals(&mut self) {
        use crate::geometry::normal::triangle_normal;

        // Accumulate normals for each vertex
        let mut normal_sums: Vec<Vector3r> = vec![Vector3r::zeros(); self.vertices.len()];
        let mut counts: Vec<usize> = vec![0; self.vertices.len()];

        for (_, face) in self.faces.iter_enumerated() {
            let a = self.vertices.position(face.vertices[0]);
            let b = self.vertices.position(face.vertices[1]);
            let c = self.vertices.position(face.vertices[2]);

            let face_normal = triangle_normal(&a, &b, &c).unwrap_or_else(|| Vector3r::z());

            for &vi in &face.vertices {
                normal_sums[vi.as_usize()] += face_normal;
                counts[vi.as_usize()] += 1;
            }
        }

        // Update vertex normals
        for (i, (sum, count)) in normal_sums.iter().zip(counts.iter()).enumerate() {
            if *count > 0 {
                let avg = sum / (*count as Real);
                let len = avg.norm();
                if len > 1e-12 {
                    self.vertices.set_normal(VertexId::new(i as u32), avg / len);
                }
            }
        }
    }
}

impl Default for IndexedMesh {
    fn default() -> Self {
        Self::new()
    }
}

// ── Builder ───────────────────────────────────────────────────

/// Ergonomic builder for constructing an `IndexedMesh`.
pub struct MeshBuilder {
    mesh: IndexedMesh,
}

impl MeshBuilder {
    /// Start building with default settings.
    pub fn new() -> Self {
        Self {
            mesh: IndexedMesh::new(),
        }
    }

    /// Start building with custom tolerance.
    pub fn with_tolerance(cell_size: Real, tolerance: Real) -> Self {
        Self {
            mesh: IndexedMesh::with_tolerance(cell_size, tolerance),
        }
    }

    /// Add a vertex by position (returns the builder for chaining).
    pub fn vertex(&mut self, pos: Point3r) -> VertexId {
        self.mesh.add_vertex_pos(pos)
    }

    /// Add a triangle from three vertex IDs.
    pub fn triangle(&mut self, v0: VertexId, v1: VertexId, v2: VertexId) -> FaceId {
        self.mesh.add_face(v0, v1, v2)
    }

    /// Add raw triangle soup — each triple is (p0, p1, p2).
    pub fn add_triangle_soup(&mut self, triangles: &[(Point3r, Point3r, Point3r)]) {
        for (a, b, c) in triangles {
            let va = self.mesh.add_vertex_pos(*a);
            let vb = self.mesh.add_vertex_pos(*b);
            let vc = self.mesh.add_vertex_pos(*c);
            self.mesh.add_face(va, vb, vc);
        }
    }

    /// Finalize the mesh (builds edges).
    pub fn build(mut self) -> IndexedMesh {
        self.mesh.rebuild_edges();
        self.mesh
    }
}

impl Default for MeshBuilder {
    fn default() -> Self {
        Self::new()
    }
}
