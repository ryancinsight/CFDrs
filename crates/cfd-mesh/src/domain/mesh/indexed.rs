use nalgebra::{Point3, Vector3};
// =========================================================================
// IndexedMesh<T> — watertight-first surface mesh, generic over precision
// =========================================================================

use crate::domain::core::index::{FaceId, RegionId, VertexId};
use crate::domain::core::scalar::Scalar;
use crate::domain::geometry::aabb::Aabb;
use crate::infrastructure::storage::attribute::AttributeStore;
use crate::infrastructure::storage::edge_store::EdgeStore;
use crate::infrastructure::storage::face_store::FaceStore;
use crate::infrastructure::storage::vertex_pool::VertexPool;
use crate::domain::topology::Cell;
use std::collections::HashMap;

/// A deduplicated, indexed triangle surface mesh — generic over scalar `T`.
///
/// | Type parameter | Precision | Tolerance |
/// |----------------|-----------|-----------|
/// | `f64` (default) | 64-bit | 1 nm |
/// | `f32`          | 32-bit | 10 µm (GPU staging) |
///
/// The default `T = f64` means all existing `IndexedMesh::new()` call-sites
/// continue to compile without any annotation.  New code may write
/// `IndexedMesh::<f32>::new()` to get single-precision geometry at zero
/// additional runtime cost.
///
/// Combines:
/// - [`VertexPool<T>`] — spatial-hash-deduplicated vertex storage
/// - `FaceStore` — indexed triangles with region tags
/// - `EdgeStore` — persistent adjacency (rebuilt on demand)
/// - `AttributeStore` — named per-face scalar channels
#[derive(Clone)]
pub struct IndexedMesh<T: Scalar = f64> {
    /// Deduplicated vertex positions and normals.
    pub vertices: VertexPool<T>,
    /// Indexed triangular faces.
    pub faces: FaceStore,
    /// Edge adjacency (lazily built from faces).
    edges: Option<EdgeStore>,
    /// Per-face scalar attributes.
    pub attributes: AttributeStore<FaceId>,
    /// Volumetric cells (for CFD support).
    pub cells: Vec<Cell>,
    /// Boundary patch names tagged by FaceId.
    pub boundary_labels: HashMap<FaceId, String>,
}

impl<T: Scalar> IndexedMesh<T> {
    /// Create an empty mesh with default millifluidic tolerances.
    pub fn new() -> Self {
        Self {
            vertices: VertexPool::default_millifluidic(),
            faces: FaceStore::new(),
            edges: None,
            attributes: AttributeStore::new(),
            cells: Vec::new(),
            boundary_labels: HashMap::new(),
        }
    }

    /// Create with explicit welding tolerance.
    pub fn with_tolerance(cell_size: T, tolerance: T) -> Self {
        Self {
            vertices: VertexPool::new(cell_size, tolerance),
            faces: FaceStore::new(),
            edges: None,
            attributes: AttributeStore::new(),
            cells: Vec::new(),
            boundary_labels: HashMap::new(),
        }
    }

    // ── Vertex operations ─────────────────────────────────────────────────

    /// Insert a vertex (deduplicated via spatial hash); returns its ID.
    pub fn add_vertex(&mut self, position: Point3<T>, normal: Vector3<T>) -> VertexId {
        self.edges = None;
        self.vertices.insert_or_weld(position, normal)
    }

    /// Insert a vertex by position only (zero normal).
    pub fn add_vertex_pos(&mut self, position: Point3<T>) -> VertexId {
        self.edges = None;
        self.vertices
            .insert_or_weld(position, Vector3::<T>::zeros())
    }

    /// Number of unique vertices.
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    // ── Face operations ───────────────────────────────────────────────────

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

    /// Flip the winding order of all faces (swap v1 <-> v2 on every triangle).
    ///
    /// Call this after building a mesh whose face-construction algorithm
    /// produces consistent *inward* normals, to obtain outward normals.
    pub fn flip_faces(&mut self) {
        self.edges = None;
        self.faces.iter_mut().for_each(|f| f.flip());
    }

    // ── Volumetric cell operations ────────────────────────────────────────

    /// Add a volumetric cell.
    pub fn add_cell(&mut self, c: Cell) {
        self.cells.push(c);
    }

    /// Number of cells.
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Immutable slice of all cells.
    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    // ── Boundary management ─────────────────────────────────────────────

    /// Label a face as a boundary with the given name.
    pub fn mark_boundary(&mut self, face_id: FaceId, label: impl Into<String>) {
        self.boundary_labels.insert(face_id, label.into());
    }

    /// Return the boundary label of a face, if any.
    pub fn boundary_label(&self, face_id: FaceId) -> Option<&str> {
        self.boundary_labels.get(&face_id).map(|s| s.as_str())
    }

    /// Return face IDs on the geometric boundary (faces belonging to exactly one cell).
    pub fn boundary_faces(&self) -> Vec<FaceId> {
        if self.cells.is_empty() {
            return self.faces.iter_enumerated().map(|(id, _)| id).collect();
        }
        let mut face_cell_count: HashMap<FaceId, usize> = HashMap::new();
        for cell in &self.cells {
            for &fv_idx in &cell.faces {
                // In IndexedMesh, Cell.faces holds FaceId cast as usize currently ? wait:
                let id = FaceId::from_usize(fv_idx);
                *face_cell_count.entry(id).or_insert(0) += 1;
            }
        }
        face_cell_count
            .into_iter()
            .filter(|&(_, count)| count == 1)
            .map(|(id, _)| id)
            .collect()
    }

    // ── Edge / adjacency access ───────────────────────────────────────────

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

    /// Immutable view of the edge store (may be stale).
    pub fn edges_ref(&self) -> Option<&EdgeStore> {
        self.edges.as_ref()
    }

    // ── Geometric queries ─────────────────────────────────────────────────

    /// Axis-aligned bounding box.
    pub fn bounding_box(&self) -> Aabb<T> {
        Aabb::from_points(self.vertices.positions())
    }

    /// Total surface area of all triangles.
    pub fn surface_area(&self) -> T {
        use crate::domain::geometry::measure;
        measure::total_surface_area(self.faces.iter_enumerated().map(|(_, f)| {
            (
                self.vertices.position(f.vertices[0]),
                self.vertices.position(f.vertices[1]),
                self.vertices.position(f.vertices[2]),
            )
        }))
    }

    /// Signed volume (positive for outward-oriented closed mesh).
    pub fn signed_volume(&self) -> T {
        use crate::domain::geometry::measure;
        measure::total_signed_volume(self.faces.iter_enumerated().map(|(_, f)| {
            (
                self.vertices.position(f.vertices[0]),
                self.vertices.position(f.vertices[1]),
                self.vertices.position(f.vertices[2]),
            )
        }))
    }

    // ── Validation ────────────────────────────────────────────────────────

    /// Check watertightness (rebuilds edges if needed).
    pub fn is_watertight(&mut self) -> bool {
        self.rebuild_edges();
        let edges = self.edges.as_ref().unwrap();
        let report = crate::application::watertight::check::check_watertight(
            &self.vertices,
            &self.faces,
            edges,
        );
        report.is_watertight
    }

    /// Run quality validation against default thresholds.
    pub fn quality_report(&self) -> crate::application::quality::validation::QualityReport {
        let validator = crate::application::quality::validation::MeshValidator::default();
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
        use crate::domain::geometry::normal::triangle_normal;

        let mut normal_sums: Vec<Vector3<T>> = vec![Vector3::<T>::zeros(); self.vertices.len()];
        let mut counts: Vec<usize> = vec![0; self.vertices.len()];

        for (_, face) in self.faces.iter_enumerated() {
            let a = self.vertices.position(face.vertices[0]);
            let b = self.vertices.position(face.vertices[1]);
            let c = self.vertices.position(face.vertices[2]);

            let face_normal = triangle_normal(a, b, c).unwrap_or_else(|| Vector3::<T>::z());

            for &vi in &face.vertices {
                normal_sums[vi.as_usize()] += face_normal;
                counts[vi.as_usize()] += 1;
            }
        }

        for (i, (sum, count)) in normal_sums.iter().zip(counts.iter()).enumerate() {
            if *count > 0 {
                let avg = *sum / <T as Scalar>::from_f64(*count as f64);
                let len = avg.norm();
                if len > <T as Scalar>::from_f64(1e-12) {
                    self.vertices.set_normal(VertexId::new(i as u32), avg / len);
                }
            }
        }
    }
}

impl<T: Scalar> Default for IndexedMesh<T> {
    fn default() -> Self {
        Self::new()
    }
}

// ── MeshBuilder<T> ────────────────────────────────────────────────────────────

/// Ergonomic builder for constructing an [`IndexedMesh<T>`].
pub struct MeshBuilder<T: Scalar = f64> {
    mesh: IndexedMesh<T>,
}

impl<T: Scalar> MeshBuilder<T> {
    /// Start building with default millifluidic tolerances.
    pub fn new() -> Self {
        Self {
            mesh: IndexedMesh::new(),
        }
    }

    /// Start building with custom tolerances.
    pub fn with_tolerance(cell_size: T, tolerance: T) -> Self {
        Self {
            mesh: IndexedMesh::with_tolerance(cell_size, tolerance),
        }
    }

    /// Add a vertex by position; returns its [`VertexId`].
    pub fn vertex(&mut self, pos: Point3<T>) -> VertexId {
        self.mesh.add_vertex_pos(pos)
    }

    /// Add a triangle from three vertex IDs.
    pub fn triangle(&mut self, v0: VertexId, v1: VertexId, v2: VertexId) -> FaceId {
        self.mesh.add_face(v0, v1, v2)
    }

    /// Add raw triangle soup — each triple is `(p0, p1, p2)`.
    pub fn add_triangle_soup(&mut self, triangles: &[(Point3<T>, Point3<T>, Point3<T>)]) {
        for (a, b, c) in triangles {
            let va = self.mesh.add_vertex_pos(*a);
            let vb = self.mesh.add_vertex_pos(*b);
            let vc = self.mesh.add_vertex_pos(*c);
            self.mesh.add_face(va, vb, vc);
        }
    }

    /// Finalise: build edges and return the mesh.
    pub fn build(mut self) -> IndexedMesh<T> {
        self.mesh.rebuild_edges();
        self.mesh
    }
}

impl<T: Scalar> Default for MeshBuilder<T> {
    fn default() -> Self {
        Self::new()
    }
}

