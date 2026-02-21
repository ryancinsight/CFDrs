//! # Mesh Types
//!
//! This module provides three mesh types:
//!
//! - **[`HalfEdgeMesh<'id>`]** — the new state-of-the-art mesh backed by a
//!   GhostCell-permissioned half-edge topology kernel with SlotMap generational
//!   keys. Use [`with_mesh`] as the entry point.
//! - **`Mesh<T>`** — legacy FEM/FVM mesh with typed vertices, faces, and cells.
//!   Kept for downstream compatibility with `cfd-3d` geometry builders.
//! - **[`IndexedMesh`]** — the watertight-first surface mesh combining
//!   `VertexPool` (spatial-hash dedup), `FaceStore`, `EdgeStore`, and
//!   `AttributeStore`. Also kept for backward compatibility.
//!
//! ## Quick Start — New API
//!
//! ```rust,ignore
//! use cfd_mesh::mesh::with_mesh;
//!
//! let result = with_mesh(|mut mesh, mut token| {
//!     let a = mesh.add_vertex([0.0, 0.0, 0.0], &token);
//!     let b = mesh.add_vertex([1.0, 0.0, 0.0], &token);
//!     let c = mesh.add_vertex([0.0, 1.0, 0.0], &token);
//!     mesh.add_triangle(a, b, c, &mut token).expect("valid triangle");
//!     mesh.vertex_count()
//! });
//! assert_eq!(result, 3);
//! ```
//!
//! ## Architecture
//!
//! ```text
//! with_mesh(|mesh, token| {
//!   │
//!   ├── GhostToken<'id>  — the single permission key for the entire mesh
//!   │                      &token  → read access to any GhostCell<'id, _>
//!   │                      &mut token → write access (exclusive)
//!   │
//!   └── HalfEdgeMesh<'id>
//!       ├── GhostSlotPool<'id, VertexKey,   VertexData>
//!       ├── GhostSlotPool<'id, HalfEdgeKey, HalfEdgeData>
//!       ├── GhostSlotPool<'id, FaceKey,     FaceData>
//!       └── SlotMap<PatchKey, BoundaryPatch>   (no aliasing → no GhostCell)
//! })

use nalgebra::{Point3, Vector3, RealField};
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
#[deprecated(note = "use `IndexedMesh` for surface meshes; `Mesh<T>` is retained only for volume/FEM tools (grid.rs, hex_to_tet.rs)")]
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

// The impl blocks below intentionally work with the deprecated Mesh<T>; suppress the
// self-referential warning so it does not pollute the compiler output for call-sites.
#[allow(deprecated)]
impl<T: Copy + RealField> Default for Mesh<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[allow(deprecated)]
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

    /// Return the axis-aligned bounding box as `(min_point, max_point)`.
    ///
    /// Returns `(Point3::origin(), Point3::origin())` for an empty mesh.
    pub fn bounds(&self) -> (Point3<T>, Point3<T>) {
        if self.vertices.is_empty() {
            return (Point3::origin(), Point3::origin());
        }
        let mut lo = self.vertices[0].position;
        let mut hi = lo;
        for v in &self.vertices[1..] {
            let p = v.position;
            for i in 0..3 {
                if p[i] < lo[i] { lo[i] = p[i]; }
                if p[i] > hi[i] { hi[i] = p[i]; }
            }
        }
        (lo, hi)
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

    /// Cast the mesh scalar type to another type `U`.
    pub fn cast<U: Copy + RealField + num_traits::FromPrimitive + 'static>(self) -> Mesh<U>
    where
        T: num_traits::ToPrimitive,
    {
        let convert_point = |p: Point3<T>| -> Point3<U> {
            Point3::new(
                U::from_f64(p.x.to_f64().unwrap()).unwrap(),
                U::from_f64(p.y.to_f64().unwrap()).unwrap(),
                U::from_f64(p.z.to_f64().unwrap()).unwrap(),
            )
        };

        let vertices = self
            .vertices
            .into_iter()
            .map(|v| Vertex {
                position: convert_point(v.position),
            })
            .collect();
        let faces = self.faces; // Faces are integer indices, no conversion needed
        let cells = self.cells; // Cells are integer indices, no conversion needed
        let boundary_labels = self.boundary_labels;
        let nodes = self.nodes.into_iter().map(convert_point).collect();
        let elements = self.elements;

        Mesh {
            vertices,
            faces,
            cells,
            boundary_labels,
            nodes,
            elements,
        }
    }
}

#[allow(deprecated)]
impl<T: Scalar> From<IndexedMesh<T>> for Mesh<T> {
    fn from(indexed: IndexedMesh<T>) -> Self {
        let mut mesh = Mesh::new();

        // Transfer vertices
        for i in 0..indexed.vertex_count() {
            let p = indexed.vertices.position(VertexId::from_usize(i));
            mesh.add_vertex(Vertex { position: *p });
        }

        // Transfer faces
        for i in 0..indexed.face_count() {
            let f = indexed.faces.get(FaceId::from_usize(i));
            let vertices: Vec<usize> = f.vertices.iter().map(|v| v.as_usize()).collect();
            // Note: IndexedMesh faces don't store normal directly in Face struct usually,
            // but for Mesh<T> we just need connectivity.
            mesh.add_face(Face { vertices });
        }

        // IndexedMesh is surface-only, so cells are empty.
        mesh
    }
}

// =========================================================================
// IndexedMesh<T> — watertight-first surface mesh, generic over precision
// =========================================================================

use crate::core::index::{VertexId, FaceId, RegionId};
use crate::core::scalar::{Real, Scalar};
use crate::storage::vertex_pool::VertexPool;
use crate::storage::face_store::FaceStore;
use crate::storage::edge_store::EdgeStore;
use crate::storage::attribute::AttributeStore;
use crate::geometry::aabb::Aabb;

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
pub struct IndexedMesh<T: Scalar = f64> {
    /// Deduplicated vertex positions and normals.
    pub vertices: VertexPool<T>,
    /// Indexed triangular faces.
    pub faces: FaceStore,
    /// Edge adjacency (lazily built from faces).
    edges: Option<EdgeStore>,
    /// Per-face scalar attributes.
    pub attributes: AttributeStore<FaceId>,
}

impl<T: Scalar> IndexedMesh<T> {
    /// Create an empty mesh with default millifluidic tolerances.
    pub fn new() -> Self {
        Self {
            vertices:   VertexPool::default_millifluidic(),
            faces:      FaceStore::new(),
            edges:      None,
            attributes: AttributeStore::new(),
        }
    }

    /// Create with explicit welding tolerance.
    pub fn with_tolerance(cell_size: T, tolerance: T) -> Self {
        Self {
            vertices:   VertexPool::new(cell_size, tolerance),
            faces:      FaceStore::new(),
            edges:      None,
            attributes: AttributeStore::new(),
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
        self.vertices.insert_or_weld(position, Vector3::<T>::zeros())
    }

    /// Number of unique vertices.
    pub fn vertex_count(&self) -> usize { self.vertices.len() }

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
    pub fn face_count(&self) -> usize { self.faces.len() }

    /// Flip the winding order of all faces (swap v1 <-> v2 on every triangle).
    ///
    /// Call this after building a mesh whose face-construction algorithm
    /// produces consistent *inward* normals, to obtain outward normals.
    pub fn flip_faces(&mut self) {
        self.edges = None;
        self.faces.iter_mut().for_each(|f| f.flip());
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
    pub fn edges_ref(&self) -> Option<&EdgeStore> { self.edges.as_ref() }

    // ── Geometric queries ─────────────────────────────────────────────────

    /// Axis-aligned bounding box.
    pub fn bounding_box(&self) -> Aabb<T> {
        Aabb::from_points(self.vertices.positions())
    }

    /// Total surface area of all triangles.
    pub fn surface_area(&self) -> T {
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
    pub fn signed_volume(&self) -> T {
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

    // ── Validation ────────────────────────────────────────────────────────

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

        let mut normal_sums: Vec<Vector3<T>> = vec![Vector3::<T>::zeros(); self.vertices.len()];
        let mut counts: Vec<usize> = vec![0; self.vertices.len()];

        for (_, face) in self.faces.iter_enumerated() {
            let a = self.vertices.position(face.vertices[0]);
            let b = self.vertices.position(face.vertices[1]);
            let c = self.vertices.position(face.vertices[2]);

            let face_normal = triangle_normal(a, b, c)
                .unwrap_or_else(|| Vector3::<T>::z());

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
    fn default() -> Self { Self::new() }
}

// ── MeshBuilder<T> ────────────────────────────────────────────────────────────

/// Ergonomic builder for constructing an [`IndexedMesh<T>`].
pub struct MeshBuilder<T: Scalar = f64> {
    mesh: IndexedMesh<T>,
}

impl<T: Scalar> MeshBuilder<T> {
    /// Start building with default millifluidic tolerances.
    pub fn new() -> Self {
        Self { mesh: IndexedMesh::new() }
    }

    /// Start building with custom tolerances.
    pub fn with_tolerance(cell_size: T, tolerance: T) -> Self {
        Self { mesh: IndexedMesh::with_tolerance(cell_size, tolerance) }
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
    fn default() -> Self { Self::new() }
}

// =============================================================================
// HalfEdgeMesh<'id> — the new GhostCell-permissioned half-edge mesh
// =============================================================================

use slotmap::SlotMap;
use crate::permission::{GhostToken, GhostCell};
use crate::storage::slotmap_pool::GhostSlotPool;
use crate::core::index::{VertexKey, HalfEdgeKey, FaceKey, PatchKey};
// Real is already imported above via `use crate::core::scalar::{Real, ...}`
use crate::topology::halfedge::{
    BoundaryPatch as HeBoundaryPatch,
    FaceData as HeFaceData,
    HalfEdgeData,
    PatchType,
    VertexData,
};
use crate::core::error::{MeshError, MeshResult};

/// The state-of-the-art half-edge mesh with GhostCell-permissioned access.
///
/// All mesh data (vertex positions, half-edge connectivity, face normals) lives
/// inside [`GhostCell`] wrappers inside [`SlotMap`]s. Reading or writing any
/// datum requires the matching [`GhostToken<'id>`].
///
/// # Entry Point
///
/// Use [`with_mesh`] to create a mesh and token together:
///
/// ```rust,ignore
/// use cfd_mesh::mesh::with_mesh;
///
/// let volume = with_mesh(|mut mesh, mut token| {
///     // build mesh here, then extract a result
///     mesh.signed_volume(&token)
/// });
/// ```
///
/// # Invariants
///
/// For every half-edge `he` in a valid mesh:
/// 1. `twin(twin(he)) == he`  — twin is an involution
/// 2. `next(prev(he)) == he`  — next and prev are inverses
/// 3. `face(next(he)) == face(he)` — all half-edges in a face loop share a face
/// 4. Face loops and vertex rings terminate in finite steps
/// 5. All keys in `HalfEdgeData` reference live entries in their respective SlotMaps
///
/// # Diagram
///
/// ```text
/// HalfEdgeMesh<'id>
/// ├── GhostSlotPool<'id, VertexKey,   VertexData>
/// │     VertexData { position: Point3<Real>, half_edge: HalfEdgeKey }
/// ├── GhostSlotPool<'id, HalfEdgeKey, HalfEdgeData>
/// │     HalfEdgeData { vertex, face, twin, next, prev }
/// ├── GhostSlotPool<'id, FaceKey,     FaceData>
/// │     FaceData { half_edge, patch, normal }
/// └── SlotMap<PatchKey, BoundaryPatch>
///       BoundaryPatch { name, patch_type }
/// ```
pub struct HalfEdgeMesh<'id> {
    vertices:   GhostSlotPool<'id, VertexKey,   VertexData>,
    half_edges: GhostSlotPool<'id, HalfEdgeKey, HalfEdgeData>,
    faces:      GhostSlotPool<'id, FaceKey,     HeFaceData>,
    /// Patches do not need GhostCell: they are identified by `PatchKey` and
    /// are never aliased in topology, so `&mut self` is sufficient for mutation.
    patches:    SlotMap<PatchKey,    HeBoundaryPatch>,
    /// Twin lookup map: (origin_vertex, tip_vertex) → HalfEdgeKey.
    /// Used during `add_triangle` to find the twin of a new half-edge.
    twin_map:   hashbrown::HashMap<(VertexKey, VertexKey), HalfEdgeKey>,
}

impl<'id> HalfEdgeMesh<'id> {
    /// Create an empty mesh. Use [`with_mesh`] instead of calling this directly.
    fn new() -> Self {
        Self {
            vertices:   GhostSlotPool::new(),
            half_edges: GhostSlotPool::new(),
            faces:      GhostSlotPool::new(),
            patches:    SlotMap::with_key(),
            twin_map:   hashbrown::HashMap::new(),
        }
    }

    // ── Vertex operations ─────────────────────────────────────────────────

    /// Add a vertex at `position` and return its key.
    ///
    /// The `half_edge` field is initially set to a placeholder and must be
    /// updated when the first face using this vertex is added.
    ///
    /// # Example
    /// ```rust,ignore
    /// let vk = mesh.add_vertex([0.0, 0.0, 0.0], &token);
    /// ```
    pub fn add_vertex(
        &mut self,
        position: impl Into<Point3<Real>>,
        _token: &GhostToken<'id>,
    ) -> VertexKey {
        let data = VertexData::new(position.into());
        self.vertices.insert(GhostCell::new(data))
    }

    /// Number of vertices.
    #[inline]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Get the position of a vertex.
    pub fn vertex_pos(
        &self,
        key: VertexKey,
        token: &GhostToken<'id>,
    ) -> Option<Point3<Real>> {
        self.vertices.get(key).map(|c| c.borrow(token).position)
    }

    /// Mutate the position of a vertex.
    pub fn set_vertex_pos(
        &self,
        key: VertexKey,
        pos: Point3<Real>,
        token: &mut GhostToken<'id>,
    ) -> MeshResult<()> {
        let cell = self.vertices.get(key).ok_or(MeshError::Other("invalid vertex key".into()))?;
        cell.borrow_mut(token).position = pos;
        Ok(())
    }

    // ── Face operations ───────────────────────────────────────────────────

    /// Add a triangular face defined by three vertex keys.
    ///
    /// Creates three half-edges (`a→b`, `b→c`, `c→a`) and their boundary
    /// sentinel twins if no adjacent face exists yet. Sets `vertex.half_edge`
    /// for each vertex to an outgoing half-edge.
    ///
    /// # Errors
    /// Returns [`MeshError::InvalidVertexRef`] if any vertex key is stale.
    /// Returns [`MeshError::NonManifoldEdge`] if adding this triangle would
    /// create a non-manifold edge (more than 2 faces per edge).
    ///
    /// # Invariants maintained
    /// After a successful call, `twin(twin(he)) == he` for all new half-edges.
    pub fn add_triangle(
        &mut self,
        a: VertexKey,
        b: VertexKey,
        c: VertexKey,
        token: &mut GhostToken<'id>,
    ) -> MeshResult<FaceKey> {
        // Validate all keys
        for &vk in &[a, b, c] {
            if !self.vertices.contains_key(vk) {
                return Err(MeshError::Other("invalid vertex key".into()));
            }
        }

        // Check for non-manifold edges.
        // A directed edge (src→dst) as an INTERIOR half-edge (face.is_some()) may
        // appear at most once.  A SENTINEL entry (face.is_none()) is a boundary
        // placeholder waiting to be replaced — adding a face over it is legal.
        for &(src, dst) in &[(a, b), (b, c), (c, a)] {
            if let Some(&existing) = self.twin_map.get(&(src, dst)) {
                let is_interior = self.half_edges
                    .get(existing)
                    .map_or(false, |cell| cell.borrow(token).face.is_some());
                if is_interior {
                    return Err(MeshError::Other(
                        "non-manifold edge: directed edge already exists".into(),
                    ));
                }
            }
        }

        // Insert placeholder half-edges (we'll fill their fields next)
        let he_ab = self.half_edges.insert(GhostCell::new(HalfEdgeData {
            vertex: b,
            face: None, // filled below
            twin: HalfEdgeKey::default(),
            next: HalfEdgeKey::default(),
            prev: HalfEdgeKey::default(),
        }));
        let he_bc = self.half_edges.insert(GhostCell::new(HalfEdgeData {
            vertex: c,
            face: None,
            twin: HalfEdgeKey::default(),
            next: HalfEdgeKey::default(),
            prev: HalfEdgeKey::default(),
        }));
        let he_ca = self.half_edges.insert(GhostCell::new(HalfEdgeData {
            vertex: a,
            face: None,
            twin: HalfEdgeKey::default(),
            next: HalfEdgeKey::default(),
            prev: HalfEdgeKey::default(),
        }));

        // Compute face normal from vertex positions
        let pa = self.vertices[a].borrow(token).position;
        let pb = self.vertices[b].borrow(token).position;
        let pc = self.vertices[c].borrow(token).position;
        let edge1 = pb - pa;
        let edge2 = pc - pa;
        let cross = edge1.cross(&edge2);
        let normal = if cross.norm() > 1e-15 {
            nalgebra::UnitVector3::new_normalize(cross)
        } else {
            nalgebra::UnitVector3::new_unchecked(nalgebra::Vector3::z())
        };

        // Create the face
        let face_key = self.faces.insert(GhostCell::new(HeFaceData::new(he_ab, normal)));

        // Wire up half-edge fields: face, next, prev
        {
            let d_ab = self.half_edges[he_ab].borrow_mut(token);
            d_ab.face = Some(face_key);
            d_ab.next = he_bc;
            d_ab.prev = he_ca;
        }
        {
            let d_bc = self.half_edges[he_bc].borrow_mut(token);
            d_bc.face = Some(face_key);
            d_bc.next = he_ca;
            d_bc.prev = he_ab;
        }
        {
            let d_ca = self.half_edges[he_ca].borrow_mut(token);
            d_ca.face = Some(face_key);
            d_ca.next = he_ab;
            d_ca.prev = he_bc;
        }

        // Wire up twins — find or create boundary sentinels
        for (src, dst, he) in [(a, b, he_ab), (b, c, he_bc), (c, a, he_ca)] {
            let reverse = (dst, src);
            if let Some(&existing_twin) = self.twin_map.get(&reverse) {
                // The reverse half-edge already exists (from a previously added face)
                self.half_edges[he].borrow_mut(token).twin = existing_twin;
                self.half_edges[existing_twin].borrow_mut(token).twin = he;
            } else {
                // Create a boundary sentinel: a half-edge with face=None pointing
                // in the reverse direction, acting as the "exterior" twin.
                let sentinel = self.half_edges.insert(GhostCell::new(HalfEdgeData {
                    vertex: src,
                    face: None,
                    twin: he,
                    next: he, // placeholder; boundary loop is stitched separately
                    prev: he,
                }));
                self.half_edges[he].borrow_mut(token).twin = sentinel;
                // Register the sentinel so the next face can claim it as its twin
                self.twin_map.insert((dst, src), sentinel);
            }
            // Register the interior half-edge
            self.twin_map.insert((src, dst), he);
        }

        // Update vertex.half_edge to point to an outgoing half-edge
        self.vertices[a].borrow_mut(token).half_edge = he_ab;
        self.vertices[b].borrow_mut(token).half_edge = he_bc;
        self.vertices[c].borrow_mut(token).half_edge = he_ca;

        // Debug-mode invariant check
        #[cfg(debug_assertions)]
        self.check_triangle_invariants(face_key, token);

        Ok(face_key)
    }

    /// Number of faces.
    #[inline]
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Number of half-edges (includes boundary sentinels; each interior edge = 2).
    #[inline]
    pub fn half_edge_count(&self) -> usize {
        self.half_edges.len()
    }

    // ── Patch operations ──────────────────────────────────────────────────

    /// Create a named boundary patch.
    pub fn add_patch(&mut self, name: impl Into<String>, patch_type: PatchType) -> PatchKey {
        self.patches.insert(HeBoundaryPatch::new(name, patch_type))
    }

    /// Assign a face to a boundary patch.
    pub fn assign_face_to_patch(
        &mut self,
        face: FaceKey,
        patch: PatchKey,
        token: &mut GhostToken<'id>,
    ) -> MeshResult<()> {
        if !self.patches.contains_key(patch) {
            return Err(MeshError::Other("invalid patch key".into()));
        }
        let cell = self.faces.get(face).ok_or(MeshError::Other("invalid face key".into()))?;
        cell.borrow_mut(token).patch = Some(patch);
        Ok(())
    }

    /// Return the [`PatchKey`] assigned to a face, if any.
    ///
    /// Returns `None` for faces that have not been assigned to any boundary
    /// patch via [`assign_face_to_patch`](Self::assign_face_to_patch).
    pub fn face_patch(&self, face: FaceKey, token: &GhostToken<'id>) -> Option<PatchKey> {
        self.faces.get(face)?.borrow(token).patch
    }

    /// Return the [`BoundaryPatch`](crate::topology::halfedge::BoundaryPatch)
    /// registered under `patch`, or `None` if the key is stale.
    pub fn patch_info(&self, patch: PatchKey) -> Option<&HeBoundaryPatch> {
        self.patches.get(patch)
    }

    // ── Traversal iterators ───────────────────────────────────────────────

    /// Iterate over all vertex keys.
    pub fn vertex_keys<'a>(&'a self) -> impl Iterator<Item = VertexKey> + use<'a, 'id> {
        self.vertices.keys()
    }

    /// Iterate over all face keys (interior faces only — no boundary sentinels).
    pub fn face_keys<'a>(&'a self) -> impl Iterator<Item = FaceKey> + use<'a, 'id> {
        self.faces.keys()
    }

    /// Iterate over the half-edges bounding a face (the face loop).
    ///
    /// Yields half-edge keys in CCW order starting from `face.half_edge`.
    ///
    /// # Panics
    /// Panics if `face_key` is invalid or the face loop is malformed.
    pub fn face_half_edges<'a>(
        &'a self,
        face_key: FaceKey,
        token: &'a GhostToken<'id>,
    ) -> Vec<HalfEdgeKey> {
        let start = self.faces[face_key].borrow(token).half_edge;
        let mut result = Vec::with_capacity(3);
        let mut current = start;
        loop {
            result.push(current);
            current = self.half_edges[current].borrow(token).next;
            if current == start { break; }
            if result.len() > 65536 {
                panic!("face loop did not close — topology corrupted");
            }
        }
        result
    }

    /// Iterate over the vertices of a face in order.
    pub fn face_vertices<'a>(
        &'a self,
        face_key: FaceKey,
        token: &'a GhostToken<'id>,
    ) -> Vec<VertexKey> {
        self.face_half_edges(face_key, token)
            .iter()
            .map(|&he| self.half_edges[he].borrow(token).vertex)
            .collect()
    }

    // ── Half-edge primitive accessors ─────────────────────────────────────

    /// Return the twin of a half-edge.
    ///
    /// # Panics
    /// Panics if `he` is not a valid half-edge key.
    #[inline]
    pub fn he_twin(&self, he: HalfEdgeKey, token: &GhostToken<'id>) -> HalfEdgeKey {
        self.half_edges[he].borrow(token).twin
    }

    /// Return the face associated with a half-edge, or `None` for boundary sentinels.
    #[inline]
    pub fn he_face(&self, he: HalfEdgeKey, token: &GhostToken<'id>) -> Option<FaceKey> {
        self.half_edges[he].borrow(token).face
    }

    /// Return the next half-edge in the face loop.
    #[inline]
    pub fn he_next(&self, he: HalfEdgeKey, token: &GhostToken<'id>) -> HalfEdgeKey {
        self.half_edges[he].borrow(token).next
    }

    /// Return the previous half-edge in the face loop.
    #[inline]
    pub fn he_prev(&self, he: HalfEdgeKey, token: &GhostToken<'id>) -> HalfEdgeKey {
        self.half_edges[he].borrow(token).prev
    }

    /// Return the tip vertex of a half-edge (the vertex it points *to*).
    #[inline]
    pub fn he_vertex(&self, he: HalfEdgeKey, token: &GhostToken<'id>) -> VertexKey {
        self.half_edges[he].borrow(token).vertex
    }

    // ── Geometric queries ─────────────────────────────────────────────────

    /// Compute the signed volume of the mesh (positive for outward orientation).
    ///
    /// Uses the divergence theorem: `V = Σ_faces (v0 · (v1 × v2)) / 6`.
    pub fn signed_volume(&self, token: &GhostToken<'id>) -> Real {
        let mut vol = Real::default();
        for face_key in self.faces.keys() {
            let verts = self.face_vertices(face_key, token);
            if verts.len() == 3 {
                if let (Some(p0), Some(p1), Some(p2)) = (
                    self.vertex_pos(verts[0], token),
                    self.vertex_pos(verts[1], token),
                    self.vertex_pos(verts[2], token),
                ) {
                    vol += p0.coords.dot(&p1.coords.cross(&p2.coords));
                }
            }
        }
        vol / 6.0
    }

    // ── Debug invariant checks ────────────────────────────────────────────

    #[cfg(debug_assertions)]
    fn check_triangle_invariants(&self, face_key: FaceKey, token: &GhostToken<'id>) {
        let half_edges = self.face_half_edges(face_key, token);
        for &he in &half_edges {
            let twin = self.half_edges[he].borrow(token).twin;
            let twin_twin = self.half_edges[twin].borrow(token).twin;
            debug_assert_eq!(he, twin_twin, "twin(twin(he)) != he for {:?}", he);
            let next = self.half_edges[he].borrow(token).next;
            let prev_of_next = self.half_edges[next].borrow(token).prev;
            debug_assert_eq!(he, prev_of_next, "prev(next(he)) != he for {:?}", he);
        }
    }
}

// ── `with_mesh` entry point ────────────────────────────────────────────────

/// The canonical entry point for the GhostCell half-edge mesh.
///
/// Introduces a fresh brand `'id`, creates an empty [`HalfEdgeMesh<'id>`] and
/// a matching [`GhostToken<'id>`], and passes both to the closure `f`. The
/// brand cannot escape the closure, ensuring the mesh and token cannot be
/// mixed with an unrelated scope.
///
/// # Example
///
/// ```rust,ignore
/// use cfd_mesh::mesh::with_mesh;
///
/// let count = with_mesh(|mut mesh, mut token| {
///     let a = mesh.add_vertex([0.0, 0.0, 0.0], &token);
///     let b = mesh.add_vertex([1.0, 0.0, 0.0], &token);
///     let c = mesh.add_vertex([0.0, 1.0, 0.0], &token);
///     mesh.add_triangle(a, b, c, &mut token).unwrap();
///     mesh.face_count()
/// });
/// assert_eq!(count, 1);
/// ```
///
/// # Type Signature
///
/// The `for<'id>` bound in the closure signature is what prevents the brand
/// from escaping: the closure must work for *any* brand, so no specific `'id`
/// can be smuggled out.
pub fn with_mesh<F, R>(f: F) -> R
where
    F: for<'id> FnOnce(HalfEdgeMesh<'id>, GhostToken<'id>) -> R,
{
    GhostToken::new(|token| f(HalfEdgeMesh::new(), token))
}

#[cfg(test)]
mod he_tests {
    use super::*;

    #[test]
    fn add_single_triangle() {
        let face_count = with_mesh(|mut mesh, mut token| {
            let a = mesh.add_vertex(Point3::new(0.0, 0.0, 0.0), &token);
            let b = mesh.add_vertex(Point3::new(1.0, 0.0, 0.0), &token);
            let c = mesh.add_vertex(Point3::new(0.0, 1.0, 0.0), &token);
            mesh.add_triangle(a, b, c, &mut token).unwrap();
            (mesh.vertex_count(), mesh.face_count())
        });
        assert_eq!(face_count, (3, 1));
    }

    #[test]
    fn twin_involution_holds() {
        with_mesh(|mut mesh, mut token| {
            let a = mesh.add_vertex(Point3::new(0.0, 0.0, 0.0), &token);
            let b = mesh.add_vertex(Point3::new(1.0, 0.0, 0.0), &token);
            let c = mesh.add_vertex(Point3::new(0.0, 1.0, 0.0), &token);
            let face = mesh.add_triangle(a, b, c, &mut token).unwrap();

            for he in mesh.face_half_edges(face, &token) {
                let twin = mesh.half_edges[he].borrow(&token).twin;
                let twin_twin = mesh.half_edges[twin].borrow(&token).twin;
                assert_eq!(he, twin_twin, "twin(twin(he)) != he for {:?}", he);
            }
        });
    }

    #[test]
    fn two_adjacent_triangles_share_interior_twin() {
        with_mesh(|mut mesh, mut token| {
            // △ A-B-C and △ B-D-C sharing edge B-C
            let a = mesh.add_vertex(Point3::new(0.0, 0.0, 0.0), &token);
            let b = mesh.add_vertex(Point3::new(1.0, 0.0, 0.0), &token);
            let c = mesh.add_vertex(Point3::new(0.5, 1.0, 0.0), &token);
            let d = mesh.add_vertex(Point3::new(1.5, 1.0, 0.0), &token);
            let f1 = mesh.add_triangle(a, b, c, &mut token).unwrap();
            let f2 = mesh.add_triangle(b, d, c, &mut token).unwrap();

            // The shared edge B→C (in f1) and C→B (in f2) should be twins
            let hes1 = mesh.face_half_edges(f1, &token);
            let hes2 = mesh.face_half_edges(f2, &token);

            // Find the he in f1 pointing to c (the B→C half-edge)
            let he_bc = hes1.iter().copied().find(|&he| {
                mesh.half_edges[he].borrow(&token).vertex == c
            });
            let he_cb = hes2.iter().copied().find(|&he| {
                mesh.half_edges[he].borrow(&token).vertex == b
            });

            if let (Some(he_bc), Some(he_cb)) = (he_bc, he_cb) {
                let twin_of_bc = mesh.half_edges[he_bc].borrow(&token).twin;
                assert_eq!(twin_of_bc, he_cb, "twin of B→C should be C→B");
                let twin_of_cb = mesh.half_edges[he_cb].borrow(&token).twin;
                assert_eq!(twin_of_cb, he_bc, "twin of C→B should be B→C");
            }
        });
    }

    #[test]
    fn patch_assignment() {
        with_mesh(|mut mesh, mut token| {
            let a = mesh.add_vertex(Point3::new(0.0, 0.0, 0.0), &token);
            let b = mesh.add_vertex(Point3::new(1.0, 0.0, 0.0), &token);
            let c = mesh.add_vertex(Point3::new(0.0, 1.0, 0.0), &token);
            let face = mesh.add_triangle(a, b, c, &mut token).unwrap();
            let patch = mesh.add_patch("inlet", PatchType::Inlet);
            mesh.assign_face_to_patch(face, patch, &mut token).unwrap();

            let assigned = mesh.faces[face].borrow(&token).patch;
            assert_eq!(assigned, Some(patch));
        });
    }
}
