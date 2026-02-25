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

    /// Repair any inconsistent face windings so that all normals point outward.
    ///
    /// Uses a manifold BFS flood from the extremal face to determine globally
    /// consistent outward orientation, then flips any inward-facing face's
    /// winding in-place (`v1 ↔ v2`). Finally calls [`recompute_normals`] to
    /// synchronise vertex normals with the repaired geometry.
    ///
    /// ## Theorem basis
    ///
    /// For any closed orientable 2-manifold M embedded in ℝ³, the face with
    /// the vertex carrying the highest X coordinate must have an outward normal
    /// with `n_x ≥ 0` (Jordan-Brouwer separation theorem applied to the +X
    /// axis half-space). BFS propagation from this seed via the half-edge
    /// adjacency graph assigns a globally consistent orientation label to every
    /// reachable face. Flipping only the faces labelled *inward* corrects the
    /// minority misclassified by Phase-4 GWN seam ambiguity in the CSG
    /// pipeline (Turk & Levoy 1994; Zhou et al. 2016 "Mesh Arrangements for
    /// Solid Geometry") without disturbing the majority that are already correct.
    ///
    /// ## Properties
    ///
    /// - **Topology-preserving**: only vertex ordering within each face changes;
    ///   no vertices or edges are created, deleted, or repositioned.
    /// - **Watertightness-preserving**: the undirected edge graph is unchanged.
    /// - **O(F + E)** time and space (same as the half-edge BFS).
    /// - **No-op** on meshes already all-outward; safe to call unconditionally.
    ///
    /// ## Disconnected components
    ///
    /// Each connected component is re-seeded from its own extremal face, so
    /// multi-component meshes (e.g. a chip body with separate channel voids)
    /// are handled correctly.
    pub fn orient_outward(&mut self) {
        use std::collections::{HashMap, VecDeque};
        use crate::domain::geometry::normal::triangle_normal;

        // Collect an owned copy of all face data so the immutable borrow ends
        // before the mutable `self.faces.iter_mut()` pass below.
        use crate::infrastructure::storage::face_store::FaceData;
        let face_list: Vec<FaceData> = self.faces.iter().copied().collect();
        let n_faces: usize = face_list.len();
        if n_faces == 0 {
            return;
        }

        // Per-face normals (None = degenerate) and max-X vertex position.
        // Both are computed from the immutable vertex pool before the later
        // mutable pass.
        let mut face_normals: Vec<Option<Vector3<T>>> = Vec::with_capacity(n_faces);
        let mut face_max_x: Vec<T> = Vec::with_capacity(n_faces);
        for face in &face_list {
            let a = self.vertices.position(face.vertices[0]);
            let b = self.vertices.position(face.vertices[1]);
            let c = self.vertices.position(face.vertices[2]);
            face_normals.push(triangle_normal(a, b, c));
            // Use explicit comparison to avoid ambiguity between
            // nalgebra::RealField::max and num_traits::Float::max.
            let ab = if a.x > b.x { a.x } else { b.x };
            face_max_x.push(if ab > c.x { ab } else { c.x });
        }

        // Directed half-edge → face index.
        let half_edge_cap: usize = n_faces * 3;
        let mut half_edge: HashMap<(VertexId, VertexId), usize> =
            HashMap::with_capacity(half_edge_cap);
        for (fi, face) in face_list.iter().enumerate() {
            let v = face.vertices;
            for k in 0..3 {
                let j = (k + 1) % 3;
                half_edge.insert((v[k], v[j]), fi);
            }
        }

        // BFS orientation labels: Some(true) = outward, Some(false) = inward.
        let mut orientation: Vec<Option<bool>> = vec![None; n_faces];

        // Outer loop handles disconnected components — each gets its own seed.
        loop {
            // Find the unvisited non-degenerate face with the maximum X vertex.
            let seed_fi = {
                let mut best_x = <T as num_traits::Float>::neg_infinity();
                let mut best: Option<usize> = None;
                for fi in 0..n_faces {
                    if orientation[fi].is_some() || face_normals[fi].is_none() {
                        continue;
                    }
                    if face_max_x[fi] > best_x {
                        best_x = face_max_x[fi];
                        best = Some(fi);
                    }
                }
                match best {
                    Some(fi) => fi,
                    None => break,
                }
            };

            // Seed orientation: the outward normal of the extremal face must
            // have a non-negative X component.
            let seed_normal = face_normals[seed_fi].unwrap();
            orientation[seed_fi] = Some(seed_normal.x >= T::zero());

            let mut queue: VecDeque<usize> = VecDeque::new();
            queue.push_back(seed_fi);

            while let Some(fi) = queue.pop_front() {
                let is_outward = orientation[fi].unwrap();
                let v = face_list[fi].vertices;
                for k in 0..3 {
                    let j = (k + 1) % 3;
                    let va = v[k];
                    let vb = v[j];
                    // Reverse edge (vb→va): consistent adjacency → same orientation.
                    if let Some(&nfi) = half_edge.get(&(vb, va)) {
                        if orientation[nfi].is_none() && face_normals[nfi].is_some() {
                            orientation[nfi] = Some(is_outward);
                            queue.push_back(nfi);
                        }
                    // Forward edge (va→vb): same direction → winding flip.
                    } else if let Some(&nfi) = half_edge.get(&(va, vb)) {
                        if orientation[nfi].is_none() && face_normals[nfi].is_some() {
                            orientation[nfi] = Some(!is_outward);
                            queue.push_back(nfi);
                        }
                    }
                }
            }
        }

        // Flip inward faces in-place (swap v1 ↔ v2).
        let mut fi = 0usize;
        for face in self.faces.iter_mut() {
            if orientation[fi] == Some(false) {
                face.flip();
            }
            fi += 1;
        }
        self.edges = None;

        // Synchronise vertex normals with the repaired winding.
        self.recompute_normals();
    }

    /// Remove all connected face components except the largest one.
    ///
    /// After a CSG Difference operation, phantom closed "islands" — small groups
    /// of faces that are locally manifold but disconnected from the main body —
    /// can appear at the seam boundary.  Each island passes `is_watertight = true`
    /// individually but inflates the Euler characteristic (χ = 4 instead of 2),
    /// corrupts the signed volume, and produces floating geometry in STL output.
    ///
    /// **Threshold:** a component is discarded when its face count is less than
    /// `max(4, largest_component_size × 5%)`.  The 5 % relative threshold with
    /// a 4-face absolute minimum reliably suppresses seam artefacts while
    /// preserving large intentional secondary bodies.
    ///
    /// **Field handling after reconstruction:**
    ///
    /// | Field             | Action                                      |
    /// |-------------------|---------------------------------------------|
    /// | `vertices`        | Rebuilt — orphaned vertices removed          |
    /// | `faces`           | Rebuilt from kept-component faces only       |
    /// | `edges`           | Set to `None` — lazily rebuilt on next call  |
    /// | `cells`           | Cleared — CSG surfaces carry no cell data    |
    /// | `attributes`      | Remapped via old→new `FaceId` translation    |
    /// | `boundary_labels` | Remapped; labels on discarded faces dropped  |
    ///
    /// # Returns
    /// Number of components discarded (`0` when the mesh was already clean).
    pub fn retain_largest_component(&mut self) -> usize {
        use crate::domain::topology::AdjacencyGraph;
        use crate::domain::topology::connectivity::connected_components;

        // Ensure edge adjacency is current.
        self.rebuild_edges();
        let edges = self.edges.as_ref().expect("edges just rebuilt");

        let adj = AdjacencyGraph::build(&self.faces, edges);
        let components = connected_components(&self.faces, &adj);

        // Fast path: already a single component.
        if components.len() <= 1 {
            return 0;
        }

        let largest_size = components.iter().map(|c| c.len()).max().unwrap_or(0);
        // Discard if face_count < max(4, largest * 0.05).
        let min_keep = ((largest_size as f64 * 0.05).ceil() as usize).max(4);

        // Single-pass over components: build a fresh mesh from kept faces only.
        let mut new_mesh = IndexedMesh::<T>::new();
        // Map old VertexId → new VertexId (None = not yet seen).
        let mut vertex_remap: Vec<Option<VertexId>> = vec![None; self.vertices.len()];
        // Map old FaceId → new FaceId for attribute/label remapping.
        let mut face_remap: HashMap<FaceId, FaceId> = HashMap::new();
        let mut discarded = 0usize;

        for component in &components {
            if component.len() < min_keep {
                discarded += 1;
                tracing::debug!(
                    "retain_largest_component: discarding {} phantom face(s) \
                     (threshold = {} faces)",
                    component.len(),
                    min_keep,
                );
                continue;
            }
            for &old_fid in component {
                let fd = *self.faces.get(old_fid);
                let mut nv = [VertexId::default(); 3];
                for (k, &vid) in fd.vertices.iter().enumerate() {
                    let idx = vid.as_usize();
                    nv[k] = *vertex_remap[idx].get_or_insert_with(|| {
                        new_mesh.add_vertex(
                            *self.vertices.position(vid),
                            *self.vertices.normal(vid),
                        )
                    });
                }
                // Guard: skip any face that collapsed under vertex welding.
                if nv[0] == nv[1] || nv[1] == nv[2] || nv[2] == nv[0] {
                    continue;
                }
                let new_fid = if fd.region == RegionId::INVALID {
                    new_mesh.add_face(nv[0], nv[1], nv[2])
                } else {
                    new_mesh.add_face_with_region(nv[0], nv[1], nv[2], fd.region)
                };
                face_remap.insert(old_fid, new_fid);
            }
        }

        if discarded == 0 {
            return 0;
        }

        // Remap per-face scalar attributes.
        let old_attrs = std::mem::replace(&mut self.attributes, AttributeStore::new());
        for channel in old_attrs.channel_names() {
            for (&old_fid, &new_fid) in &face_remap {
                if let Some(val) = old_attrs.get(channel, old_fid) {
                    new_mesh.attributes.set(channel, new_fid, val);
                }
            }
        }

        // Remap boundary labels.
        let old_labels = std::mem::replace(&mut self.boundary_labels, HashMap::new());
        new_mesh.boundary_labels = old_labels
            .into_iter()
            .filter_map(|(old_fid, label)| {
                face_remap.get(&old_fid).map(|&nf| (nf, label))
            })
            .collect();

        // Swap stores in-place.
        self.vertices        = new_mesh.vertices;
        self.faces           = new_mesh.faces;
        self.edges           = None;         // stale; lazily rebuilt on next use
        self.cells           = Vec::new();   // CSG surfaces carry no volumetric cells
        self.attributes      = new_mesh.attributes;
        self.boundary_labels = new_mesh.boundary_labels;

        tracing::debug!(
            "retain_largest_component: removed {} component(s); {} faces remain",
            discarded,
            self.faces.len(),
        );
        discarded
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

