//! Core mesh data structure with topology and connectivity

use crate::topology::{Cell, Edge, ElementType, Face, Vertex};
use nalgebra::RealField;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Main mesh data structure
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct Mesh<T: RealField + Copy> {
    /// Vertices in the mesh
    vertices: Vec<Vertex<T>>,
    /// Edges connecting vertices
    edges: Vec<Edge>,
    /// Faces formed by edges
    faces: Vec<Face>,
    /// Cells formed by faces
    cells: Vec<Cell>,
    /// Boundary markers for faces
    boundary_markers: HashMap<usize, String>,
    /// Partition ID (rank) that owns this mesh (for distributed meshes)
    pub partition_id: Option<usize>,
}

impl<T: RealField + Copy> Mesh<T> {
    /// Create an empty mesh
    #[must_use]
    pub fn new() -> Self {
        Self {
            vertices: Vec::new(),
            edges: Vec::new(),
            faces: Vec::new(),
            cells: Vec::new(),
            boundary_markers: HashMap::new(),
            partition_id: None,
        }
    }

    /// Set the partition ID for this mesh
    pub fn set_partition_id(&mut self, partition_id: usize) {
        self.partition_id = Some(partition_id);
    }

    /// Add a vertex to the mesh
    pub fn add_vertex(&mut self, vertex: Vertex<T>) -> usize {
        let idx = self.vertices.len();
        self.vertices.push(vertex);
        idx
    }

    /// Add an edge to the mesh
    pub fn add_edge(&mut self, edge: Edge) -> usize {
        let idx = self.edges.len();
        self.edges.push(edge);
        idx
    }

    /// Add a face to the mesh
    pub fn add_face(&mut self, face: Face) -> usize {
        let idx = self.faces.len();
        self.faces.push(face);
        idx
    }

    /// Add a cell to the mesh
    pub fn add_cell(&mut self, cell: Cell) -> usize {
        let idx = self.cells.len();
        self.cells.push(cell);
        idx
    }

    /// Mark a face as boundary with a label
    pub fn mark_boundary(&mut self, face_idx: usize, label: String) {
        self.boundary_markers.insert(face_idx, label);
    }

    /// Get the bounding box of the mesh
    #[must_use]
    pub fn bounds(&self) -> (nalgebra::Point3<T>, nalgebra::Point3<T>) {
        if self.vertices.is_empty() {
            return (nalgebra::Point3::origin(), nalgebra::Point3::origin());
        }

        let mut min_x = self.vertices[0].position.x;
        let mut max_x = self.vertices[0].position.x;
        let mut min_y = self.vertices[0].position.y;
        let mut max_y = self.vertices[0].position.y;
        let mut min_z = self.vertices[0].position.z;
        let mut max_z = self.vertices[0].position.z;

        for v in &self.vertices[1..] {
            if v.position.x < min_x { min_x = v.position.x; }
            if v.position.x > max_x { max_x = v.position.x; }
            if v.position.y < min_y { min_y = v.position.y; }
            if v.position.y > max_y { max_y = v.position.y; }
            if v.position.z < min_z { min_z = v.position.z; }
            if v.position.z > max_z { max_z = v.position.z; }
        }

        (
            nalgebra::Point3::new(min_x, min_y, min_z),
            nalgebra::Point3::new(max_x, max_y, max_z),
        )
    }

    /// Get vertex by index
    #[must_use]
    pub fn vertex(&self, idx: usize) -> Option<&Vertex<T>> {
        self.vertices.get(idx)
    }

    /// Get edge by index
    #[must_use]
    pub fn edge(&self, idx: usize) -> Option<&Edge> {
        self.edges.get(idx)
    }

    /// Get face by index
    #[must_use]
    pub fn face(&self, idx: usize) -> Option<&Face> {
        self.faces.get(idx)
    }

    /// Get cell by index
    #[must_use]
    pub fn cell(&self, idx: usize) -> Option<&Cell> {
        self.cells.get(idx)
    }

    /// Number of vertices
    #[must_use]
    pub fn vertex_count(&self) -> usize {
        self.vertices.len()
    }

    /// Number of edges
    #[must_use]
    pub fn edge_count(&self) -> usize {
        self.edges.len()
    }

    /// Number of faces
    #[must_use]
    pub fn face_count(&self) -> usize {
        self.faces.len()
    }

    /// Number of cells
    #[must_use]
    pub fn cell_count(&self) -> usize {
        self.cells.len()
    }

    /// Get all topologically external face indices (referenced by exactly one cell)
    #[must_use]
    pub fn external_faces(&self) -> Vec<usize> {
        let mut face_cell_count = vec![0; self.faces.len()];
        for cell in &self.cells {
            for &f_idx in &cell.faces {
                if f_idx < face_cell_count.len() {
                    face_cell_count[f_idx] += 1;
                }
            }
        }
        
        face_cell_count
            .into_iter()
            .enumerate()
            .filter(|&(_, count)| count == 1)
            .map(|(idx, _)| idx)
            .collect()
    }

    /// Get all currently marked boundary face indices
    #[must_use]
    pub fn marked_boundary_faces(&self) -> Vec<usize> {
        self.boundary_markers.keys().copied().collect()
    }

    /// Backwards compatibility: currently marked or external (used by logic that expects to find all boundaries)
    #[must_use]
    pub fn boundary_faces(&self) -> Vec<usize> {
        let mut result: std::collections::HashSet<usize> = self.external_faces().into_iter().collect();
        for &idx in self.boundary_markers.keys() {
            result.insert(idx);
        }
        let mut v: Vec<usize> = result.into_iter().collect();
        v.sort_unstable();
        v
    }

    /// Get boundary label for a face
    pub fn boundary_label(&self, face_idx: usize) -> Option<&str> {
        self.boundary_markers.get(&face_idx).map(String::as_str)
    }

    /// Get all cells
    #[must_use]
    pub fn cells(&self) -> &[Cell] {
        &self.cells
    }

    /// Get all vertices  
    #[must_use]
    pub fn vertices(&self) -> &[Vertex<T>] {
        &self.vertices
    }

    /// Get mutable access to all vertices
    pub fn vertices_mut(&mut self) -> &mut [Vertex<T>] {
        &mut self.vertices
    }

    /// Get all edges
    #[must_use]
    pub fn edges(&self) -> &[Edge] {
        &self.edges
    }

    /// Get all faces
    #[must_use]
    pub fn faces(&self) -> &[Face] {
        &self.faces
    }

    /// Get faces of a cell as an iterator
    pub fn element_faces<'a>(&'a self, cell: &'a Cell) -> impl Iterator<Item = &'a Face> + 'a {
        cell.faces
            .iter()
            .filter_map(move |&face_idx| self.face(face_idx))
    }

    /// Get faces of a cell (allocating version for compatibility)
    #[deprecated(note = "Use element_faces() iterator for zero-copy access")]
    #[must_use]
    pub fn get_element_faces<'a>(&'a self, cell: &'a Cell) -> Vec<&'a Face> {
        self.element_faces(cell).collect()
    }

    /// Get vertices of a cell as an iterator (zero-copy)
    pub fn element_vertices<'a>(
        &'a self,
        cell: &'a Cell,
    ) -> impl Iterator<Item = &'a Vertex<T>> + 'a {
        use std::collections::HashSet;

        // Collect unique indices first (necessary for deduplication)
        let mut vertex_indices = HashSet::new();
        for &face_idx in &cell.faces {
            if let Some(face) = self.face(face_idx) {
                vertex_indices.extend(&face.vertices);
            }
        }

        // Return iterator over vertices
        vertex_indices
            .into_iter()
            .filter_map(move |idx| self.vertex(idx))
    }

    /// Get vertices of a cell (allocating version for compatibility)
    #[deprecated(note = "Use element_vertices() iterator for zero-copy access")]
    #[must_use]
    pub fn get_element_vertices<'a>(&'a self, cell: &'a Cell) -> Vec<&'a Vertex<T>> {
        self.element_vertices(cell).collect()
    }

    /// Get canonical ordered vertices for element quality metrics (sorted by vertex index)
    /// Invariant: sequential idx order matches standard element ordering (v0 apex/base[0], etc.)
    /// Lit: Required for Knupp Jacobian, Shewchuk linear elements
    pub fn ordered_element_vertices<'a>(&'a self, cell: &'a Cell) -> Vec<&'a Vertex<T>> {
        use std::collections::HashSet;
        let mut vset: HashSet<usize> = HashSet::new();
        for &face_idx in &cell.faces {
            if let Some(face) = self.face(face_idx) {
                vset.extend(face.vertices.iter().copied());
            }
        }
        let mut vidx: Vec<usize> = vset.into_iter().collect();
        vidx.sort_unstable();
        vidx.into_iter().map(|idx| self.vertex(idx).expect("Valid vertex index")).collect()
    }

    /// Merge another mesh into this one, welding vertices within `tolerance`.
    ///
    /// Vertices from `other` that are within `tolerance` of an existing vertex
    /// are fused; otherwise they are appended. Edges, faces, cells, and boundary
    /// markers are remapped accordingly.
    pub fn merge(&mut self, other: &Mesh<T>, tolerance: T) {
        let _base_vertex_count = self.vertices.len();
        let _base_edge_count = self.edges.len();
        let _base_face_count = self.faces.len();

        // Build a vertex index map: other_idx -> self_idx
        let mut vertex_map: Vec<usize> = Vec::with_capacity(other.vertices.len());
        for ov in &other.vertices {
            // Try to find an existing vertex within tolerance
            let mut matched = None;
            for (si, sv) in self.vertices.iter().enumerate() {
                if sv.distance_to(ov) < tolerance {
                    matched = Some(si);
                    break;
                }
            }
            match matched {
                Some(idx) => vertex_map.push(idx),
                None => {
                    let idx = self.vertices.len();
                    self.vertices.push(ov.clone());
                    vertex_map.push(idx);
                }
            }
        }

        // Remap and append edges
        for oe in &other.edges {
            let mut new_edge = oe.clone();
            new_edge.start = vertex_map[new_edge.start];
            new_edge.end = vertex_map[new_edge.end];
            self.edges.push(new_edge);
        }

        // Remap faces and deduplicate
        let mut face_map = Vec::with_capacity(other.faces.len());
        
        // Helper to get sorted vertices for comparison
        fn get_sorted_vertices<T: RealField + Copy>(face: &Face) -> Vec<usize> {
            let mut v = face.vertices.clone();
            v.sort_unstable();
            v
        }

        for of in &other.faces {
            let mut new_face = of.clone();
            new_face.vertices = new_face.vertices.iter().map(|&v| vertex_map[v]).collect();
            
            // Check for duplicate
            let new_sorted = get_sorted_vertices::<T>(&new_face);
            let mut found_idx = None;
            
            // Naive search (optimization: use a hashmap for large meshes)
            for (i, existing_face) in self.faces.iter().enumerate() {
                if existing_face.vertices.len() == new_face.vertices.len() {
                    let existing_sorted = get_sorted_vertices::<T>(existing_face);
                    if new_sorted == existing_sorted {
                        found_idx = Some(i);
                        break;
                    }
                }
            }
            
            match found_idx {
                Some(idx) => face_map.push(idx),
                None => {
                    let idx = self.faces.len();
                    self.faces.push(new_face);
                    face_map.push(idx);
                }
            }
        }

        // Remap and append cells
        for oc in &other.cells {
            let mut new_cell = oc.clone();
            new_cell.faces = new_cell.faces.iter().map(|&f| face_map[f]).collect();
            self.cells.push(new_cell);
        }

        // Remap and append boundary markers
        for (&face_idx, label) in &other.boundary_markers {
            // Only copy boundary marker if it's still a valid boundary face
            // (If it was merged, it might now be internal)
            // But usually we just want to preserve the label if the face exists
            // Since we remapped face_idx via face_map, we use that.
            self.boundary_markers.insert(face_map[face_idx], label.clone());
        }
    }

    /// Check mesh validity
    pub fn validate(&self) -> Result<(), String> {
        // Check vertex indices in edges
        for edge in &self.edges {
            if edge.start >= self.vertices.len() || edge.end >= self.vertices.len() {
                return Err(format!("Edge references invalid vertex: {edge:?}"));
            }
        }

        // Check vertex indices in faces
        for face in &self.faces {
            for &v in &face.vertices {
                if v >= self.vertices.len() {
                    return Err(format!("Face references invalid vertex: {v}"));
                }
            }
        }

        // Check face indices in cells
        for cell in &self.cells {
            for &f in &cell.faces {
                if f >= self.faces.len() {
                    return Err(format!("Cell references invalid face: {f}"));
                }
            }
        }

        Ok(())
    }

    /// Stitch another mesh onto this one at the specified boundary interface.
    ///
    /// This topologically fuses `other` to `self` by identifying vertices at the
    /// interface between `self_face_label` and `other_face_label`.
    pub fn stitch(
        &mut self, 
        other: &Mesh<T>, 
        self_label: &str, 
        other_label: &str,
        tolerance: T
    ) -> std::result::Result<(), String> {
        // 1. Identify interface vertices on both sides
        let self_interface_faces: Vec<usize> = self.boundary_markers.iter()
            .filter(|&(_, l)| l == self_label)
            .map(|(&f, _)| f)
            .collect();

        if self_interface_faces.is_empty() {
             return Err(format!("Self mesh has no faces with label '{self_label}'"));
        }

        let other_interface_faces: Vec<usize> = other.boundary_markers.iter()
            .filter(|&(_, l)| l == other_label)
            .map(|(&f, _)| f)
            .collect();
        
        if other_interface_faces.is_empty() {
            return Err(format!("Other mesh has no faces with label '{other_label}'"));
        }

        // Collect unique vertices at the interfaces
        let mut self_interface_verts = std::collections::HashSet::new();
        for &f_idx in &self_interface_faces {
            if let Some(face) = self.face(f_idx) {
                for &v in &face.vertices {
                    self_interface_verts.insert(v);
                }
            }
        }

        let mut other_interface_verts = std::collections::HashSet::new();
        for &f_idx in &other_interface_faces {
            if let Some(face) = other.face(f_idx) {
                for &v in &face.vertices {
                    other_interface_verts.insert(v);
                }
            }
        }
        
        if self_interface_verts.is_empty() {
             return Err("Found boundary faces but no vertices for self interface".to_string());
        }

        // 2. Build Vertex Mapping
        let mut vertex_map = vec![usize::MAX; other.vertices.len()];
        let mut mapped_count = 0;

        for &ov_idx in &other_interface_verts {
            let ov = &other.vertices[ov_idx];
            let mut found = None;
            let mut min_dist = tolerance; // Only look for vertices closer than tolerance
            
            for &sv_idx in &self_interface_verts {
                let sv = &self.vertices[sv_idx];
                let dist = sv.distance_to(ov);
                if dist < min_dist {
                    min_dist = dist;
                    found = Some(sv_idx);
                }
            }

            match found {
                Some(sv_idx) => {
                    vertex_map[ov_idx] = sv_idx;
                    mapped_count += 1;
                },
                None => {
                    return Err(format!(
                        "Failed to match interface vertex at {:?} (tolerance {:?})", 
                        ov.position, tolerance
                    ));
                }
            }
        }

        if mapped_count != other_interface_verts.len() {
             return Err(format!(
                 "Not all interface vertices mapped! mapped={} / total={}", 
                 mapped_count, other_interface_verts.len()
             ));
        }

        // For non-interface vertices, append to self
        for (i, v) in other.vertices.iter().enumerate() {
            if vertex_map[i] == usize::MAX {
                let new_idx = self.vertices.len();
                self.vertices.push(v.clone());
                vertex_map[i] = new_idx;
            }
        }

        // 3. Import Topology
        for e in &other.edges {
            let mut new_edge = e.clone();
            new_edge.start = vertex_map[e.start];
            new_edge.end = vertex_map[e.end];
            self.edges.push(new_edge);
        }

        // Build lookup of existing faces by canonical vertex set so stitched
        // interface faces are shared between adjacent cells (not duplicated).
        let mut canonical_face_map: std::collections::HashMap<Vec<usize>, usize> =
            std::collections::HashMap::new();
        for (i, face) in self.faces.iter().enumerate() {
            let mut key = face.vertices.clone();
            key.sort_unstable();
            canonical_face_map.insert(key, i);
        }

        let mut face_map = vec![usize::MAX; other.faces.len()];
        for (i, f) in other.faces.iter().enumerate() {
             let mut new_face = f.clone();
             new_face.vertices = f.vertices.iter().map(|&v| vertex_map[v]).collect();
             let mut key = new_face.vertices.clone();
             key.sort_unstable();

             if let Some(&existing_idx) = canonical_face_map.get(&key) {
                 face_map[i] = existing_idx;
             } else {
                 let idx = self.faces.len();
                 self.faces.push(new_face);
                 canonical_face_map.insert(key, idx);
                 face_map[i] = idx;
             }
        }

        for c in &other.cells {
            let mut new_cell = c.clone();
            new_cell.faces = c.faces.iter().map(|&f| face_map[f]).collect();
            self.cells.push(new_cell);
        }

        for (&f_idx, label) in &other.boundary_markers {
            let new_f_idx = face_map[f_idx];
            if label != other_label {
                self.boundary_markers.insert(new_f_idx, label.clone());
            }
        }
        
        for &f_idx in &self_interface_faces {
            self.boundary_markers.remove(&f_idx);
        }

        Ok(())
    }

    /// Compute mesh statistics
    #[must_use]
    #[allow(clippy::field_reassign_with_default)] // Clear, sequential initialization pattern
    pub fn statistics(&self) -> MeshStatistics {
        let mut stats = MeshStatistics::default();
        stats.vertex_count = self.vertices.len();
        stats.edge_count = self.edges.len();
        stats.face_count = self.faces.len();
        stats.cell_count = self.cells.len();
        stats.boundary_face_count = self.boundary_markers.len();

        // Count element types
        for cell in &self.cells {
            match cell.element_type {
                ElementType::Tetrahedron => stats.tetrahedra += 1,
                ElementType::Hexahedron => stats.hexahedra += 1,
                ElementType::Pyramid => stats.pyramids += 1,
                ElementType::Prism => stats.prisms += 1,
                _ => {}
            }
        }

        stats
    }
}

impl<T: RealField + Copy> Default for Mesh<T> {
    fn default() -> Self {
        Self::new()
    }
}

/// Mesh statistics
#[derive(Debug, Default, Clone)]
pub struct MeshStatistics {
    /// Number of vertices
    pub vertex_count: usize,
    /// Number of edges
    pub edge_count: usize,
    /// Number of faces
    pub face_count: usize,
    /// Number of cells
    pub cell_count: usize,
    /// Number of boundary faces
    pub boundary_face_count: usize,
    /// Number of tetrahedral cells
    pub tetrahedra: usize,
    /// Number of hexahedral cells
    pub hexahedra: usize,
    /// Number of pyramid cells
    pub pyramids: usize,
    /// Number of prism cells
    pub prisms: usize,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_mesh_operations() {
        let mut mesh = Mesh::<f64>::new();

        // Add vertices
        let v0 = mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex(Vertex::from_coords(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex(Vertex::from_coords(0.0, 1.0, 0.0));

        assert_eq!(mesh.vertex_count(), 3);

        // Add edge
        mesh.add_edge(Edge::new(v0, v1));
        assert_eq!(mesh.edge_count(), 1);

        // Add face
        mesh.add_face(Face::triangle(v0, v1, v2));
        assert_eq!(mesh.face_count(), 1);

        // Validate mesh
        assert!(mesh.validate().is_ok());
    }

    #[test]
    fn test_mesh_topology() {
        let mut mesh = Mesh::<f64>::new();

        // Create a simple tetrahedron
        mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(1.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(0.0, 1.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 1.0));

        // Add faces
        mesh.add_face(Face::triangle(0, 1, 2));
        mesh.add_face(Face::triangle(0, 1, 3));
        mesh.add_face(Face::triangle(0, 2, 3));
        mesh.add_face(Face::triangle(1, 2, 3));

        // Add cell
        mesh.add_cell(Cell::tetrahedron(0, 1, 2, 3));

        assert_eq!(mesh.vertex_count(), 4);
        assert_eq!(mesh.face_count(), 4);
        assert_eq!(mesh.cell_count(), 1);

        // Check statistics
        let stats = mesh.statistics();
        assert_eq!(stats.tetrahedra, 1);
        assert_eq!(stats.hexahedra, 0);
    }

    #[test]
    fn test_boundary_marking() {
        let mut mesh = Mesh::<f64>::new();

        // Add a simple face
        mesh.add_vertex(Vertex::from_coords(0.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(1.0, 0.0, 0.0));
        mesh.add_vertex(Vertex::from_coords(0.0, 1.0, 0.0));
        let face_idx = mesh.add_face(Face::triangle(0, 1, 2));

        // Mark as boundary
        mesh.mark_boundary(face_idx, "inlet".to_string());

        assert_eq!(mesh.boundary_faces().len(), 1);
        assert_eq!(mesh.boundary_label(face_idx), Some("inlet"));
    }

    #[test]
    fn test_mesh_stitching() {
        // Create two 1x1x1 cubes
        // Cube 1: [0,1] x [0,1] x [0,1]
        let mut mesh1 = crate::grid::StructuredGridBuilder::<f64>::new(1, 1, 1).build().unwrap();
        // Cube 2: [1,2] x [0,1] x [0,1] (shifted by 1 in x)
        let mut mesh2 = crate::grid::StructuredGridBuilder::<f64>::new(1, 1, 1)
            .with_bounds(((1.0, 2.0), (0.0, 1.0), (0.0, 1.0)))
            .build()
            .unwrap();

        // Mark interface faces
        // Mesh1: x=1 face (interface to mesh2)
        for f in 0..mesh1.face_count() {
            let face = mesh1.face(f).unwrap();
            let mut all_x1 = true;
            for &v in &face.vertices {
                if (mesh1.vertex(v).unwrap().position.x - 1.0).abs() > 1e-5 {
                    all_x1 = false;
                }
            }
            if all_x1 { mesh1.mark_boundary(f, "interface_1".to_string()); }
        }

        // Mesh2: x=1 face (interface to mesh1)
        for f in 0..mesh2.face_count() {
            let face = mesh2.face(f).unwrap();
            let mut all_x1 = true;
            for &v in &face.vertices {
                if (mesh2.vertex(v).unwrap().position.x - 1.0).abs() > 1e-5 {
                    all_x1 = false;
                }
            }
            if all_x1 { mesh2.mark_boundary(f, "interface_2".to_string()); }
        }

        // Stitch
        mesh1.stitch(&mesh2, "interface_1", "interface_2", 1e-4).unwrap();

        // Verify
        // 2 cubes = 2 cells
        assert_eq!(mesh1.cell_count(), 2);
        
        // Vertices:
        // Cube 1 has 8. Cube 2 has 8.
        // They share 4 vertices (the face at x=1).
        // Total should be 8 + 8 - 4 = 12.
        assert_eq!(mesh1.vertex_count(), 12, "Vertices should be stitched (merged)");
        
        // Ensure "interface_1" label is removed
        assert!(mesh1.boundary_markers.values().all(|l| l != "interface_1"));

        // The stitched interface should be internal (shared by two cells), not
        // counted as an external boundary.
        assert_eq!(mesh1.external_faces().len(), 10, "Joined cubes should have 10 external faces");
    }
}
