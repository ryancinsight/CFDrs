//! Problem definition for FEM

use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::domain::core::index::FaceId;
use cfd_mesh::IndexedMesh;
use nalgebra::{RealField, Vector3};
use std::collections::HashMap;

/// Problem definition for 3D incompressible flow using FEM
#[derive(Clone)]
pub struct StokesFlowProblem<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Computational mesh
    pub mesh: IndexedMesh<T>,
    /// Fluid properties
    pub fluid: ConstantPropertyFluid<T>,
    /// Boundary conditions mapped by node index
    pub boundary_conditions: HashMap<usize, BoundaryCondition<T>>,
    /// Body force (e.g., gravity)
    pub body_force: Option<Vector3<T>>,
    /// Spatially varying viscosity (per element)
    pub element_viscosities: Option<Vec<T>>,
    /// Number of corner nodes (original vertices before P2 conversion)
    /// Used for pressure degrees of freedom in Taylor-Hood
    pub n_corner_nodes: usize,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> StokesFlowProblem<T> {
    /// Create a new Stokes flow problem
    pub fn new(
        mesh: IndexedMesh<T>,
        fluid: ConstantPropertyFluid<T>,
        boundary_conditions: HashMap<usize, BoundaryCondition<T>>,
        n_corner_nodes: usize,
    ) -> Self {
        Self {
            mesh,
            fluid,
            boundary_conditions,
            body_force: None,
            element_viscosities: None,
            n_corner_nodes,
        }
    }

    /// Set body force (e.g., gravity)
    pub fn with_body_force(mut self, force: Vector3<T>) -> Self {
        self.body_force = Some(force);
        self
    }

    /// Validate problem setup
    pub fn validate(&self) -> Result<()> {
        // Check that all boundary nodes have boundary conditions
        let boundary_nodes = self.get_boundary_nodes();
        let missing_bcs: Vec<usize> = boundary_nodes
            .into_iter()
            .filter(|&node| !self.boundary_conditions.contains_key(&node))
            .collect();

        if !missing_bcs.is_empty() {
            return Err(Error::InvalidConfiguration(format!(
                "Missing boundary conditions for nodes: {missing_bcs:?}"
            )));
        }

        Ok(())
    }

    /// Get all boundary node indices
    ///
    /// Boundary nodes are vertices that belong to boundary faces.
    /// A boundary face is one that is referenced by only one cell,
    /// or is explicitly marked as a boundary face.
    ///
    /// # Algorithm
    ///
    /// 1. Identify boundary faces:
    ///    - Faces marked explicitly in boundary_markers
    ///    - Faces belonging to exactly one cell (external boundaries)
    /// 2. Collect unique vertices from all boundary faces
    ///
    /// # Returns
    ///
    /// Vector of unique node indices on the boundary
    ///
    /// # References
    ///
    /// - Zienkiewicz & Taylor (2000): The Finite Element Method, Vol 1
    /// - Hughes (2000): The Finite Element Method - boundary topology
    pub fn get_boundary_nodes(&self) -> Vec<usize> {
        use std::collections::{HashMap, HashSet};

        // 1. Geometric Deduplication Check
        // Count how many times each face (defined by sorted vertices) appears in cells
        // Key: Sorted list of vertex indices (representing a face)
        // Value: Count
        let mut face_geo_count: HashMap<Vec<usize>, usize> = HashMap::new();

        for cell in self.mesh.cells.iter() {
            for &face_idx in &cell.faces {
                if face_idx < self.mesh.face_count() {
                    let face = self.mesh.faces.get(FaceId::from_usize(face_idx));
                    let mut verts: Vec<usize> = face.vertices.iter().map(|v| v.as_usize()).collect();
                    verts.sort_unstable();
                    *face_geo_count.entry(verts).or_insert(0) += 1;
                }
            }
        }

        let mut boundary_vertices: HashSet<usize> = HashSet::new();

        // Add vertices from geometrically unique faces (count == 1)
        for (verts, count) in &face_geo_count {
            if *count == 1 {
                for &v_idx in verts {
                    boundary_vertices.insert(v_idx);
                }
            }
        }

        // 2. Add vertices from explicitly marked boundary faces
        // We use boundary_labels directly to avoid relying on topological boundary_faces()
        for (face_id, _) in &self.mesh.boundary_labels {
             let face = self.mesh.faces.get(*face_id);
             for &v_id in &face.vertices {
                 boundary_vertices.insert(v_id.as_usize());
             }
        }

        // Convert to sorted vector for deterministic output
        let mut result: Vec<usize> = boundary_vertices.into_iter().collect();
        result.sort_unstable();
        result
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_mesh::IndexedMesh;
    use cfd_mesh::domain::core::index::VertexId;
    use cfd_mesh::domain::topology::{Cell, Face};
    use nalgebra::Point3;

    /// Create a simple tetrahedral mesh for testing
    ///
    /// Mesh structure:
    /// ```text
    ///       v3 (top)
    ///       /|\
    ///      / | \
    ///     /  |  \
    ///   v0---+---v2
    ///     \  |  /
    ///      \ | /
    ///       \|/
    ///       v1 (bottom)
    /// ```
    ///
    /// One tetrahedron with 4 triangular faces, all are boundary faces
    fn create_test_tet_mesh() -> IndexedMesh<f64> {
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3::new(0.5, 1.0, 0.0));
        let v3 = mesh.add_vertex_pos(Point3::new(0.5, 0.5, 1.0));

        let f0 = mesh.add_face(v0, v1, v2).0;
        let f1 = mesh.add_face(v0, v1, v3).0;
        let f2 = mesh.add_face(v1, v2, v3).0;
        let f3 = mesh.add_face(v2, v0, v3).0;

        mesh.add_cell(Cell::tetrahedron(f0 as usize, f1 as usize, f2 as usize, f3 as usize));
        mesh
    }

    /// Create a two-tet mesh to test internal face detection
    ///
    /// Mesh structure: two tetrahedra sharing a face
    /// ```text
    ///    v3     v4
    ///    /|\   /|
    ///   / | \ / |
    ///  v0-+-v2--+
    ///   \ | / \ |
    ///    \|/   \|
    ///    v1     v1
    /// ```
    fn create_two_tet_mesh() -> IndexedMesh<f64> {
        let mut mesh = IndexedMesh::new();
        let v0 = mesh.add_vertex_pos(Point3::new(0.0, 0.0, 0.0));
        let v1 = mesh.add_vertex_pos(Point3::new(1.0, 0.0, 0.0));
        let v2 = mesh.add_vertex_pos(Point3::new(0.5, 1.0, 0.0));
        let v3 = mesh.add_vertex_pos(Point3::new(0.5, 0.5, 1.0));
        let v4 = mesh.add_vertex_pos(Point3::new(1.5, 0.5, 1.0));

        let f0 = mesh.add_face(v0, v1, v2).0;
        let f1 = mesh.add_face(v0, v1, v3).0;
        let f2 = mesh.add_face(v1, v2, v3).0; // shared
        let f3 = mesh.add_face(v2, v0, v3).0;

        let f4 = mesh.add_face(v1, v2, v4).0;
        let f5 = mesh.add_face(v1, v3, v4).0;
        let f6 = mesh.add_face(v2, v3, v4).0;

        mesh.add_cell(Cell::tetrahedron(f0 as usize, f1 as usize, f2 as usize, f3 as usize));
        mesh.add_cell(Cell::tetrahedron(f2 as usize, f4 as usize, f5 as usize, f6 as usize));
        mesh
    }

    #[test]
    fn test_boundary_detection_single_tet() {
        let mesh = create_test_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();
        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);
        let boundary_nodes = problem.get_boundary_nodes();
        assert_eq!(boundary_nodes.len(), 4);
    }

    #[test]
    fn test_boundary_detection_two_tets() {
        let mesh = create_two_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();
        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 5);
        let boundary_nodes = problem.get_boundary_nodes();
        assert_eq!(boundary_nodes.len(), 5);
    }

    #[test]
    fn test_boundary_detection_with_markers() {
        let mut mesh = create_test_tet_mesh();
        use cfd_mesh::domain::core::index::FaceId;
        mesh.mark_boundary(FaceId::from_usize(0), "inlet".to_string());
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();
        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);
        let boundary_nodes = problem.get_boundary_nodes();
        assert_eq!(boundary_nodes.len(), 4);
    }

    #[test]
    fn test_validate_with_missing_boundary_conditions() {
        let mesh = create_test_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();
        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);
        assert!(problem.validate().is_err());
    }

    #[test]
    fn test_validate_with_complete_boundary_conditions() {
        let mesh = create_test_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let mut boundary_conditions = HashMap::new();
        for i in 0..4 {
            boundary_conditions.insert(i, BoundaryCondition::Dirichlet { value: 0.0, component_values: None });
        }
        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);
        assert!(problem.validate().is_ok());
    }

    #[test]
    fn test_boundary_nodes_sorted() {
        let mesh = create_test_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();
        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);
        let boundary_nodes = problem.get_boundary_nodes();
        let mut sorted = boundary_nodes.clone();
        sorted.sort_unstable();
        assert_eq!(boundary_nodes, sorted);
    }

    /// Create a star mesh (12 tets) around a center node, with duplicated faces
    fn create_bad_star_mesh() -> IndexedMesh<f64> {
        let mut mesh = IndexedMesh::new();

        // Vertices
        let c = mesh.add_vertex_pos(Point3::new(0.5, 0.5, 0.5)); // Center (index 0)
        let v000 = mesh.add_vertex_pos(Point3::new(0.0, 0.0, 0.0)); // 1
        let v100 = mesh.add_vertex_pos(Point3::new(1.0, 0.0, 0.0)); // 2
        let v110 = mesh.add_vertex_pos(Point3::new(1.0, 1.0, 0.0)); // 3
        let v010 = mesh.add_vertex_pos(Point3::new(0.0, 1.0, 0.0)); // 4
        let v001 = mesh.add_vertex_pos(Point3::new(0.0, 0.0, 1.0)); // 5
        let v101 = mesh.add_vertex_pos(Point3::new(1.0, 0.0, 1.0)); // 6
        let v111 = mesh.add_vertex_pos(Point3::new(1.0, 1.0, 1.0)); // 7
        let v011 = mesh.add_vertex_pos(Point3::new(0.0, 1.0, 1.0)); // 8

        let corners = [v000, v100, v110, v010, v001, v101, v111, v011];
        // Indices relative to corners array:
        // Bottom: 0,1,2,3 -> (0,0,0), (1,0,0), (1,1,0), (0,1,0).
        // Top: 4,5,6,7 -> (0,0,1), (1,0,1), (1,1,1), (0,1,1).

        let surface_tris = vec![
            // Bottom (z=0)
            (corners[0], corners[3], corners[2]), (corners[0], corners[2], corners[1]),
            // Top (z=1)
            (corners[4], corners[5], corners[6]), (corners[4], corners[6], corners[7]),
            // Front (y=0)
            (corners[0], corners[1], corners[5]), (corners[1], corners[5], corners[6]), // No, 1-2-6-5 logic?
            // Let's use simpler explicit coords to be sure.
            // v000, v100, v101, v001. -> 0,1,5,4.
            (corners[0], corners[1], corners[5]), (corners[0], corners[5], corners[4]),
            // Back (y=1)
            // v010, v110, v111, v011. -> 3,2,6,7.
            (corners[3], corners[2], corners[6]), (corners[3], corners[6], corners[7]),
            // Left (x=0)
            // v000, v001, v011, v010. -> 0,4,7,3.
            (corners[0], corners[4], corners[7]), (corners[0], corners[7], corners[3]),
            // Right (x=1)
            // v100, v110, v111, v101. -> 1,2,6,5.
            (corners[1], corners[2], corners[6]), (corners[1], corners[6], corners[5]),
        ];

        for (v1, v2, v3) in surface_tris {
            // Create tet connecting Center to this triangle
            // Intentionally create NEW faces for each tet
            let f0 = mesh.add_face(c, v1, v2);
            let f1 = mesh.add_face(c, v2, v3);
            let f2 = mesh.add_face(c, v3, v1);
            let f3 = mesh.add_face(v1, v2, v3); // Surface face

            mesh.add_cell(Cell::tetrahedron(f0.as_usize(), f1.as_usize(), f2.as_usize(), f3.as_usize()));
        }

        mesh
    }

    #[test]
    fn test_internal_node_is_not_boundary_with_bad_topology() {
        // Create a mesh with duplicated faces
        let mesh = create_bad_star_mesh();

        // Node 0 is the center node
        let center_node = 0;

        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 9);
        let boundary_nodes = problem.get_boundary_nodes();

        // With geometric deduplication, the internal node (0) should NOT be in boundary nodes.
        // It is shared by 12 tets. The faces connected to it are internal.

        assert!(!boundary_nodes.contains(&center_node),
            "Internal node {} should NOT be identified as boundary node, even with duplicated faces. Boundary nodes found: {:?}", center_node, boundary_nodes);
    }
}
