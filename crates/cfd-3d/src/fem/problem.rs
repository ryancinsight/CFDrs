//! Problem definition for FEM

use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::physics::fluid::ConstantPropertyFluid;
use cfd_mesh::mesh::Mesh;
use nalgebra::{RealField, Vector3};
use std::collections::HashMap;

/// Problem definition for 3D incompressible flow using FEM
#[derive(Debug, Clone)]
pub struct StokesFlowProblem<T: RealField + Copy> {
    /// Computational mesh
    pub mesh: Mesh<T>,
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

impl<T: RealField + Copy> StokesFlowProblem<T> {
    /// Create a new Stokes flow problem
    pub fn new(
        mesh: Mesh<T>,
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
    fn get_boundary_nodes(&self) -> Vec<usize> {
        use std::collections::HashSet;

        // Only consider topologically external faces (referenced by exactly one cell).
        // This avoids flagging internal nodes as boundary nodes even if an internal face is marked.
        let external_faces = self.mesh.external_faces();

        // Collect unique vertices from all external faces
        let mut boundary_vertices: HashSet<usize> = HashSet::new();
        for &face_idx in &external_faces {
            if let Some(face) = self.mesh.face(face_idx) {
                boundary_vertices.extend(&face.vertices);
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
    use cfd_mesh::topology::{Cell, Face};
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
    fn create_test_tet_mesh() -> Mesh<f64> {
        let mut mesh = Mesh::new();

        // Add 4 vertices
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(0.0, 0.0, 0.0)));
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(1.0, 0.0, 0.0)));
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(0.5, 1.0, 0.0)));
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(0.5, 0.5, 1.0)));

        // Add 4 triangular faces
        let f0 = mesh.add_face(Face::triangle(0, 1, 2)); // bottom
        let f1 = mesh.add_face(Face::triangle(0, 1, 3)); // front
        let f2 = mesh.add_face(Face::triangle(1, 2, 3)); // right
        let f3 = mesh.add_face(Face::triangle(2, 0, 3)); // left

        // Add 1 tetrahedral cell
        mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));

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
    fn create_two_tet_mesh() -> Mesh<f64> {
        let mut mesh = Mesh::new();

        // Add 5 vertices
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(0.0, 0.0, 0.0)));
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(1.0, 0.0, 0.0)));
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(0.5, 1.0, 0.0)));
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(0.5, 0.5, 1.0)));
        mesh.add_vertex(cfd_mesh::topology::Vertex::new(Point3::new(1.5, 0.5, 1.0)));

        // First tetrahedron faces
        let f0 = mesh.add_face(Face::triangle(0, 1, 2)); // bottom tet1
        let f1 = mesh.add_face(Face::triangle(0, 1, 3)); // front tet1
        let f2 = mesh.add_face(Face::triangle(1, 2, 3)); // shared face (internal)
        let f3 = mesh.add_face(Face::triangle(2, 0, 3)); // left tet1

        // Second tetrahedron faces (shares f2)
        let f4 = mesh.add_face(Face::triangle(1, 2, 4)); // right tet2
        let f5 = mesh.add_face(Face::triangle(1, 3, 4)); // bottom tet2
        let f6 = mesh.add_face(Face::triangle(2, 3, 4)); // top tet2

        // Add cells
        mesh.add_cell(Cell::tetrahedron(f0, f1, f2, f3));
        mesh.add_cell(Cell::tetrahedron(f2, f4, f5, f6));

        mesh
    }

    #[test]
    fn test_boundary_detection_single_tet() {
        // Single tetrahedron - all 4 vertices should be on boundary
        let mesh = create_test_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);
        let boundary_nodes = problem.get_boundary_nodes();

        // All 4 vertices should be boundary nodes
        assert_eq!(boundary_nodes.len(), 4);
        assert_eq!(boundary_nodes, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_boundary_detection_two_tets() {
        // Two tetrahedra sharing a face - shared vertices should still be boundary
        // because they are on external faces
        let mesh = create_two_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 5);
        let boundary_nodes = problem.get_boundary_nodes();

        // All 5 vertices are on boundary (even shared ones appear on external faces)
        assert_eq!(boundary_nodes.len(), 5);
        assert!(boundary_nodes.contains(&0));
        assert!(boundary_nodes.contains(&1));
        assert!(boundary_nodes.contains(&2));
        assert!(boundary_nodes.contains(&3));
        assert!(boundary_nodes.contains(&4));
    }

    #[test]
    fn test_boundary_detection_with_markers() {
        // Test that explicitly marked boundary faces are detected
        let mut mesh = create_test_tet_mesh();

        // Mark face 0 as "inlet"
        mesh.mark_boundary(0, "inlet".to_string());

        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);
        let boundary_nodes = problem.get_boundary_nodes();

        // All 4 vertices should still be boundary nodes
        assert_eq!(boundary_nodes.len(), 4);
        assert_eq!(boundary_nodes, vec![0, 1, 2, 3]);
    }

    #[test]
    fn test_validate_with_missing_boundary_conditions() {
        // Test validation catches missing boundary conditions
        let mesh = create_test_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new(); // Empty - no BCs defined

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);

        // Should fail validation - boundary nodes exist but no BCs
        let result = problem.validate();
        assert!(result.is_err());

        if let Err(Error::InvalidConfiguration(msg)) = result {
            assert!(msg.contains("Missing boundary conditions"));
            assert!(msg.contains("nodes:"));
        }
    }

    #[test]
    fn test_validate_with_complete_boundary_conditions() {
        // Test validation passes when all boundary nodes have BCs
        let mesh = create_test_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();

        // Add BCs for all 4 boundary nodes (simple Dirichlet)
        let mut boundary_conditions = HashMap::new();
        for i in 0..4 {
            boundary_conditions.insert(
                i,
                BoundaryCondition::Dirichlet {
                    value: 0.0,
                    component_values: None,
                },
            );
        }

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);

        // Should pass validation
        let result = problem.validate();
        assert!(result.is_ok());
    }

    #[test]
    fn test_boundary_nodes_sorted() {
        // Test that boundary nodes are returned in sorted order
        let mesh = create_test_tet_mesh();
        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 4);
        let boundary_nodes = problem.get_boundary_nodes();

        // Check sorted order
        let mut sorted = boundary_nodes.clone();
        sorted.sort_unstable();
        assert_eq!(boundary_nodes, sorted);
    }
    #[test]
    fn test_interior_node_ignored() {
        use cfd_mesh::topology::{Cell, Face, Vertex};
        use nalgebra::Point3;

        // Create a 4-tet mesh by subdividing a large tetrahedron with a central node.
        // Vertices:
        // 0: (0,0,0)
        // 1: (1,0,0)
        // 2: (0,1,0)
        // 3: (0,0,1)
        // 4: (0.2, 0.2, 0.2) - Interior node

        let mut mesh = Mesh::new();
        mesh.add_vertex(Vertex::new(Point3::new(0.0, 0.0, 0.0)));
        mesh.add_vertex(Vertex::new(Point3::new(1.0, 0.0, 0.0)));
        mesh.add_vertex(Vertex::new(Point3::new(0.0, 1.0, 0.0)));
        mesh.add_vertex(Vertex::new(Point3::new(0.0, 0.0, 1.0)));
        mesh.add_vertex(Vertex::new(Point3::new(0.2, 0.2, 0.2)));

        // Tet 1: 0,1,2,4 (Base 0,1,2, Peak 4)
        // Faces: (0,1,2), (0,1,4), (1,2,4), (2,0,4)
        let f0 = mesh.add_face(Face::triangle(0, 1, 2)); // External
        let f_014 = mesh.add_face(Face::triangle(0, 1, 4)); // Internal
        let f_124 = mesh.add_face(Face::triangle(1, 2, 4)); // Internal
        let f_204 = mesh.add_face(Face::triangle(2, 0, 4)); // Internal
        mesh.add_cell(Cell::tetrahedron(f0, f_014, f_124, f_204));

        // Tet 2: 0,1,3,4 (Base 0,1,3, Peak 4)
        // Faces: (0,1,3), (0,1,4)*, (1,3,4), (3,0,4)
        let f1 = mesh.add_face(Face::triangle(0, 1, 3)); // External
        // f_014 is shared
        let f_134 = mesh.add_face(Face::triangle(1, 3, 4)); // Internal
        let f_304 = mesh.add_face(Face::triangle(3, 0, 4)); // Internal
        mesh.add_cell(Cell::tetrahedron(f1, f_014, f_134, f_304));

        // Tet 3: 0,2,3,4 (Base 0,2,3, Peak 4)
        // Faces: (0,2,3), (0,2,4)*, (2,3,4), (3,0,4)*
        let f2 = mesh.add_face(Face::triangle(0, 2, 3)); // External
        // f_204 is shared
        let f_234 = mesh.add_face(Face::triangle(2, 3, 4)); // Internal
        // f_304 is shared
        mesh.add_cell(Cell::tetrahedron(f2, f_204, f_234, f_304));

        // Tet 4: 1,2,3,4 (Base 1,2,3, Peak 4)
        // Faces: (1,2,3), (1,2,4)*, (2,3,4)*, (1,3,4)*
        let f3 = mesh.add_face(Face::triangle(1, 2, 3)); // External
        // f_124 is shared
        // f_234 is shared
        // f_134 is shared
        mesh.add_cell(Cell::tetrahedron(f3, f_124, f_234, f_134));

        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();
        let problem = StokesFlowProblem::new(mesh.clone(), fluid.clone(), boundary_conditions.clone(), 5);

        let boundary_nodes = problem.get_boundary_nodes();

        // Node 4 is interior. It should NOT be in boundary_nodes.
        assert!(!boundary_nodes.contains(&4), "Interior node 4 should not be boundary");
        assert!(boundary_nodes.contains(&0));
        assert!(boundary_nodes.contains(&1));
        assert!(boundary_nodes.contains(&2));
        assert!(boundary_nodes.contains(&3));

        // Now mark an internal face as boundary.
        // Pick f_014. It connects 0, 1, 4.
        let mut mesh_marked = mesh.clone();
        mesh_marked.mark_boundary(f_014, "internal_boundary".to_string());

        let problem_marked = StokesFlowProblem::new(mesh_marked, fluid, boundary_conditions, 5);
        let boundary_nodes_marked = problem_marked.get_boundary_nodes();

        // This assertion will FAIL with buggy implementation if we allow internal markers to add nodes
        // It should PASS after fix (ignoring internal markers).
        assert!(!boundary_nodes_marked.contains(&4), "Interior node 4 should not be included even if internal face is marked");
    }
}
