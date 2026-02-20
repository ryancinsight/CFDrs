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
    /// Boundary nodes are vertices that belong to TRUE boundary faces.
    /// A true boundary face is one that:
    /// 1. Is referenced by exactly one cell (topological boundary).
    /// 2. Does NOT have a geometric duplicate (internal interface with non-unified topology).
    ///
    /// # Algorithm
    ///
    /// 1. Identify topological boundary faces (count == 1).
    /// 2. Filter out faces that are geometrically identical (duplicate faces).
    ///    This handles disjoint meshes where internal faces are duplicated.
    /// 3. Collect unique vertices from remaining faces.
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

        // 1. Identify topological boundary faces (referenced by exactly one cell)
        let mut face_cell_count: HashMap<usize, usize> = HashMap::new();
        for cell in self.mesh.cells() {
            for &face_idx in &cell.faces {
                *face_cell_count.entry(face_idx).or_insert(0) += 1;
            }
        }

        let topo_boundary_faces: Vec<usize> = face_cell_count
            .into_iter()
            .filter(|&(_, count)| count == 1)
            .map(|(face_idx, _)| face_idx)
            .collect();

        // 2. Identify and filter duplicate faces (geometric coincidence)
        // Count how many times each set of vertices appears among the topological boundary faces.
        // If a face appears >= 2 times (e.g., one for Cell A, one for Cell B), it is an internal interface.
        let mut geom_face_count: HashMap<Vec<usize>, usize> = HashMap::new();

        for &f_idx in &topo_boundary_faces {
            if let Some(face) = self.mesh.face(f_idx) {
                let mut v = face.vertices.clone();
                v.sort_unstable(); // Sort to ensure consistent key for geometric comparison
                *geom_face_count.entry(v).or_insert(0) += 1;
            }
        }

        // Collect true boundary faces (count == 1)
        let mut true_boundary_faces: HashSet<usize> = HashSet::new();

        // Include explicitly marked faces? (Currently, marked faces are usually subsets of topo_boundary_faces)
        // If we want to support internal baffles (marked internal faces), we should check markers.
        // But for standard validation, we assume we only care about the hull.
        // For robustness, if a face is marked, we might consider it a boundary regardless?
        // Let's stick to the geometric hull logic to fix the "interior node" bug.

        for &f_idx in &topo_boundary_faces {
            if let Some(face) = self.mesh.face(f_idx) {
                let mut v = face.vertices.clone();
                v.sort_unstable();
                // Only include if it appears exactly once geometrically
                if geom_face_count.get(&v) == Some(&1) {
                    true_boundary_faces.insert(f_idx);
                }
            }
        }

        // 3. Collect unique vertices from all TRUE boundary faces
        let mut boundary_vertices: HashSet<usize> = HashSet::new();
        for &face_idx in &true_boundary_faces {
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
    fn test_internal_node_is_not_boundary_with_unique_faces() {
        use cfd_mesh::topology::Vertex;
        // Construct a mesh with a central node (0) surrounded by 4 nodes (1,2,3,4).
        // This forms a large tetrahedron subdivided into 4 smaller tets meeting at 0.
        // Node 0 is STRICTLY internal. It should NOT be returned by get_boundary_nodes().

        let mut mesh = Mesh::new();
        mesh.add_vertex(Vertex::new(Point3::new(0.0, 0.0, 0.0))); // 0: Center
        mesh.add_vertex(Vertex::new(Point3::new(1.0, 1.0, 1.0))); // 1
        mesh.add_vertex(Vertex::new(Point3::new(1.0, -1.0, -1.0))); // 2
        mesh.add_vertex(Vertex::new(Point3::new(-1.0, 1.0, -1.0))); // 3
        mesh.add_vertex(Vertex::new(Point3::new(-1.0, -1.0, 1.0))); // 4

        // Faces involving 0 (Internal) - UNIQUE FACES (Shared)
        let f_012 = mesh.add_face(Face::triangle(0, 1, 2));
        let f_023 = mesh.add_face(Face::triangle(0, 2, 3));
        let f_031 = mesh.add_face(Face::triangle(0, 3, 1));
        let f_014 = mesh.add_face(Face::triangle(0, 1, 4));
        let f_042 = mesh.add_face(Face::triangle(0, 4, 2)); // 2-4 edge
        let f_034 = mesh.add_face(Face::triangle(0, 3, 4));

        // Faces on the hull (External)
        let f_123 = mesh.add_face(Face::triangle(1, 2, 3));
        let f_142 = mesh.add_face(Face::triangle(1, 4, 2));
        let f_134 = mesh.add_face(Face::triangle(1, 3, 4));
        let f_234 = mesh.add_face(Face::triangle(2, 3, 4));

        // Cells
        // Tet 1: (0,1,2,3) uses f_012, f_023, f_031, f_123
        mesh.add_cell(Cell::tetrahedron(f_012, f_023, f_031, f_123));

        // Tet 2: (0,1,4,2) uses f_014, f_042, f_012 (shared), f_142
        mesh.add_cell(Cell::tetrahedron(f_014, f_042, f_012, f_142));

        // Tet 3: (0,1,3,4) uses f_031 (shared), f_014 (shared), f_034, f_134
        mesh.add_cell(Cell::tetrahedron(f_031, f_014, f_034, f_134));

        // Tet 4: (0,2,3,4) uses f_023 (shared), f_034 (shared), f_042 (shared), f_234
        mesh.add_cell(Cell::tetrahedron(f_023, f_034, f_042, f_234));

        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let boundary_conditions = HashMap::new();

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 5);
        let boundary_nodes = problem.get_boundary_nodes();

        // Node 0 is internal. It should NOT be in boundary_nodes.
        assert!(!boundary_nodes.contains(&0), "Node 0 (internal) should not be flagged as boundary");

        // Nodes 1,2,3,4 are external. They MUST be in boundary_nodes.
        assert!(boundary_nodes.contains(&1));
        assert!(boundary_nodes.contains(&2));
        assert!(boundary_nodes.contains(&3));
        assert!(boundary_nodes.contains(&4));
    }

    #[test]
    fn test_internal_node_is_boundary_with_duplicate_faces_fix() {
        use cfd_mesh::topology::Vertex;
        // Same geometry, but internal faces are DUPLICATED.
        // This simulates a mesh where shared faces are not unified.
        // BEFORE FIX: Node 0 was flagged as boundary because duplicate faces had count=1.
        // AFTER FIX: Node 0 should NOT be flagged as boundary because duplicate faces are detected geometrically.

        let mut mesh = Mesh::new();
        mesh.add_vertex(Vertex::new(Point3::new(0.0, 0.0, 0.0))); // 0
        mesh.add_vertex(Vertex::new(Point3::new(1.0, 1.0, 1.0))); // 1
        mesh.add_vertex(Vertex::new(Point3::new(1.0, -1.0, -1.0))); // 2
        mesh.add_vertex(Vertex::new(Point3::new(-1.0, 1.0, -1.0))); // 3
        mesh.add_vertex(Vertex::new(Point3::new(-1.0, -1.0, 1.0))); // 4

        // Internal faces - DUPLICATED for each tet
        // Tet 1 internal faces
        let f_012_t1 = mesh.add_face(Face::triangle(0, 1, 2));
        let f_023_t1 = mesh.add_face(Face::triangle(0, 2, 3));
        let f_031_t1 = mesh.add_face(Face::triangle(0, 3, 1));

        // Tet 2 internal faces
        let f_014_t2 = mesh.add_face(Face::triangle(0, 1, 4));
        let f_042_t2 = mesh.add_face(Face::triangle(0, 4, 2));
        let f_012_t2 = mesh.add_face(Face::triangle(0, 1, 2)); // Duplicate of f_012_t1

        // Tet 3 internal faces
        let f_031_t3 = mesh.add_face(Face::triangle(0, 3, 1)); // Duplicate of f_031_t1
        let f_014_t3 = mesh.add_face(Face::triangle(0, 1, 4)); // Duplicate of f_014_t2
        let f_034_t3 = mesh.add_face(Face::triangle(0, 3, 4));

        // Tet 4 internal faces
        let f_023_t4 = mesh.add_face(Face::triangle(0, 2, 3)); // Duplicate of f_023_t1
        let f_034_t4 = mesh.add_face(Face::triangle(0, 3, 4)); // Duplicate of f_034_t3
        let f_042_t4 = mesh.add_face(Face::triangle(0, 4, 2)); // Duplicate of f_042_t2

        // External faces (count=1 naturally)
        let f_123 = mesh.add_face(Face::triangle(1, 2, 3));
        let f_142 = mesh.add_face(Face::triangle(1, 4, 2));
        let f_134 = mesh.add_face(Face::triangle(1, 3, 4));
        let f_234 = mesh.add_face(Face::triangle(2, 3, 4));

        // Cells
        mesh.add_cell(Cell::tetrahedron(f_012_t1, f_023_t1, f_031_t1, f_123));
        mesh.add_cell(Cell::tetrahedron(f_014_t2, f_042_t2, f_012_t2, f_142));
        mesh.add_cell(Cell::tetrahedron(f_031_t3, f_014_t3, f_034_t3, f_134));
        mesh.add_cell(Cell::tetrahedron(f_023_t4, f_034_t4, f_042_t4, f_234));

        let fluid = ConstantPropertyFluid::<f64>::water_20c().unwrap();
        let mut boundary_conditions = HashMap::new();

        // Assign BCs to external nodes (1,2,3,4) so we can pass validate()
        for i in 1..=4 {
            boundary_conditions.insert(i, BoundaryCondition::Dirichlet { value: 0.0, component_values: None });
        }

        let problem = StokesFlowProblem::new(mesh, fluid, boundary_conditions, 5);
        let boundary_nodes = problem.get_boundary_nodes();

        // Node 0 should NOT be flagged as boundary, because its faces are geometrically paired.
        assert!(!boundary_nodes.contains(&0), "Node 0 (internal) should NOT be flagged as boundary despite duplicate faces");

        // External nodes should be flagged
        assert!(boundary_nodes.contains(&1));
        assert!(boundary_nodes.contains(&2));
        assert!(boundary_nodes.contains(&3));
        assert!(boundary_nodes.contains(&4));

        // Validation should now PASS
        assert!(problem.validate().is_ok());
    }
}
