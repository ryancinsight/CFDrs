//! Exact Volumetric Rasterization
//!
//! Converts a boundary-representation `IndexedMesh` into a Symmetric Sparse Voxel DAG (SSVDAG).
//!
//! ## Algorithm
//!
//! The octree is built by pure depth-first subdivision: every cell with `depth < max_depth`
//! is unconditionally split into 8 sub-octants.  At the leaf level (`depth == max_depth`) the
//! exact Generalized Winding Number (GWN) evaluated at the cell centre determines whether the
//! cell is solid or empty.
//!
//! An early-exit optimisation prunes cells that do not overlap the mesh's axis-aligned
//! bounding box at all: such cells cannot be interior to any closed surface, so GWN = 0.
//!
//! ## Correctness Guarantees
//!
//! * **Topological correctness** — GWN is ±1 inside a watertight 2-manifold and 0 outside
//!   (Jacobson et al. 2013).  No floating-point epsilon displacements are required.
//! * **Isomorphic DAG compression** — identical subtrees share a single node, giving O(depth)
//!   memory for uniform regions rather than O(8^depth).
//! * **Exact subdivision** — dividing an IEEE 754 float by 2.0 only decrements the exponent
//!   (no mantissa rounding until subnormal underflow), preserving geometric fidelity.

use crate::application::csg::arrangement::classify::gwn;
use crate::domain::core::scalar::{Point3r, Scalar};
use crate::domain::geometry::aabb::Aabb;
use crate::domain::mesh::indexed::IndexedMesh;

use super::core::{DagIndex, SparseVoxelOctree, SvoNode};

impl SparseVoxelOctree {
    /// Exact volumetric rasterization of an `IndexedMesh` into a new `SparseVoxelOctree`.
    ///
    /// The resulting octree is built recursively up to `max_depth` levels.
    /// `mesh` should be a closed, watertight 2-manifold for correct GWN inside/outside
    /// classification.
    ///
    /// # Arguments
    ///
    /// * `mesh`      — The boundary-representation surface mesh to rasterize.
    /// * `max_depth` — Maximum octree depth.  Leaf voxel edge length ≈ mesh_extent / 2^max_depth.
    pub fn from_mesh<T: Scalar>(mesh: &IndexedMesh<T>, max_depth: u8) -> Self {
        // Degenerate mesh guard.
        let root_aabb_t = mesh.bounding_box();
        if root_aabb_t.volume() == <T as Scalar>::from_f64(0.0) {
            return Self::new(Aabb::new(
                Point3r::new(-1.0, -1.0, -1.0),
                Point3r::new(1.0, 1.0, 1.0),
            ));
        }

        // Helper: losslessly convert a T scalar to f64.
        let to_f = |v: T| v.to_f64().unwrap_or(0.0);

        // Expand root by 1% to avoid exactly splitting surface triangles at the octree
        // root boundary (which would make the root-level GWN ambiguous).
        let dx = to_f((root_aabb_t.max.x - root_aabb_t.min.x) * <T as Scalar>::from_f64(0.01));
        let dy = to_f((root_aabb_t.max.y - root_aabb_t.min.y) * <T as Scalar>::from_f64(0.01));
        let dz = to_f((root_aabb_t.max.z - root_aabb_t.min.z) * <T as Scalar>::from_f64(0.01));

        // Unpadded mesh bounding box in f64 (for the early-exit optimisation).
        let mesh_bb = Aabb::new(
            Point3r::new(to_f(root_aabb_t.min.x), to_f(root_aabb_t.min.y), to_f(root_aabb_t.min.z)),
            Point3r::new(to_f(root_aabb_t.max.x), to_f(root_aabb_t.max.y), to_f(root_aabb_t.max.z)),
        );

        // Padded root domain for the octree.
        let expanded_root = Aabb::new(
            Point3r::new(to_f(root_aabb_t.min.x) - dx, to_f(root_aabb_t.min.y) - dy, to_f(root_aabb_t.min.z) - dz),
            Point3r::new(to_f(root_aabb_t.max.x) + dx, to_f(root_aabb_t.max.y) + dy, to_f(root_aabb_t.max.z) + dz),
        );

        let mut svo = Self::new(expanded_root.clone());
        let new_root = svo.rasterize_recursive(mesh, &expanded_root, &mesh_bb, 0, max_depth);
        svo.root_index = new_root;
        svo
    }

    /// Recursive octree subdivision with exact GWN leaf classification.
    ///
    /// ## Design invariants
    ///
    /// 1. **Correctness before speed**: every cell with `depth < max_depth` is subdivided
    ///    unconditionally.  The only early exit is for cells that lie entirely outside the
    ///    mesh bounding box (guaranteed to be exterior).
    /// 2. **Leaf classification**: GWN at the cell centre determines solid (`|wn| > 0.5`)
    ///    or empty (`|wn| ≤ 0.5`).  Van Oosterom–Strackee solid-angle formula (1983).
    /// 3. **Isomorphic compression**: 8 identical children collapse to a single Leaf node.
    fn rasterize_recursive<T: Scalar>(
        &mut self,
        mesh: &IndexedMesh<T>,
        current_aabb: &Aabb,
        mesh_bb: &Aabb,
        depth: u8,
        max_depth: u8,
    ) -> DagIndex {
        // Early exit: cell entirely outside the mesh bounding box → definitely empty.
        if !current_aabb.intersects(mesh_bb) {
            return self.intern_node(SvoNode::Leaf(false));
        }

        // Leaf level: evaluate exact GWN at the cell centre.
        if depth >= max_depth {
            let center = current_aabb.center();
            let p = nalgebra::Point3::new(
                <T as Scalar>::from_f64(center.x),
                <T as Scalar>::from_f64(center.y),
                <T as Scalar>::from_f64(center.z),
            );
            let wn = gwn(&p, mesh.faces.as_slice(), &mesh.vertices);
            let is_solid = wn.to_f64().unwrap_or(0.0).abs() > 0.5;
            return self.intern_node(SvoNode::Leaf(is_solid));
        }

        // Internal: unconditionally subdivide into 8 equal sub-octants.
        let sub_aabbs = Self::exact_subdivide(current_aabb);
        let mut children = [DagIndex(0); 8];
        for i in 0..8 {
            children[i] =
                self.rasterize_recursive(mesh, &sub_aabbs[i], mesh_bb, depth + 1, max_depth);
        }

        // Isomorphic DAG compression: all-same-Leaf children collapse to one Leaf.
        if children.iter().all(|&c| c == children[0]) {
            if let SvoNode::Leaf(v) = self.nodes[children[0].0 as usize] {
                return self.intern_node(SvoNode::Leaf(v));
            }
        }

        self.intern_node(SvoNode::Internal(children))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::geometry::primitives::cube::Cube;
    use crate::domain::geometry::primitives::PrimitiveMesh;

    /// At max_depth=3 the expanded root is (±0.51)³.
    /// Each leaf side = 1.02/8 = 0.1275 mm.
    /// Corner leaf centre = 0.51 − 0.1275/2 = 0.446 < 0.5  → inside cube.
    /// Therefore ALL 512 potential leaves are classified solid, and isomorphic
    /// compression collapses the entire DAG to a single Leaf(true).
    #[test]
    fn shallow_rasterize_compresses_to_solid_leaf() {
        let mesh = Cube::centred(1.0).build().unwrap();
        let svo = SparseVoxelOctree::from_mesh(&mesh, 3);

        println!("depth=3 node count: {}", svo.nodes.len());
        match &svo.nodes[svo.root_index.0 as usize] {
            SvoNode::Leaf(solid) => {
                assert!(*solid, "All depth-3 leaf centres are inside the cube → must be solid");
            }
            SvoNode::Internal(_) => {
                panic!("At depth=3 every leaf centre is inside the cube; \
                        isomorphic compression should yield a single Leaf(true)");
            }
        }
    }

    /// At max_depth=6 the expanded root is (±0.51)³.
    /// Each leaf side = 1.02/64 = 0.015938 mm.
    /// Corner leaf centre = 0.51 − 0.015938/2 ≈ 0.502 > 0.5  → OUTSIDE cube.
    /// Therefore at least one corner leaf is empty, and the root must be Internal.
    ///
    /// This test proves that subdivision is actually occurring (not just evaluating
    /// GWN at the root centre) and that exterior cells are correctly classified.
    #[test]
    fn deep_rasterize_has_interior_root_with_empty_corners() {
        let mesh = Cube::centred(1.0).build().unwrap();
        let svo = SparseVoxelOctree::from_mesh(&mesh, 6);

        println!("depth=6 node count: {}", svo.nodes.len());
        match &svo.nodes[svo.root_index.0 as usize] {
            SvoNode::Internal(children) => {
                assert_eq!(children.len(), 8, "Internal node must have exactly 8 children");
                // Verify DAG compression is working: node count must be far below 8^6=262144
                assert!(
                    svo.nodes.len() < 1000,
                    "DAG compression should deduplicate identical subtrees; got {} nodes",
                    svo.nodes.len()
                );
            }
            SvoNode::Leaf(_) => {
                panic!(
                    "At depth=6, corner leaf centres (≈0.502) are outside the cube; \
                     the root must be Internal with some empty leaf descendants"
                );
            }
        }
    }

    /// GWN sanity check: interior point gives |wn| ≈ 1, exterior gives |wn| ≈ 0.
    #[test]
    fn gwn_interior_exterior_classification() {
        let mesh = Cube::centred(1.0).build().unwrap();

        // Interior: (0,0,0) is the cube centre → winding number ≈ +1.
        let p_in = nalgebra::Point3::new(0.0_f64, 0.0, 0.0);
        let wn_in = gwn(&p_in, mesh.faces.as_slice(), &mesh.vertices);
        assert!(
            wn_in.abs() > 0.5,
            "GWN at cube centre must be near ±1, got {wn_in}"
        );

        // Exterior: (5,5,5) is far outside → winding number ≈ 0.
        let p_out = nalgebra::Point3::new(5.0_f64, 5.0, 5.0);
        let wn_out = gwn(&p_out, mesh.faces.as_slice(), &mesh.vertices);
        assert!(
            wn_out.abs() < 0.1,
            "GWN at far exterior must be near 0, got {wn_out}"
        );
    }
}
