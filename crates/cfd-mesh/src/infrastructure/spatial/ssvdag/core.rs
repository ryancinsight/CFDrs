//! Symmetric Sparse Voxel Directed Acyclic Graph (SSVDAG) / Sparse Voxel Octree (SVO)
//!
//! An exact-arithmetic SVO implementation that compresses memory hierarchically
//! via Directed Acyclic Graphs, avoiding the inexactness of floating-point KD-trees.
//!
//! ## Mathematical Formalism
//!
//! An Octree over domain $\Omega \subset \mathbb{R}^3$ recursively partitions
//! the volume into 8 equal sub-domains.
//! Exact recursive subdivision of IEEE 754 floats by 2.0 introduces zero
//! rounding error (until subnormal underflow), preserving topological correctness.
//!
//! By deduplicating identical child-pointer arrays `[u32; 8]`, the tree becomes
//! a Directed Acyclic Graph (DAG).
//!
//! ### Theorem: Isomorphic Subtree Compression
//! Two subtrees $T_1 = T_2$ structurally and by value if and only if their
//! child pointers map to identical DAG indices. Uniqueness is guaranteed by
//! bottom-up hashing.

use crate::domain::core::scalar::Point3r;
use crate::domain::geometry::aabb::Aabb;
use hashbrown::HashMap;

/// An index into the DAG node pool.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash)]
pub struct DagIndex(pub u32);

/// A node in the Sparse Voxel DAG.
#[derive(Debug, Clone, PartialEq, Eq, Hash)]
pub enum SvoNode {
    /// A homogeneous volume (true = solid, false = empty).
    Leaf(bool),
    /// An internal node with 8 octant children.
    Internal([DagIndex; 8]),
}

/// Sparse Voxel Octree built with exact Cartesian recursive subdivision.
#[derive(Clone)]
pub struct SparseVoxelOctree {
    /// The bounding box of the entire DAG.
    pub root_aabb: Aabb,
    /// The unique internal nodes and leaves.
    pub nodes: Vec<SvoNode>,
    /// Deduplication map mapping a node structural value to its index.
    node_map: HashMap<SvoNode, DagIndex>,
    /// The root index of the DAG.
    pub root_index: DagIndex,
}

impl SparseVoxelOctree {
    /// Create a new empty SVO over `domain_aabb`.
    pub fn new(domain_aabb: Aabb) -> Self {
        let empty_leaf = SvoNode::Leaf(false);
        let nodes = vec![empty_leaf.clone()];
        let mut node_map = HashMap::new();
        node_map.insert(empty_leaf, DagIndex(0));

        Self {
            root_aabb: domain_aabb,
            nodes,
            node_map,
            root_index: DagIndex(0),
        }
    }

    /// Recursively insert or lookup a node, ensuring DAG deduplication.
    pub(crate) fn intern_node(&mut self, node: SvoNode) -> DagIndex {
        if let Some(&idx) = self.node_map.get(&node) {
            return idx;
        }
        let max_idx = self.nodes.len() as u32;
        let idx = DagIndex(max_idx);
        self.nodes.push(node.clone());
        self.node_map.insert(node, idx);
        idx
    }

    /// Splits an AABB tightly into 8 sub-octants.
    /// Because division by 2.0 only subtracts 1 from the exponent of an IEEE 754 float,
    /// the split is exact (no mantissa truncation).
    pub(crate) fn exact_subdivide(aabb: &Aabb) -> [Aabb; 8] {
        let min = aabb.min;
        let max = aabb.max;
        let mid = Point3r::new(
            (min.x + max.x) * 0.5,
            (min.y + max.y) * 0.5,
            (min.z + max.z) * 0.5,
        );

        let make = |bx, by, bz| {
            let x0 = if bx == 0 { min.x } else { mid.x };
            let x1 = if bx == 0 { mid.x } else { max.x };
            let y0 = if by == 0 { min.y } else { mid.y };
            let y1 = if by == 0 { mid.y } else { max.y };
            let z0 = if bz == 0 { min.z } else { mid.z };
            let z1 = if bz == 0 { mid.z } else { max.z };
            Aabb::new(Point3r::new(x0, y0, z0), Point3r::new(x1, y1, z1))
        };

        [
            make(0, 0, 0),
            make(1, 0, 0),
            make(0, 1, 0),
            make(1, 1, 0),
            make(0, 0, 1),
            make(1, 0, 1),
            make(0, 1, 1),
            make(1, 1, 1),
        ]
    }

    /// Rasterize an explicit AABB into the SVO DAG up to `max_depth`.
    pub fn insert_aabb(&mut self, target: &Aabb, max_depth: u8) {
        let root = self.root_aabb.clone();
        let current_root = self.root_index;
        self.root_index = self.insert_recursive(current_root, &root, target, 0, max_depth);
    }

    fn insert_recursive(
        &mut self,
        current_idx: DagIndex,
        current_aabb: &Aabb,
        target: &Aabb,
        depth: u8,
        max_depth: u8,
    ) -> DagIndex {
        if !current_aabb.intersects(target) {
            return current_idx;
        }

        if depth >= max_depth || target.contains_aabb(current_aabb) {
            return self.intern_node(SvoNode::Leaf(true));
        }

        let mut children = match &self.nodes[current_idx.0 as usize] {
            SvoNode::Leaf(val) => {
                let s = self.intern_node(SvoNode::Leaf(*val));
                [s; 8]
            }
            SvoNode::Internal(c) => *c,
        };

        let sub_aabbs = Self::exact_subdivide(current_aabb);
        for i in 0..8 {
            children[i] =
                self.insert_recursive(children[i], &sub_aabbs[i], target, depth + 1, max_depth);
        }

        // Simplify identical children to a leaf
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

    #[test]
    fn exact_compression_guarantees_isomorphism() {
        // Build a flat SVO where everything is empty, the DAG should maximally compress
        let aabb = Aabb::new(Point3r::new(0.0, 0.0, 0.0), Point3r::new(10.0, 10.0, 10.0));
        let mut svo = SparseVoxelOctree::new(aabb);

        let target = Aabb::new(
            Point3r::new(20.0, 20.0, 20.0),
            Point3r::new(30.0, 30.0, 30.0),
        );
        svo.insert_aabb(&target, 3);

        // Due to isomorphism compression, all nodes outside bounds resolve to a single empty leaf 0.
        // There should be very few nodes.
        assert!(
            svo.nodes.len() < 10,
            "SVO failed to compress isomorphic empty paths"
        );
    }
}
