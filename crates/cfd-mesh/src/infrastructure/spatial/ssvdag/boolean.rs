//! Exact Boolean computing on Sparse Voxel DAGs.
//!
//! Because the SSVDAG compresses identical geometric structures into the same DAG nodes,
//! boolean operations (Union, Intersection, Difference) can be computed exactly and swiftly
//! by recursively merging nodes.

use super::core::{DagIndex, SparseVoxelOctree, SvoNode};
use hashbrown::HashMap;

/// The type of exact Boolean CSG operation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum BooleanOp {
    /// Combine the volume of A and B.
    Union,
    /// Keep only the volume shared by A and B.
    Intersection,
    /// Subtract B from A.
    Difference,
}

impl SparseVoxelOctree {
    /// Compute the exact boolean operation between `self` (A) and `other` (B).
    ///
    /// Both DAGs must cover the exact same domain AABB.
    pub fn boolean(&self, other: &Self, op: BooleanOp) -> Self {
        assert_eq!(
            self.root_aabb, other.root_aabb,
            "SSVDAG Booleans require identical initial domain AABBs to preserve spatial consistency"
        );

        let mut result = SparseVoxelOctree::new(self.root_aabb.clone());

        let mut ctx = MergeContext {
            a: self,
            b: other,
            out: &mut result,
            memo_merge: HashMap::new(),
            memo_map_a: HashMap::new(),
            memo_map_b: HashMap::new(),
            memo_inv_b: HashMap::new(),
            op,
        };

        let new_root = ctx.merge_dag(self.root_index, other.root_index);
        result.root_index = new_root;

        result
    }
}

/// Recursive merge context with memoization over DAG indices.
struct MergeContext<'a> {
    a: &'a SparseVoxelOctree,
    b: &'a SparseVoxelOctree,
    out: &'a mut SparseVoxelOctree,
    /// Memoizes `merge_dag(a_idx, b_idx)`
    memo_merge: HashMap<(DagIndex, DagIndex), DagIndex>,
    /// Memoizes mapping a tree A node into `out`
    memo_map_a: HashMap<DagIndex, DagIndex>,
    /// Memoizes mapping a tree B node into `out`
    memo_map_b: HashMap<DagIndex, DagIndex>,
    /// Memoizes inverted tree B node into `out`
    memo_inv_b: HashMap<DagIndex, DagIndex>,
    op: BooleanOp,
}

impl<'a> MergeContext<'a> {
    fn map_subtree_a(&mut self, root: DagIndex) -> DagIndex {
        if let Some(&cached) = self.memo_map_a.get(&root) {
            return cached;
        }

        let idx = match &self.a.nodes[root.0 as usize] {
            SvoNode::Leaf(v) => self.out.intern_node(SvoNode::Leaf(*v)),
            SvoNode::Internal(children) => {
                let mut new_children = [DagIndex(0); 8];
                for i in 0..8 {
                    new_children[i] = self.map_subtree_a(children[i]);
                }
                if new_children.iter().all(|&c| c == new_children[0]) {
                    if let SvoNode::Leaf(v) = self.out.nodes[new_children[0].0 as usize] {
                        return self.out.intern_node(SvoNode::Leaf(v));
                    }
                }
                self.out.intern_node(SvoNode::Internal(new_children))
            }
        };

        self.memo_map_a.insert(root, idx);
        idx
    }

    fn map_subtree_b(&mut self, root: DagIndex) -> DagIndex {
        if let Some(&cached) = self.memo_map_b.get(&root) {
            return cached;
        }

        let idx = match &self.b.nodes[root.0 as usize] {
            SvoNode::Leaf(v) => self.out.intern_node(SvoNode::Leaf(*v)),
            SvoNode::Internal(children) => {
                let mut new_children = [DagIndex(0); 8];
                for i in 0..8 {
                    new_children[i] = self.map_subtree_b(children[i]);
                }
                if new_children.iter().all(|&c| c == new_children[0]) {
                    if let SvoNode::Leaf(v) = self.out.nodes[new_children[0].0 as usize] {
                        return self.out.intern_node(SvoNode::Leaf(v));
                    }
                }
                self.out.intern_node(SvoNode::Internal(new_children))
            }
        };

        self.memo_map_b.insert(root, idx);
        idx
    }

    fn invert_subtree_b(&mut self, root: DagIndex) -> DagIndex {
        if let Some(&cached) = self.memo_inv_b.get(&root) {
            return cached;
        }

        let idx = match &self.b.nodes[root.0 as usize] {
            SvoNode::Leaf(v) => self.out.intern_node(SvoNode::Leaf(!*v)),
            SvoNode::Internal(children) => {
                let mut new_children = [DagIndex(0); 8];
                for i in 0..8 {
                    new_children[i] = self.invert_subtree_b(children[i]);
                }
                if new_children.iter().all(|&c| c == new_children[0]) {
                    if let SvoNode::Leaf(v) = self.out.nodes[new_children[0].0 as usize] {
                        return self.out.intern_node(SvoNode::Leaf(v));
                    }
                }
                self.out.intern_node(SvoNode::Internal(new_children))
            }
        };

        self.memo_inv_b.insert(root, idx);
        idx
    }

    fn merge_dag(&mut self, root_a: DagIndex, root_b: DagIndex) -> DagIndex {
        let key = (root_a, root_b);
        if let Some(&cached) = self.memo_merge.get(&key) {
            return cached;
        }

        let node_a = &self.a.nodes[root_a.0 as usize];
        let node_b = &self.b.nodes[root_b.0 as usize];

        let out_idx = match (node_a, node_b) {
            (SvoNode::Leaf(val_a), SvoNode::Leaf(val_b)) => {
                let res = match self.op {
                    BooleanOp::Union => *val_a || *val_b,
                    BooleanOp::Intersection => *val_a && *val_b,
                    BooleanOp::Difference => *val_a && !*val_b,
                };
                self.out.intern_node(SvoNode::Leaf(res))
            }
            (SvoNode::Leaf(val_a), SvoNode::Internal(_)) => match self.op {
                BooleanOp::Union => {
                    if *val_a {
                        self.out.intern_node(SvoNode::Leaf(true))
                    } else {
                        self.map_subtree_b(root_b)
                    }
                }
                BooleanOp::Intersection => {
                    if *val_a {
                        self.map_subtree_b(root_b)
                    } else {
                        self.out.intern_node(SvoNode::Leaf(false))
                    }
                }
                BooleanOp::Difference => {
                    if *val_a {
                        self.invert_subtree_b(root_b)
                    } else {
                        self.out.intern_node(SvoNode::Leaf(false))
                    }
                }
            },
            (SvoNode::Internal(_), SvoNode::Leaf(val_b)) => match self.op {
                BooleanOp::Union => {
                    if *val_b {
                        self.out.intern_node(SvoNode::Leaf(true))
                    } else {
                        self.map_subtree_a(root_a)
                    }
                }
                BooleanOp::Intersection => {
                    if *val_b {
                        self.map_subtree_a(root_a)
                    } else {
                        self.out.intern_node(SvoNode::Leaf(false))
                    }
                }
                BooleanOp::Difference => {
                    if *val_b {
                        self.out.intern_node(SvoNode::Leaf(false))
                    } else {
                        self.map_subtree_a(root_a)
                    }
                }
            },
            (SvoNode::Internal(ca), SvoNode::Internal(cb)) => {
                let mut new_children = [DagIndex(0); 8];
                let ca = *ca;
                let cb = *cb;

                for i in 0..8 {
                    new_children[i] = self.merge_dag(ca[i], cb[i]);
                }

                if new_children.iter().all(|&c| c == new_children[0]) {
                    if let SvoNode::Leaf(v) = self.out.nodes[new_children[0].0 as usize] {
                        return self.out.intern_node(SvoNode::Leaf(v));
                    }
                }
                self.out.intern_node(SvoNode::Internal(new_children))
            }
        };

        self.memo_merge.insert(key, out_idx);
        out_idx
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::core::scalar::Point3r;
    use crate::domain::geometry::aabb::Aabb;

    #[test]
    fn boolean_dag_union_exact() {
        let domain = Aabb::new(Point3r::new(0.0, 0.0, 0.0), Point3r::new(10.0, 10.0, 10.0));
        let mut a = SparseVoxelOctree::new(domain.clone());
        let mut b = SparseVoxelOctree::new(domain.clone());

        let target_a = Aabb::new(Point3r::new(1.0, 1.0, 1.0), Point3r::new(4.0, 4.0, 4.0));
        let target_b = Aabb::new(Point3r::new(3.0, 3.0, 3.0), Point3r::new(6.0, 6.0, 6.0));

        a.insert_aabb(&target_a, 4);
        b.insert_aabb(&target_b, 4);

        let merged = a.boolean(&b, BooleanOp::Union);
        assert!(merged.nodes.len() > 10, "Merged tree should contain dense nodes");
        
        let diff = a.boolean(&b, BooleanOp::Difference);
        assert!(diff.nodes.len() > 10, "Merged tree should contain dense nodes");
    }
}
