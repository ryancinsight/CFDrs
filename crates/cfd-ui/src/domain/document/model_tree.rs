//! Model tree node — hierarchical view of the scene for UI display.

/// A node in the model tree view.
#[derive(Clone, Debug)]
pub struct ModelTreeNode {
    /// Scene graph node index this tree node references.
    pub scene_index: usize,
    /// Display name in the tree.
    pub label: String,
    /// Whether this node is expanded in the tree view.
    pub expanded: bool,
    /// Child tree nodes.
    pub children: Vec<ModelTreeNode>,
}

impl ModelTreeNode {
    /// Create a leaf tree node.
    #[must_use]
    pub fn leaf(scene_index: usize, label: String) -> Self {
        Self {
            scene_index,
            label,
            expanded: false,
            children: Vec::new(),
        }
    }

    /// Create a group tree node.
    #[must_use]
    pub fn group(scene_index: usize, label: String, children: Vec<Self>) -> Self {
        Self {
            scene_index,
            label,
            expanded: true,
            children,
        }
    }
}
