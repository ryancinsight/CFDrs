//! Model tree panel — hierarchical scene browser.

use crate::domain::document::model_tree::ModelTreeNode;
use crate::domain::scene::graph::SceneGraph;

/// Build a model tree from the scene graph for UI display.
#[must_use]
pub fn build_model_tree(scene: &SceneGraph) -> Vec<ModelTreeNode> {
    let root = scene.root();
    scene
        .children(root)
        .iter()
        .filter_map(|&idx| {
            let node = scene.node(idx)?;
            Some(build_subtree(scene, idx, node.name.clone()))
        })
        .collect()
}

fn build_subtree(scene: &SceneGraph, idx: usize, label: String) -> ModelTreeNode {
    let children = scene
        .children(idx)
        .iter()
        .filter_map(|&child_idx| {
            let child = scene.node(child_idx)?;
            Some(build_subtree(scene, child_idx, child.name.clone()))
        })
        .collect();

    if scene.children(idx).is_empty() {
        ModelTreeNode::leaf(idx, label)
    } else {
        ModelTreeNode::group(idx, label, children)
    }
}
