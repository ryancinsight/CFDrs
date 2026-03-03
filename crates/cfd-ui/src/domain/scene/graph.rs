//! Scene graph — tree of named scene nodes with transform hierarchy.

use crate::domain::scene::camera::OrbitalCamera;
use nalgebra::Isometry3;

/// Handle referencing a mesh stored in the project document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct MeshHandle(pub usize);

/// Handle referencing a schematic stored in the project document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct SchematicHandle(pub usize);

/// Handle referencing simulation results stored in the project document.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub struct ResultHandle(pub usize);

/// The payload carried by a scene node.
#[derive(Clone, Debug)]
pub enum SceneEntity {
    /// A triangulated mesh.
    Mesh(MeshHandle),
    /// A 2D schematic layout.
    Schematic(SchematicHandle),
    /// Simulation result field data.
    SimulationResult(ResultHandle),
    /// A grouping node with no geometry.
    Group,
}

/// A node in the scene graph.
#[derive(Clone, Debug)]
pub struct SceneNode {
    /// Human-readable name shown in the model tree.
    pub name: String,
    /// The entity this node represents.
    pub entity: SceneEntity,
    /// Local-to-parent transform.
    pub transform: Isometry3<f64>,
    /// Visibility toggle.
    pub visible: bool,
    /// Indices of child nodes in `SceneGraph::nodes`.
    pub children: Vec<usize>,
    /// Index of parent node (None for root).
    pub parent: Option<usize>,
}

/// A flat-array scene graph with parent/child relationships.
#[derive(Clone, Debug)]
pub struct SceneGraph {
    nodes: Vec<SceneNode>,
    root: usize,
    camera: OrbitalCamera,
}

impl SceneGraph {
    /// Create an empty scene graph with a root group node and default camera.
    #[must_use]
    pub fn new() -> Self {
        let root_node = SceneNode {
            name: "Root".to_owned(),
            entity: SceneEntity::Group,
            transform: Isometry3::identity(),
            visible: true,
            children: Vec::new(),
            parent: None,
        };
        Self {
            nodes: vec![root_node],
            root: 0,
            camera: OrbitalCamera::default(),
        }
    }

    /// Add a node as a child of the root. Returns the node index.
    pub fn add_node(&mut self, name: String, entity: SceneEntity) -> usize {
        let idx = self.nodes.len();
        self.nodes.push(SceneNode {
            name,
            entity,
            transform: Isometry3::identity(),
            visible: true,
            children: Vec::new(),
            parent: Some(self.root),
        });
        self.nodes[self.root].children.push(idx);
        idx
    }

    /// Remove a node by index. Orphans its children (they become invisible).
    pub fn remove_node(&mut self, idx: usize) {
        if idx == self.root || idx >= self.nodes.len() {
            return;
        }
        if let Some(parent) = self.nodes[idx].parent {
            self.nodes[parent].children.retain(|&c| c != idx);
        }
        self.nodes[idx].visible = false;
        self.nodes[idx].children.clear();
    }

    /// Immutable access to a node.
    #[must_use]
    pub fn node(&self, idx: usize) -> Option<&SceneNode> {
        self.nodes.get(idx)
    }

    /// Mutable access to a node.
    pub fn node_mut(&mut self, idx: usize) -> Option<&mut SceneNode> {
        self.nodes.get_mut(idx)
    }

    /// The root node index.
    #[must_use]
    pub fn root(&self) -> usize {
        self.root
    }

    /// The camera.
    #[must_use]
    pub fn camera(&self) -> &OrbitalCamera {
        &self.camera
    }

    /// Mutable camera access.
    pub fn camera_mut(&mut self) -> &mut OrbitalCamera {
        &mut self.camera
    }

    /// Iterate over all visible mesh handles in the scene.
    pub fn visible_meshes(&self) -> impl Iterator<Item = (usize, MeshHandle)> + '_ {
        self.nodes.iter().enumerate().filter_map(|(i, node)| {
            if node.visible {
                if let SceneEntity::Mesh(h) = &node.entity {
                    return Some((i, *h));
                }
            }
            None
        })
    }

    /// Number of nodes (including root).
    #[must_use]
    pub fn node_count(&self) -> usize {
        self.nodes.len()
    }

    /// Children of a node.
    #[must_use]
    pub fn children(&self, idx: usize) -> &[usize] {
        self.nodes
            .get(idx)
            .map(|n| n.children.as_slice())
            .unwrap_or(&[])
    }
}

impl Default for SceneGraph {
    fn default() -> Self {
        Self::new()
    }
}
