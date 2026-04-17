//! Project document — top-level container for a design session.

use cfd_mesh::IndexedMesh;
use crate::domain::scene::graph::{MeshHandle, SceneGraph};
use serde::{Deserialize, Serialize};

/// Metadata associated with a saved project.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct ProjectMetadata {
    /// Human-readable project name.
    pub name: String,
    /// Optional description.
    pub description: String,
    /// Author name.
    pub author: String,
}

impl Default for ProjectMetadata {
    fn default() -> Self {
        Self {
            name: "Untitled Project".to_owned(),
            description: String::new(),
            author: String::new(),
        }
    }
}

/// The project document holds all state for a single design session.
pub struct ProjectDocument {
    /// Metadata about this project.
    pub metadata: ProjectMetadata,
    /// The scene graph containing all nodes.
    pub scene: SceneGraph,
    /// Arena-allocated meshes referenced by `MeshHandle`.
    meshes: Vec<Option<IndexedMesh<f64>>>,
    /// Whether the document has unsaved changes.
    dirty: bool,
}

impl ProjectDocument {
    /// Create a new empty project.
    #[must_use]
    pub fn new() -> Self {
        Self {
            metadata: ProjectMetadata::default(),
            scene: SceneGraph::new(),
            meshes: Vec::new(),
            dirty: false,
        }
    }

    /// Store a mesh and return its handle.
    pub fn add_mesh(&mut self, mesh: IndexedMesh<f64>) -> MeshHandle {
        let handle = MeshHandle(self.meshes.len());
        self.meshes.push(Some(mesh));
        self.dirty = true;
        handle
    }

    /// Remove a mesh by handle (sets slot to None).
    pub fn remove_mesh(&mut self, handle: MeshHandle) {
        if handle.0 < self.meshes.len() {
            self.meshes[handle.0] = None;
            self.dirty = true;
        }
    }

    /// Access a mesh by handle.
    #[must_use]
    pub fn mesh(&self, handle: MeshHandle) -> Option<&IndexedMesh<f64>> {
        self.meshes.get(handle.0).and_then(|m| m.as_ref())
    }

    /// Mutable access to a mesh by handle.
    pub fn mesh_mut(&mut self, handle: MeshHandle) -> Option<&mut IndexedMesh<f64>> {
        self.dirty = true;
        self.meshes.get_mut(handle.0).and_then(|m| m.as_mut())
    }

    /// Number of live (non-removed) meshes.
    #[must_use]
    pub fn mesh_count(&self) -> usize {
        self.meshes.iter().filter(|m| m.is_some()).count()
    }

    /// Whether there are unsaved changes.
    #[must_use]
    pub fn is_dirty(&self) -> bool {
        self.dirty
    }

    /// Mark the document as saved.
    pub fn mark_saved(&mut self) {
        self.dirty = false;
    }

    /// Mark the document as modified.
    pub fn mark_dirty(&mut self) {
        self.dirty = true;
    }

    /// Iterate over all live mesh slots.
    pub(crate) fn mesh_slots(&self) -> impl Iterator<Item = (usize, &IndexedMesh<f64>)> {
        self.meshes
            .iter()
            .enumerate()
            .filter_map(|(slot, mesh)| mesh.as_ref().map(|mesh| (slot, mesh)))
    }

    /// Rebuild a document from persisted state.
    pub(crate) fn from_parts(
        metadata: ProjectMetadata,
        scene: SceneGraph,
        meshes: Vec<Option<IndexedMesh<f64>>>,
        dirty: bool,
    ) -> Self {
        Self {
            metadata,
            scene,
            meshes,
            dirty,
        }
    }
}

impl Default for ProjectDocument {
    fn default() -> Self {
        Self::new()
    }
}
