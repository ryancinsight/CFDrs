//! Mesh import commands — STL, OBJ, PLY, and auto-dispatch by extension.

use crate::domain::document::history::UndoableCommand;
use crate::domain::document::project::ProjectDocument;
use crate::domain::scene::graph::{MeshHandle, SceneEntity};
use std::path::PathBuf;

/// Command to import an STL file into the project.
pub struct ImportStlCommand {
    path: PathBuf,
    created_handle: Option<MeshHandle>,
    created_node: Option<usize>,
}

impl ImportStlCommand {
    /// Create a new STL import command.
    #[must_use]
    pub fn new(path: PathBuf) -> Self {
        Self {
            path,
            created_handle: None,
            created_node: None,
        }
    }
}

impl UndoableCommand for ImportStlCommand {
    fn execute(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        let file = std::fs::File::open(&self.path)?;
        let reader = std::io::BufReader::new(file);
        let mesh = cfd_mesh::infrastructure::io::stl::read_stl(reader)
            .map_err(|e| anyhow::anyhow!("{e}"))?;

        let name = self.path.file_stem().map_or_else(
            || "Imported STL".to_owned(),
            |s| s.to_string_lossy().into_owned(),
        );

        let handle = doc.add_mesh(mesh);
        let node_idx = doc.scene.add_node(name, SceneEntity::Mesh(handle));
        self.created_handle = Some(handle);
        self.created_node = Some(node_idx);
        Ok(())
    }

    fn undo(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        if let Some(node_idx) = self.created_node {
            doc.scene.remove_node(node_idx);
        }
        if let Some(handle) = self.created_handle {
            doc.remove_mesh(handle);
        }
        self.created_handle = None;
        self.created_node = None;
        Ok(())
    }

    fn description(&self) -> &'static str {
        "Import STL"
    }
}

/// Command to import any supported mesh file (STL, OBJ, PLY).
///
/// Detects the format from the file extension and dispatches to the
/// appropriate reader in `cfd_mesh::infrastructure::io`.
pub struct ImportMeshCommand {
    path: PathBuf,
    created_handle: Option<MeshHandle>,
    created_node: Option<usize>,
}

impl ImportMeshCommand {
    #[must_use]
    pub fn new(path: PathBuf) -> Self {
        Self {
            path,
            created_handle: None,
            created_node: None,
        }
    }
}

impl UndoableCommand for ImportMeshCommand {
    fn execute(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        let ext = self
            .path
            .extension()
            .and_then(|e| e.to_str())
            .unwrap_or("")
            .to_ascii_lowercase();

        let mesh = match ext.as_str() {
            "stl" => {
                let file = std::fs::File::open(&self.path)?;
                let reader = std::io::BufReader::new(file);
                cfd_mesh::infrastructure::io::stl::read_stl(reader)
                    .map_err(|e| anyhow::anyhow!("{e}"))?
            }
            "obj" => {
                let file = std::fs::File::open(&self.path)?;
                let reader = std::io::BufReader::new(file);
                cfd_mesh::infrastructure::io::obj::read_obj(reader)
                    .map_err(|e| anyhow::anyhow!("{e}"))?
            }
            "ply" => {
                let file = std::fs::File::open(&self.path)?;
                let reader = std::io::BufReader::new(file);
                cfd_mesh::infrastructure::io::ply::read_ply(reader)
                    .map_err(|e| anyhow::anyhow!("{e}"))?
            }
            _ => return Err(anyhow::anyhow!("unsupported file format: .{ext}")),
        };

        let name = self.path.file_stem().map_or_else(
            || "Imported Mesh".to_owned(),
            |s| s.to_string_lossy().into_owned(),
        );

        let handle = doc.add_mesh(mesh);
        let node_idx = doc.scene.add_node(name, SceneEntity::Mesh(handle));
        self.created_handle = Some(handle);
        self.created_node = Some(node_idx);
        Ok(())
    }

    fn undo(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        if let Some(node_idx) = self.created_node {
            doc.scene.remove_node(node_idx);
        }
        if let Some(handle) = self.created_handle {
            doc.remove_mesh(handle);
        }
        self.created_handle = None;
        self.created_node = None;
        Ok(())
    }

    fn description(&self) -> &'static str {
        "Import Mesh"
    }
}
