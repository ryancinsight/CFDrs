//! Project document persistence for `cfd-ui`.
//!
//! # Theorem — Project Roundtrip Preservation
//!
//! For any saved project whose scene graph references only persisted mesh slots
//! and group nodes, `load_project_document(save_project_document(doc))`
//! reconstructs an equivalent document with identical metadata, camera state,
//! mesh geometry, boundary labels, scalar face attributes, and scene-node
//! hierarchy.
//!
//! **Proof sketch**: the snapshot format stores every persisted degree of
//! freedom explicitly: mesh vertices, normals, faces, region tags, boundary
//! labels, volumetric cells, scalar channels, scene-node transforms, parent /
//! child indices, and orbital-camera parameters. Loading rebuilds the exact
//! indexed-mesh and scene-graph structures without geometric welding, so all
//! identifiers and connectivity are preserved. The loader rejects unsupported
//! scene entities and broken parent/child relationships, preventing silent
//! state loss. Validation is exercised by
//! `project_roundtrip_preserves_scene_and_mesh_state`.

use std::fs;
use std::path::Path;

use anyhow::{Context, Result};
use cfd_mesh::domain::core::{FaceId, RegionId, VertexId};
use cfd_mesh::domain::topology::{Cell, ElementType};
use cfd_mesh::IndexedMesh;
use nalgebra::{Isometry3, Point3, Quaternion, Translation3, UnitQuaternion, Vector3};
use serde::{Deserialize, Serialize};

use crate::domain::document::project::{ProjectDocument, ProjectMetadata};
use crate::domain::scene::camera::OrbitalCamera;
use crate::domain::scene::graph::{MeshHandle, SceneEntity, SceneGraph, SceneNode};

#[derive(Debug, Clone, Serialize, Deserialize)]
struct ProjectSnapshot {
    metadata: ProjectMetadata,
    scene: SceneGraphSnapshot,
    meshes: Vec<Option<MeshSnapshot>>,
    dirty: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SceneGraphSnapshot {
    root: usize,
    camera: OrbitalCameraSnapshot,
    nodes: Vec<SceneNodeSnapshot>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct OrbitalCameraSnapshot {
    target: [f64; 3],
    distance: f64,
    azimuth: f64,
    elevation: f64,
    fov_y: f64,
    near: f64,
    far: f64,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct SceneNodeSnapshot {
    name: String,
    entity: SceneEntitySnapshot,
    translation: [f64; 3],
    rotation_xyzw: [f64; 4],
    visible: bool,
    children: Vec<usize>,
    parent: Option<usize>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
#[serde(tag = "kind", rename_all = "snake_case")]
enum SceneEntitySnapshot {
    Mesh { slot: usize },
    Group,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct MeshSnapshot {
    cell_size: f64,
    tolerance: Option<f64>,
    vertices: Vec<VertexSnapshot>,
    faces: Vec<FaceSnapshot>,
    cells: Vec<CellSnapshot>,
    boundary_labels: Vec<BoundaryLabelSnapshot>,
    face_attributes: Vec<FaceAttributeChannelSnapshot>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct VertexSnapshot {
    position: [f64; 3],
    normal: [f64; 3],
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FaceSnapshot {
    vertices: [u32; 3],
    region: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct CellSnapshot {
    faces: Vec<usize>,
    element_type: CellElementTypeSnapshot,
    vertex_ids: Vec<usize>,
}

#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
#[serde(rename_all = "snake_case")]
enum CellElementTypeSnapshot {
    Tetrahedron,
    Hexahedron,
    Triangle,
    Quadrilateral,
    Wedge,
    Pyramid,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct BoundaryLabelSnapshot {
    face: u32,
    label: String,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FaceAttributeChannelSnapshot {
    name: String,
    entries: Vec<FaceAttributeEntrySnapshot>,
}

#[derive(Debug, Clone, Serialize, Deserialize)]
struct FaceAttributeEntrySnapshot {
    face: u32,
    value: f64,
}

/// Save a project document to a deterministic JSON snapshot.
pub fn save_project_document(path: &Path, document: &mut ProjectDocument) -> Result<()> {
    let snapshot = ProjectSnapshot::from_document(document)?;
    if let Some(parent) = path.parent() {
        fs::create_dir_all(parent)
            .with_context(|| format!("failed to create project directory {}", parent.display()))?;
    }
    let json = serde_json::to_string_pretty(&snapshot)
        .context("failed to serialize project snapshot to JSON")?;
    fs::write(path, json)
        .with_context(|| format!("failed to write project file {}", path.display()))?;
    document.mark_saved();
    Ok(())
}

/// Load a project document from a JSON snapshot.
pub fn load_project_document(path: &Path) -> Result<ProjectDocument> {
    let json = fs::read_to_string(path)
        .with_context(|| format!("failed to read project file {}", path.display()))?;
    let snapshot: ProjectSnapshot = serde_json::from_str(&json)
        .with_context(|| format!("failed to parse project JSON {}", path.display()))?;
    snapshot.into_document()
}

impl ProjectSnapshot {
    fn from_document(document: &ProjectDocument) -> Result<Self> {
        let scene = SceneGraphSnapshot::from_scene(document)?;
        let mut meshes = Vec::new();
        let live_meshes = document.mesh_slots().collect::<Vec<_>>();
        let slot_count = live_meshes
            .iter()
            .map(|(slot, _)| slot + 1)
            .max()
            .unwrap_or(0);
        meshes.resize(slot_count, None);
        for (slot, mesh) in live_meshes {
            meshes[slot] = Some(MeshSnapshot::from_mesh(mesh));
        }
        Ok(Self {
            metadata: document.metadata.clone(),
            scene,
            meshes,
            dirty: document.is_dirty(),
        })
    }

    fn into_document(self) -> Result<ProjectDocument> {
        let meshes = self
            .meshes
            .into_iter()
            .map(|mesh| mesh.map(MeshSnapshot::into_mesh))
            .collect::<Vec<_>>();
        let scene = self.scene.into_scene(&meshes)?;
        Ok(ProjectDocument::from_parts(
            self.metadata,
            scene,
            meshes,
            self.dirty,
        ))
    }
}

impl SceneGraphSnapshot {
    fn from_scene(document: &ProjectDocument) -> Result<Self> {
        let scene = &document.scene;
        let nodes = scene
            .nodes()
            .iter()
            .map(SceneNodeSnapshot::from_node)
            .collect::<Result<Vec<_>>>()?;
        Ok(Self {
            root: scene.root(),
            camera: OrbitalCameraSnapshot::from_camera(scene.camera()),
            nodes,
        })
    }

    fn into_scene(self, meshes: &[Option<IndexedMesh<f64>>]) -> Result<SceneGraph> {
        let nodes = self
            .nodes
            .into_iter()
            .map(|node| node.into_node(meshes))
            .collect::<Result<Vec<_>>>()?;
        SceneGraph::from_parts(nodes, self.root, self.camera.into_camera())
    }
}

impl OrbitalCameraSnapshot {
    fn from_camera(camera: &OrbitalCamera) -> Self {
        Self {
            target: [camera.target.x, camera.target.y, camera.target.z],
            distance: camera.distance,
            azimuth: camera.azimuth,
            elevation: camera.elevation,
            fov_y: camera.fov_y,
            near: camera.near,
            far: camera.far,
        }
    }

    fn into_camera(self) -> OrbitalCamera {
        OrbitalCamera {
            target: Point3::new(self.target[0], self.target[1], self.target[2]),
            distance: self.distance,
            azimuth: self.azimuth,
            elevation: self.elevation,
            fov_y: self.fov_y,
            near: self.near,
            far: self.far,
        }
    }
}

impl SceneNodeSnapshot {
    fn from_node(node: &SceneNode) -> Result<Self> {
        let translation = node.transform.translation.vector;
        let quaternion = node.transform.rotation.quaternion();
        Ok(Self {
            name: node.name.clone(),
            entity: SceneEntitySnapshot::from_entity(&node.entity)?,
            translation: [translation.x, translation.y, translation.z],
            rotation_xyzw: [quaternion.i, quaternion.j, quaternion.k, quaternion.w],
            visible: node.visible,
            children: node.children.clone(),
            parent: node.parent,
        })
    }

    fn into_node(self, meshes: &[Option<IndexedMesh<f64>>]) -> Result<SceneNode> {
        let rotation = UnitQuaternion::from_quaternion(Quaternion::new(
            self.rotation_xyzw[3],
            self.rotation_xyzw[0],
            self.rotation_xyzw[1],
            self.rotation_xyzw[2],
        ));
        let translation = Translation3::new(
            self.translation[0],
            self.translation[1],
            self.translation[2],
        );
        Ok(SceneNode {
            name: self.name,
            entity: self.entity.into_entity(meshes)?,
            transform: Isometry3::from_parts(translation, rotation),
            visible: self.visible,
            children: self.children,
            parent: self.parent,
        })
    }
}

impl SceneEntitySnapshot {
    fn from_entity(entity: &SceneEntity) -> Result<Self> {
        match entity {
            SceneEntity::Mesh(handle) => Ok(Self::Mesh { slot: handle.0 }),
            SceneEntity::Group => Ok(Self::Group),
            SceneEntity::Schematic(_) => {
                anyhow::bail!("project persistence does not support schematic scene nodes yet")
            }
            SceneEntity::SimulationResult(_) => anyhow::bail!(
                "project persistence does not support simulation-result scene nodes yet"
            ),
        }
    }

    fn into_entity(self, meshes: &[Option<IndexedMesh<f64>>]) -> Result<SceneEntity> {
        match self {
            Self::Mesh { slot } => {
                let Some(Some(_)) = meshes.get(slot) else {
                    anyhow::bail!("scene references missing mesh slot {slot}");
                };
                Ok(SceneEntity::Mesh(MeshHandle(slot)))
            }
            Self::Group => Ok(SceneEntity::Group),
        }
    }
}

impl MeshSnapshot {
    fn from_mesh(mesh: &IndexedMesh<f64>) -> Self {
        let mut boundary_labels = mesh
            .boundary_labels
            .iter()
            .map(|(face, label)| BoundaryLabelSnapshot {
                face: face.0,
                label: label.clone(),
            })
            .collect::<Vec<_>>();
        boundary_labels.sort_by(|a, b| a.face.cmp(&b.face).then_with(|| a.label.cmp(&b.label)));

        let mut face_attributes = mesh
            .attributes
            .channel_names()
            .into_iter()
            .map(|name| {
                let mut entries = mesh
                    .attributes
                    .iter_channel(name)
                    .expect("channel listed by name must exist")
                    .map(|(face, value)| FaceAttributeEntrySnapshot {
                        face: face.0,
                        value,
                    })
                    .collect::<Vec<_>>();
                entries.sort_by_key(|entry| entry.face);
                FaceAttributeChannelSnapshot {
                    name: name.to_owned(),
                    entries,
                }
            })
            .collect::<Vec<_>>();
        face_attributes.sort_by(|a, b| a.name.cmp(&b.name));

        Self {
            cell_size: mesh.vertices.cell_size(),
            tolerance: mesh.vertices.tolerance(),
            vertices: mesh
                .vertices
                .iter()
                .map(|(_, vertex)| VertexSnapshot {
                    position: [vertex.position.x, vertex.position.y, vertex.position.z],
                    normal: [vertex.normal.x, vertex.normal.y, vertex.normal.z],
                })
                .collect(),
            faces: mesh
                .faces
                .iter()
                .map(|face| FaceSnapshot {
                    vertices: [face.vertices[0].0, face.vertices[1].0, face.vertices[2].0],
                    region: face.region.0,
                })
                .collect(),
            cells: mesh.cells().iter().map(CellSnapshot::from_cell).collect(),
            boundary_labels,
            face_attributes,
        }
    }

    fn into_mesh(self) -> IndexedMesh<f64> {
        let mut mesh = if let Some(tolerance) = self.tolerance {
            IndexedMesh::with_tolerance(self.cell_size, tolerance)
        } else {
            IndexedMesh::with_cell_size(self.cell_size)
        };

        for vertex in self.vertices {
            mesh.add_vertex_unique(
                Point3::new(vertex.position[0], vertex.position[1], vertex.position[2]),
                Vector3::new(vertex.normal[0], vertex.normal[1], vertex.normal[2]),
            );
        }

        for face in self.faces {
            mesh.add_face_with_region(
                VertexId::new(face.vertices[0]),
                VertexId::new(face.vertices[1]),
                VertexId::new(face.vertices[2]),
                RegionId::new(face.region),
            );
        }

        for cell in self.cells {
            mesh.add_cell(cell.into_cell());
        }

        for boundary in self.boundary_labels {
            mesh.mark_boundary(FaceId(boundary.face), boundary.label);
        }

        for channel in self.face_attributes {
            for entry in channel.entries {
                mesh.attributes
                    .set(&channel.name, FaceId(entry.face), entry.value);
            }
        }

        mesh.rebuild_edges();
        mesh
    }
}

impl CellSnapshot {
    fn from_cell(cell: &Cell) -> Self {
        Self {
            faces: cell.faces.clone(),
            element_type: CellElementTypeSnapshot::from(cell.element_type),
            vertex_ids: cell.vertex_ids.clone(),
        }
    }

    fn into_cell(self) -> Cell {
        Cell {
            faces: self.faces,
            element_type: self.element_type.into(),
            vertex_ids: self.vertex_ids,
        }
    }
}

impl From<ElementType> for CellElementTypeSnapshot {
    fn from(value: ElementType) -> Self {
        match value {
            ElementType::Tetrahedron => Self::Tetrahedron,
            ElementType::Hexahedron => Self::Hexahedron,
            ElementType::Triangle => Self::Triangle,
            ElementType::Quadrilateral => Self::Quadrilateral,
            ElementType::Wedge => Self::Wedge,
            ElementType::Pyramid => Self::Pyramid,
        }
    }
}

impl From<CellElementTypeSnapshot> for ElementType {
    fn from(value: CellElementTypeSnapshot) -> Self {
        match value {
            CellElementTypeSnapshot::Tetrahedron => Self::Tetrahedron,
            CellElementTypeSnapshot::Hexahedron => Self::Hexahedron,
            CellElementTypeSnapshot::Triangle => Self::Triangle,
            CellElementTypeSnapshot::Quadrilateral => Self::Quadrilateral,
            CellElementTypeSnapshot::Wedge => Self::Wedge,
            CellElementTypeSnapshot::Pyramid => Self::Pyramid,
        }
    }
}

#[cfg(test)]
mod tests {
    use std::time::{SystemTime, UNIX_EPOCH};

    use super::*;
    use crate::domain::scene::graph::SceneEntity;

    fn temp_project_path(label: &str) -> std::path::PathBuf {
        let nonce = SystemTime::now()
            .duration_since(UNIX_EPOCH)
            .expect("clock before unix epoch")
            .as_nanos();
        std::env::temp_dir().join(format!("cfd-ui-{label}-{nonce}.cfdproj"))
    }

    #[test]
    fn project_roundtrip_preserves_scene_and_mesh_state() {
        let path = temp_project_path("roundtrip");

        let mut document = ProjectDocument::new();
        document.metadata.name = "Roundtrip".to_owned();
        document.metadata.author = "CFDrs".to_owned();

        let mut mesh = IndexedMesh::with_tolerance(1e-3, 1e-3);
        let v0 = mesh.add_vertex_unique(Point3::new(0.0, 0.0, 0.0), Vector3::z());
        let v1 = mesh.add_vertex_unique(Point3::new(1.0, 0.0, 0.0), Vector3::z());
        let v2 = mesh.add_vertex_unique(Point3::new(0.0, 1.0, 0.0), Vector3::z());
        let face = mesh.add_face_with_region(v0, v1, v2, RegionId::new(7));
        mesh.mark_boundary(face, "inlet");
        mesh.attributes.set("quality", face, 0.875);
        mesh.rebuild_edges();

        let handle = document.add_mesh(mesh);
        let node = document
            .scene
            .add_node("triangle".to_owned(), SceneEntity::Mesh(handle));
        document.scene.camera_mut().target = Point3::new(1.0, 2.0, 3.0);
        document.scene.camera_mut().distance = 9.0;
        document
            .scene
            .node_mut(node)
            .expect("node exists")
            .transform = Isometry3::translation(1.5, -2.0, 0.25);

        save_project_document(&path, &mut document).expect("save succeeds");
        let loaded = load_project_document(&path).expect("load succeeds");

        assert_eq!(loaded.metadata.name, "Roundtrip");
        assert_eq!(loaded.metadata.author, "CFDrs");
        assert_eq!(loaded.scene.node_count(), 2);
        assert_eq!(loaded.scene.camera().target, Point3::new(1.0, 2.0, 3.0));
        assert_eq!(loaded.scene.camera().distance, 9.0);
        assert_eq!(loaded.mesh_count(), 1);

        let loaded_mesh = loaded.mesh(MeshHandle(0)).expect("mesh slot preserved");
        assert_eq!(loaded_mesh.vertex_count(), 3);
        assert_eq!(loaded_mesh.face_count(), 1);
        assert_eq!(loaded_mesh.boundary_label(FaceId(0)), Some("inlet"));
        assert_eq!(
            loaded_mesh.attributes.get("quality", FaceId(0)),
            Some(0.875)
        );

        let _ = fs::remove_file(path);
    }

    #[test]
    fn save_rejects_unpersisted_scene_entities() {
        let mut document = ProjectDocument::new();
        document.scene.add_node(
            "schematic".to_owned(),
            SceneEntity::Schematic(crate::domain::scene::graph::SchematicHandle(0)),
        );

        let path = temp_project_path("unsupported");
        let error = save_project_document(&path, &mut document)
            .err()
            .expect("unsupported entity must fail");
        assert!(error.to_string().contains("schematic scene nodes"));
    }
}
