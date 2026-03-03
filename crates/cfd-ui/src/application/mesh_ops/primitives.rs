//! Primitive mesh creation commands.

use crate::domain::document::history::UndoableCommand;
use crate::domain::document::project::ProjectDocument;
use crate::domain::scene::graph::{MeshHandle, SceneEntity};
use cfd_mesh::IndexedMesh;
use cfd_mesh::domain::geometry::primitives::PrimitiveError;

/// Specification for creating a primitive mesh.
#[derive(Clone, Debug)]
pub enum PrimitiveSpec {
    Cube { width: f64, height: f64, depth: f64 },
    Cylinder { radius: f64, height: f64, segments: usize },
    Sphere { radius: f64, segments: usize, stacks: usize },
    Cone { radius: f64, height: f64, segments: usize },
    Torus { major_radius: f64, minor_radius: f64, major_segments: usize, minor_segments: usize },
    Pipe { outer_radius: f64, inner_radius: f64, height: f64, segments: usize },
    Ellipsoid { semi_x: f64, semi_y: f64, semi_z: f64, segments: usize, stacks: usize },
    Capsule { radius: f64, cylinder_height: f64, segments: usize, hemisphere_stacks: usize },
    Frustum { bottom_radius: f64, top_radius: f64, height: f64, segments: usize },
    GeodesicSphere { radius: f64, frequency: usize },
    RoundedCube { width: f64, height: f64, depth: f64, corner_radius: f64, corner_segments: usize },
    Elbow { tube_radius: f64, bend_radius: f64, bend_angle: f64, tube_segments: usize, arc_segments: usize },
}

impl PrimitiveSpec {
    /// Build the primitive mesh from this specification.
    pub fn build(&self) -> Result<IndexedMesh<f64>, PrimitiveError> {
        use cfd_mesh::domain::geometry::primitives::*;
        use nalgebra::Point3;
        let origin = Point3::<f64>::origin();
        match self {
            Self::Cube { width, height, depth } => {
                (Cube { origin, width: *width, height: *height, depth: *depth }).build()
            }
            Self::Cylinder { radius, height, segments } => {
                (Cylinder { base_center: origin, radius: *radius, height: *height, segments: *segments }).build()
            }
            Self::Sphere { radius, segments, stacks } => {
                (UvSphere { center: origin, radius: *radius, segments: *segments, stacks: *stacks }).build()
            }
            Self::Cone { radius, height, segments } => {
                (Cone { base_center: origin, radius: *radius, height: *height, segments: *segments }).build()
            }
            Self::Torus { major_radius, minor_radius, major_segments, minor_segments } => {
                (Torus { major_radius: *major_radius, minor_radius: *minor_radius, major_segments: *major_segments, minor_segments: *minor_segments }).build()
            }
            Self::Pipe { outer_radius, inner_radius, height, segments } => {
                (Pipe { base_center: origin, outer_radius: *outer_radius, inner_radius: *inner_radius, height: *height, segments: *segments }).build()
            }
            Self::Ellipsoid { semi_x, semi_y, semi_z, segments, stacks } => {
                (Ellipsoid { center: origin, semi_x: *semi_x, semi_y: *semi_y, semi_z: *semi_z, segments: *segments, stacks: *stacks }).build()
            }
            Self::Capsule { radius, cylinder_height, segments, hemisphere_stacks } => {
                (Capsule { center: origin, radius: *radius, cylinder_height: *cylinder_height, segments: *segments, hemisphere_stacks: *hemisphere_stacks }).build()
            }
            Self::Frustum { bottom_radius, top_radius, height, segments } => {
                (Frustum { base_center: origin, bottom_radius: *bottom_radius, top_radius: *top_radius, height: *height, segments: *segments }).build()
            }
            Self::GeodesicSphere { radius, frequency } => {
                (GeodesicSphere { center: origin, radius: *radius, frequency: *frequency }).build()
            }
            Self::RoundedCube { width, height, depth, corner_radius, corner_segments } => {
                (RoundedCube { origin, width: *width, height: *height, depth: *depth, corner_radius: *corner_radius, corner_segments: *corner_segments }).build()
            }
            Self::Elbow { tube_radius, bend_radius, bend_angle, tube_segments, arc_segments } => {
                (Elbow { tube_radius: *tube_radius, bend_radius: *bend_radius, bend_angle: *bend_angle, tube_segments: *tube_segments, arc_segments: *arc_segments }).build()
            }
        }
    }

    /// Human-readable name for this primitive type.
    #[must_use]
    pub fn type_name(&self) -> &str {
        match self {
            Self::Cube { .. } => "Cube",
            Self::Cylinder { .. } => "Cylinder",
            Self::Sphere { .. } => "Sphere",
            Self::Cone { .. } => "Cone",
            Self::Torus { .. } => "Torus",
            Self::Pipe { .. } => "Pipe",
            Self::Ellipsoid { .. } => "Ellipsoid",
            Self::Capsule { .. } => "Capsule",
            Self::Frustum { .. } => "Frustum",
            Self::GeodesicSphere { .. } => "Geodesic Sphere",
            Self::RoundedCube { .. } => "Rounded Cube",
            Self::Elbow { .. } => "Elbow",
        }
    }
}

/// Command to create a primitive mesh and add it to the scene.
pub struct CreatePrimitiveCommand {
    spec: PrimitiveSpec,
    name: String,
    created_handle: Option<MeshHandle>,
    created_node: Option<usize>,
}

impl CreatePrimitiveCommand {
    /// Create a new primitive creation command.
    #[must_use]
    pub fn new(spec: PrimitiveSpec, name: String) -> Self {
        Self {
            spec,
            name,
            created_handle: None,
            created_node: None,
        }
    }
}

impl UndoableCommand for CreatePrimitiveCommand {
    fn execute(&mut self, doc: &mut ProjectDocument) -> anyhow::Result<()> {
        let mesh = self.spec.build().map_err(|e| anyhow::anyhow!("{e}"))?;
        let handle = doc.add_mesh(mesh);
        let node_idx = doc.scene.add_node(
            self.name.clone(),
            SceneEntity::Mesh(handle),
        );
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

    fn description(&self) -> &str {
        "Create Primitive"
    }
}
