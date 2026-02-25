//! Region and boundary-label mapping for the blueprint pipeline.

use cfd_schematics::NodeKind;
use cfd_schematics::domain::therapy_metadata::TherapyZone;

use crate::domain::core::index::{FaceId, RegionId};
use crate::domain::mesh::IndexedMesh;
use crate::domain::topology::halfedge::PatchType;

/// Maps domain concepts to `RegionId` and boundary-patch labels.
pub struct RegionMap;

impl RegionMap {
    /// Map a `TherapyZone` variant to the corresponding `RegionId`.
    ///
    /// - `CancerTarget`  → `RegionId(10)`
    /// - `HealthyBypass` → `RegionId(11)`
    /// - `MixedFlow`     → `RegionId(12)`
    pub fn therapy_zone_to_region(zone: &TherapyZone) -> RegionId {
        match zone {
            TherapyZone::CancerTarget => RegionId::from_usize(10),
            TherapyZone::HealthyBypass => RegionId::from_usize(11),
            TherapyZone::MixedFlow => RegionId::from_usize(12),
        }
    }

    /// Map a `NodeKind` to a human-readable boundary label.
    ///
    /// - `Inlet`  → `"inlet"`
    /// - `Outlet` → `"outlet"`
    /// - others   → `"wall"`
    pub fn node_kind_to_label(kind: NodeKind) -> &'static str {
        match kind {
            NodeKind::Inlet => "inlet",
            NodeKind::Outlet => "outlet",
            _ => "wall",
        }
    }

    /// Map a `NodeKind` to a CFD boundary `PatchType`.
    ///
    /// - `Inlet`  → `PatchType::Inlet`
    /// - `Outlet` → `PatchType::Outlet`
    /// - others   → `PatchType::Wall`
    pub fn node_kind_to_patch_type(kind: NodeKind) -> PatchType {
        match kind {
            NodeKind::Inlet => PatchType::Inlet,
            NodeKind::Outlet => PatchType::Outlet,
            _ => PatchType::Wall,
        }
    }
}

/// Associates face IDs with boundary-patch labels for bulk application to a mesh.
#[derive(Debug, Default)]
pub struct BoundaryLabelMap {
    /// (face_id, label) pairs accumulated for later application.
    pub entries: Vec<(FaceId, String)>,
}

impl BoundaryLabelMap {
    /// Create an empty map.
    pub fn new() -> Self {
        Self::default()
    }

    /// Record a boundary label for a face.
    pub fn insert(&mut self, face_id: FaceId, label: impl Into<String>) {
        self.entries.push((face_id, label.into()));
    }

    /// Apply all recorded labels to `mesh` via `IndexedMesh::mark_boundary`.
    pub fn apply_to(&self, mesh: &mut IndexedMesh) {
        for (face_id, label) in &self.entries {
            mesh.mark_boundary(*face_id, label.clone());
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use cfd_schematics::domain::therapy_metadata::TherapyZone;

    #[test]
    fn cancer_target_maps_to_10() {
        assert_eq!(
            RegionMap::therapy_zone_to_region(&TherapyZone::CancerTarget),
            RegionId::from_usize(10)
        );
    }

    #[test]
    fn healthy_bypass_maps_to_11() {
        assert_eq!(
            RegionMap::therapy_zone_to_region(&TherapyZone::HealthyBypass),
            RegionId::from_usize(11)
        );
    }

    #[test]
    fn mixed_flow_maps_to_12() {
        assert_eq!(
            RegionMap::therapy_zone_to_region(&TherapyZone::MixedFlow),
            RegionId::from_usize(12)
        );
    }

    #[test]
    fn inlet_node_label() {
        assert_eq!(RegionMap::node_kind_to_label(NodeKind::Inlet), "inlet");
    }

    #[test]
    fn outlet_patch_type() {
        assert_eq!(
            RegionMap::node_kind_to_patch_type(NodeKind::Outlet),
            PatchType::Outlet
        );
    }

    #[test]
    fn apply_to_marks_faces() {
        let mut mesh = IndexedMesh::new();
        use crate::domain::core::scalar::{Point3r, Vector3r};
        let v0 = mesh.add_vertex(Point3r::new(0.0, 0.0, 0.0), Vector3r::new(0.0, 0.0, -1.0));
        let v1 = mesh.add_vertex(Point3r::new(1.0, 0.0, 0.0), Vector3r::new(0.0, 0.0, -1.0));
        let v2 = mesh.add_vertex(Point3r::new(0.0, 1.0, 0.0), Vector3r::new(0.0, 0.0, -1.0));
        let fid = mesh.add_face(v0, v1, v2);

        let v3 = mesh.add_vertex(Point3r::new(0.0, 0.0, 1.0), Vector3r::new(0.0, 0.0, 1.0));
        let v4 = mesh.add_vertex(Point3r::new(1.0, 0.0, 1.0), Vector3r::new(0.0, 0.0, 1.0));
        let v5 = mesh.add_vertex(Point3r::new(0.0, 1.0, 1.0), Vector3r::new(0.0, 0.0, 1.0));
        let fid2 = mesh.add_face(v3, v4, v5);

        let mut blm = BoundaryLabelMap::new();
        blm.insert(fid, "inlet");
        blm.insert(fid2, "outlet");
        blm.apply_to(&mut mesh);

        assert_eq!(mesh.boundary_label(fid), Some("inlet"));
        assert_eq!(mesh.boundary_label(fid2), Some("outlet"));
    }
}
