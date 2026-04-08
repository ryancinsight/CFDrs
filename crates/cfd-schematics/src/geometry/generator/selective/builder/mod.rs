//! Channel-network assembler for the four selective-tree topology variants.

mod cct;
mod cif;
mod dtcv;
mod tbt;

use std::collections::HashMap;

use super::super::super::types::Point2D;
use super::path_geometry::{build_serpentine_lobe_path, infer_serpentine_shape, polyline_length_mm};
use super::{CenterSerpentinePathSpec, SelectiveTreeRequest};
use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::geometry::metadata::{
    ChannelVisualRole, GeometryAuthoringProvenance, JunctionFamily, VenturiGeometryMetadata,
};

pub(super) struct SelectiveTreeBuilder {
    name: String,
    box_dims: (f64, f64),
    nodes: Vec<NodeSpec>,
    channels: Vec<ChannelSpec>,
    node_ids: HashMap<String, crate::domain::model::NodeId>,
}

impl SelectiveTreeBuilder {
    pub(super) fn new(name: String, box_dims: (f64, f64)) -> Self {
        Self {
            name,
            box_dims,
            nodes: Vec::new(),
            channels: Vec::new(),
            node_ids: HashMap::new(),
        }
    }

    pub(super) fn finish(self) -> NetworkBlueprint {
        let (w, h) = self.box_dims;
        let mut blueprint = NetworkBlueprint {
            name: self.name,
            box_dims: self.box_dims,
            nodes: self.nodes,
            channels: self.channels,
            render_hints: None,
            box_outline: vec![
                ((0.0, 0.0), (w, 0.0)),
                ((w, 0.0), (w, h)),
                ((w, h), (0.0, h)),
                ((0.0, h), (0.0, 0.0)),
            ],
            topology: None,
            lineage: None,
            metadata: None,
            geometry_authored: true,
        };
        blueprint.insert_metadata(GeometryAuthoringProvenance::selective_wrapper());
        blueprint
    }

    fn add_node(
        &mut self,
        name: &str,
        point: Point2D,
        kind: NodeKind,
        family: Option<JunctionFamily>,
    ) {
        let mut node = NodeSpec::new_at(name, kind, point);
        if let Some(junction_family) = family {
            node.junction_geometry = Some(crate::geometry::metadata::JunctionGeometryMetadata {
                junction_family,
                branch_angles_deg: Vec::new(),
                merge_angles_deg: Vec::new(),
            });
        }
        let id_str = node.id.clone();
        self.nodes.push(node);
        self.node_ids.insert(name.to_string(), id_str);
    }

    fn add_channel(
        &mut self,
        name: &str,
        from: &str,
        to: &str,
        path: Vec<Point2D>,
        width_m: f64,
        length_m: f64,
        height_m: f64,
        role: ChannelVisualRole,
        zone: crate::domain::therapy_metadata::TherapyZone,
        shape: Option<ChannelShape>,
        venturi: Option<VenturiGeometryMetadata>,
    ) {
        let final_length_m = if length_m > 0.0 {
            length_m
        } else if path.len() >= 2 {
            polyline_length_mm(&path) * 1.0e-3
        } else {
            length_m
        };
        let final_shape = match shape {
            Some(ChannelShape::Serpentine {
                segments,
                bend_radius_m,
                wave_type,
            }) => {
                if bend_radius_m == 0.0 {
                    path.first()
                        .zip(path.last())
                        .map_or(ChannelShape::Straight, |(start, end)| {
                            infer_serpentine_shape(&path, *start, *end, width_m * 1.0e3)
                        })
                } else {
                    ChannelShape::Serpentine {
                        segments,
                        bend_radius_m,
                        wave_type,
                    }
                }
            }
            Some(other) => other,
            None => ChannelShape::Straight,
        };

        let mut ch = ChannelSpec::new_pipe_rect(
            name,
            self.node_ids[from].0.clone(),
            self.node_ids[to].0.clone(),
            final_length_m,
            width_m,
            height_m,
            0.0,
            0.0,
        );
        ch.path = path;
        ch.channel_shape = final_shape;
        ch.visual_role = Some(role);
        ch.therapy_zone = Some(zone);
        ch.venturi_geometry = venturi;
        self.channels.push(ch);
    }

    fn y_mid(&self) -> f64 {
        self.box_dims.1 * 0.5
    }

    fn leaf_triplet(&self, center_y: f64, gap: f64) -> [f64; 3] {
        [center_y - gap, center_y, center_y + gap]
    }

    fn serpentine_path(
        &self,
        start: Point2D,
        end: Point2D,
        width_m: f64,
        spec: CenterSerpentinePathSpec,
        amplitude_mm: f64,
    ) -> Vec<Point2D> {
        build_serpentine_lobe_path(start, end, width_m, spec, amplitude_mm)
    }
}