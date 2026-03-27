use super::super::types::{Point2D, SplitType};
use crate::config::{ChannelTypeConfig, GeometryConfig};
use crate::domain::model::{ChannelShape, ChannelSpec, NetworkBlueprint, NodeKind, NodeSpec};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::{
    BlueprintRenderHints, ChannelVenturiSpec, ChannelVisualRole, GeometryAuthoringProvenance,
    IncrementalFiltrationParams, JunctionFamily, JunctionGeometryMetadata, MetadataContainer,
    VenturiGeometryMetadata,
};
use crate::topology::{
    BlueprintTopologyFactory, BlueprintTopologySpec, SplitKind, TreatmentActuationMode,
};
use std::cmp::Ordering;
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CenterSerpentinePathSpec {
    pub segments: usize,
    pub bend_radius_m: f64,
    /// Waveform type for the serpentine path.  Defaults to Sine when
    /// constructed from legacy topology specs that omit this field.
    pub wave_type: crate::topology::SerpentineWaveType,
}

use crate::geometry::builders::ChannelExt;

#[derive(Debug, Clone, Copy)]
struct PendingVenturiPath {
    channel_idx: usize,
    start: Option<Point2D>,
    end: Option<Point2D>,
    preferred_y: f64,
    fallback_length_m: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SelectiveTreeTopology {
    CascadeCenterTrifurcation {
        n_levels: usize,
        center_frac: f64,
        venturi_treatment_enabled: bool,
        center_serpentine: Option<super::selective::CenterSerpentinePathSpec>,
    },
    IncrementalFiltrationTriBi {
        n_pretri: usize,
        pretri_center_frac: f64,
        terminal_tri_center_frac: f64,
        bi_treat_frac: f64,
        venturi_treatment_enabled: bool,
        center_serpentine: Option<super::selective::CenterSerpentinePathSpec>,
        outlet_tail_length_m: f64,
    },
    TriBiTriSelective {
        first_center_frac: f64,
        bi_treat_frac: f64,
        second_center_frac: f64,
    },
    DoubleTrifurcationCif {
        split1_center_frac: f64,
        split2_center_frac: f64,
        center_throat_count: u8,
        inter_throat_spacing_m: f64,
    },
}

#[derive(Debug, Clone, PartialEq)]
pub struct SelectiveTreeRequest {
    pub name: String,
    pub box_dims_mm: (f64, f64),
    pub trunk_length_m: f64,
    pub branch_length_m: f64,
    pub hybrid_branch_length_m: f64,
    pub main_width_m: f64,
    pub throat_width_m: f64,
    pub throat_length_m: f64,
    pub channel_height_m: f64,
    pub topology: SelectiveTreeTopology,
}

pub fn create_selective_tree_geometry(request: &SelectiveTreeRequest) -> NetworkBlueprint {
    let mut builder = SelectiveTreeBuilder::new(request.name.clone(), request.box_dims_mm);
    match &request.topology {
        SelectiveTreeTopology::CascadeCenterTrifurcation {
            n_levels,
            center_frac,
            venturi_treatment_enabled,
            center_serpentine,
        } => builder.build_cct(
            *n_levels,
            *center_frac,
            *venturi_treatment_enabled,
            *center_serpentine,
            request,
        ),
        SelectiveTreeTopology::IncrementalFiltrationTriBi {
            n_pretri,
            pretri_center_frac,
            terminal_tri_center_frac,
            bi_treat_frac,
            venturi_treatment_enabled,
            center_serpentine,
            outlet_tail_length_m,
        } => builder.build_cif(
            *n_pretri,
            *pretri_center_frac,
            *terminal_tri_center_frac,
            *bi_treat_frac,
            *venturi_treatment_enabled,
            *center_serpentine,
            *outlet_tail_length_m,
            request,
        ),
        SelectiveTreeTopology::TriBiTriSelective {
            first_center_frac,
            bi_treat_frac,
            second_center_frac,
        } => builder.build_tbt(
            *first_center_frac,
            *bi_treat_frac,
            *second_center_frac,
            request,
        ),
        SelectiveTreeTopology::DoubleTrifurcationCif {
            split1_center_frac,
            split2_center_frac,
            center_throat_count,
            inter_throat_spacing_m,
        } => builder.build_dtcv(
            *split1_center_frac,
            *split2_center_frac,
            *center_throat_count,
            *inter_throat_spacing_m,
            request,
        ),
    }
    builder.finish()
}

struct SelectiveTreeBuilder {
    name: String,
    box_dims: (f64, f64),
    nodes: Vec<NodeSpec>,
    channels: Vec<ChannelSpec>,
    node_ids: HashMap<String, crate::domain::model::NodeId>,
}

impl SelectiveTreeBuilder {
    fn new(name: String, box_dims: (f64, f64)) -> Self {
        Self {
            name,
            box_dims,
            nodes: Vec::new(),
            channels: Vec::new(),
            node_ids: HashMap::new(),
        }
    }

    fn finish(self) -> NetworkBlueprint {
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

    fn build_cct(
        &mut self,
        n_levels: usize,
        center_frac: f64,
        venturi: bool,
        center_serp: Option<super::selective::CenterSerpentinePathSpec>,
        req: &SelectiveTreeRequest,
    ) {
        use crate::domain::therapy_metadata::TherapyZone;
        let y_mid = self.y_mid();
        self.add_node("inlet", (0.0, y_mid), NodeKind::Inlet, None);
        self.add_node(
            "periph_merge",
            (104.0, y_mid),
            NodeKind::Junction,
            Some(JunctionFamily::Merge),
        );
        self.add_node(
            "outlet_merge",
            (116.0, y_mid),
            NodeKind::Junction,
            Some(JunctionFamily::Merge),
        );
        self.add_node("throat_in", (68.0, y_mid), NodeKind::Junction, None);
        self.add_node("throat_out", (76.0, y_mid), NodeKind::Junction, None);
        self.add_node("outlet", (self.box_dims.0, y_mid), NodeKind::Outlet, None);
        for lv in 0..n_levels {
            self.add_node(
                &format!("split_lv{lv}"),
                (18.0 + lv as f64 * 12.0, y_mid),
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            );
        }
        self.add_channel(
            "inlet_section",
            "inlet",
            "split_lv0",
            vec![(0.0, y_mid), (18.0, y_mid)],
            req.main_width_m,
            0.0,
            req.channel_height_m,
            ChannelVisualRole::Trunk,
            TherapyZone::MixedFlow,
            None,
            None,
        );
        let mut center_width = req.main_width_m;
        let side_frac = (1.0 - center_frac) * 0.5;
        for lv in 0..n_levels {
            let split = format!("split_lv{lv}");
            let next = if lv + 1 < n_levels {
                format!("split_lv{}", lv + 1)
            } else {
                "throat_in".to_string()
            };
            let off = 12.0 + 5.0 * (n_levels.saturating_sub(lv + 1)) as f64;
            let side_w = (center_width * side_frac).max(req.throat_width_m);
            let ctr_w = (center_width * center_frac).max(req.throat_width_m);
            for (id, y) in [("L", y_mid - off), ("R", y_mid + off)] {
                let start_pt = self
                    .nodes
                    .iter()
                    .find(|n| n.id == self.node_ids[&split])
                    .unwrap()
                    .point;
                self.add_channel(
                    &format!("{id}_lv{lv}"),
                    &split,
                    "periph_merge",
                    vec![
                        start_pt,
                        (34.0 + lv as f64 * 10.0, y),
                        (90.0, y),
                        (104.0, y_mid),
                    ],
                    side_w,
                    0.0,
                    req.channel_height_m,
                    ChannelVisualRole::PeripheralBypass,
                    TherapyZone::HealthyBypass,
                    None,
                    None,
                );
            }
            let start = self
                .nodes
                .iter()
                .find(|n| n.id == self.node_ids[&split])
                .unwrap()
                .point;
            let end = self
                .nodes
                .iter()
                .find(|n| n.id == self.node_ids[&next])
                .unwrap()
                .point;
            let path = if let Some(spec) = center_serp {
                self.serpentine_path(start, end, ctr_w, spec, 1.6)
            } else {
                vec![start, end]
            };
            self.add_channel(
                &format!("center_lv{lv}"),
                &split,
                &next,
                path,
                ctr_w,
                0.0,
                req.channel_height_m,
                ChannelVisualRole::CenterTreatment,
                TherapyZone::CancerTarget,
                center_serp.map(|_| ChannelShape::Serpentine {
                    segments: 2,
                    bend_radius_m: 0.0,
                    wave_type: crate::topology::SerpentineWaveType::Sine,
                }),
                None,
            );
            center_width = ctr_w;
        }
        let zone = TherapyZone::CancerTarget;
        if venturi {
            self.add_channel(
                "throat_section",
                "throat_in",
                "throat_out",
                vec![(68.0, y_mid), (72.0, y_mid), (76.0, y_mid)],
                req.throat_width_m,
                req.throat_length_m,
                req.channel_height_m,
                ChannelVisualRole::VenturiThroat,
                zone,
                None,
                Some(VenturiGeometryMetadata {
                    throat_width_m: req.throat_width_m,
                    throat_height_m: req.channel_height_m,
                    throat_length_m: req.throat_length_m,
                    inlet_width_m: center_width,
                    outlet_width_m: center_width,
                    convergent_half_angle_deg: 15.0,
                    divergent_half_angle_deg: 15.0,
                    throat_position: 0.5,
                }),
            );
        } else {
            let path = if let Some(spec) = center_serp {
                self.serpentine_path((68.0, y_mid), (76.0, y_mid), center_width, spec, 1.5)
            } else {
                vec![(68.0, y_mid), (76.0, y_mid)]
            };
            self.add_channel(
                "treatment_section",
                "throat_in",
                "throat_out",
                path,
                center_width,
                if center_serp.is_some() {
                    0.0
                } else {
                    req.branch_length_m.max(req.throat_length_m)
                },
                req.channel_height_m,
                ChannelVisualRole::CenterTreatment,
                zone,
                center_serp.map(|_| ChannelShape::Serpentine {
                    segments: 2,
                    bend_radius_m: 0.0,
                    wave_type: crate::topology::SerpentineWaveType::Sine,
                }),
                None,
            );
        }
        self.add_channel(
            "center_to_merge",
            "throat_out",
            "outlet_merge",
            vec![(76.0, y_mid), (116.0, y_mid)],
            center_width,
            req.branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::MergeCollector,
            TherapyZone::MixedFlow,
            None,
            None,
        );
        self.add_channel(
            "periph_to_merge",
            "periph_merge",
            "outlet_merge",
            vec![(104.0, y_mid), (116.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::MergeCollector,
            TherapyZone::HealthyBypass,
            None,
            None,
        );
        self.add_channel(
            "trunk_out",
            "outlet_merge",
            "outlet",
            vec![(116.0, y_mid), (self.box_dims.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::Trunk,
            TherapyZone::MixedFlow,
            None,
            None,
        );
    }

    fn build_cif(
        &mut self,
        n_pretri: usize,
        pretri_center_frac: f64,
        terminal_tri_center_frac: f64,
        bi_treat_frac: f64,
        venturi: bool,
        center_serp: Option<super::selective::CenterSerpentinePathSpec>,
        outlet_tail_length_m: f64,
        req: &SelectiveTreeRequest,
    ) {
        use crate::domain::therapy_metadata::TherapyZone;
        let y_mid = self.y_mid();
        for (name, x, kind, family) in [
            ("inlet", 0.0, NodeKind::Inlet, None),
            (
                "periph_merge",
                104.0,
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            ),
            (
                "outlet_merge",
                116.0,
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            ),
            (
                "hy_tri",
                52.0,
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            ),
            (
                "hy_bi",
                64.0,
                NodeKind::Junction,
                Some(JunctionFamily::Bifurcation),
            ),
            ("throat_in", 72.0, NodeKind::Junction, None),
            ("throat_out", 80.0, NodeKind::Junction, None),
            ("outlet", self.box_dims.0, NodeKind::Outlet, None),
        ] {
            self.add_node(name, (x, y_mid), kind, family);
        }
        for lv in 0..n_pretri {
            self.add_node(
                &format!("split_lv{lv}"),
                (18.0 + lv as f64 * 11.0, y_mid),
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            );
        }
        self.add_channel(
            "inlet_section",
            "inlet",
            "split_lv0",
            vec![(0.0, y_mid), (18.0, y_mid)],
            req.main_width_m,
            0.0,
            req.channel_height_m,
            ChannelVisualRole::Trunk,
            TherapyZone::MixedFlow,
            None,
            None,
        );

        if let Some(ch) = self.channels.last_mut() {
            ch.add_metadata(IncrementalFiltrationParams {
                n_pretri: n_pretri as u8,
                pretri_center_frac,
                terminal_tri_center_frac,
                bi_treat_frac,
                outlet_tail_length_m,
            });
        }
        let mut center_width = req.main_width_m;
        for lv in 0..n_pretri {
            let split = format!("split_lv{lv}");
            let next = if lv + 1 < n_pretri {
                format!("split_lv{}", lv + 1)
            } else {
                "hy_tri".to_string()
            };
            let side_off = 10.0 + 4.0 * (n_pretri.saturating_sub(lv + 1)) as f64;
            let t = if n_pretri > 1 {
                lv as f64 / (n_pretri - 1) as f64
            } else {
                0.0
            };
            let level_frac =
                pretri_center_frac + t * (terminal_tri_center_frac - pretri_center_frac);
            let side_w = (center_width * (1.0 - level_frac) * 0.5).max(req.throat_width_m);
            let ctr_w = (center_width * level_frac).max(req.throat_width_m);
            for (id, y) in [("L", y_mid - side_off), ("R", y_mid + side_off)] {
                let start_pt = self
                    .nodes
                    .iter()
                    .find(|n| n.id == self.node_ids[&split])
                    .unwrap()
                    .point;
                self.add_channel(
                    &format!("{id}_lv{lv}"),
                    &split,
                    "periph_merge",
                    vec![
                        start_pt,
                        (34.0 + lv as f64 * 9.0, y),
                        (50.0 + lv as f64 * 8.0, y),
                        (92.0, y),
                        (104.0, y_mid),
                    ],
                    side_w,
                    req.branch_length_m,
                    req.channel_height_m,
                    ChannelVisualRole::PeripheralBypass,
                    TherapyZone::HealthyBypass,
                    None,
                    None,
                );
            }
            let start = self
                .nodes
                .iter()
                .find(|n| n.id == self.node_ids[&split])
                .unwrap()
                .point;
            let end = self
                .nodes
                .iter()
                .find(|n| n.id == self.node_ids[&next])
                .unwrap()
                .point;
            let path = if let Some(spec) = center_serp {
                self.serpentine_path(start, end, ctr_w, spec, 1.4)
            } else {
                vec![start, end]
            };
            self.add_channel(
                &format!("center_lv{lv}"),
                &split,
                &next,
                path,
                ctr_w,
                req.branch_length_m,
                req.channel_height_m,
                ChannelVisualRole::CenterTreatment,
                TherapyZone::CancerTarget,
                center_serp.map(|spec| ChannelShape::Serpentine {
                    segments: spec.segments,
                    bend_radius_m: spec.bend_radius_m,
                    wave_type: spec.wave_type,
                }),
                None,
            );
            center_width = ctr_w;
        }
        let tri_side =
            (center_width * (1.0 - terminal_tri_center_frac) * 0.5).max(req.throat_width_m);
        let tri_ctr = (center_width * terminal_tri_center_frac).max(req.throat_width_m);
        for (id, y) in [("hy_tri_L", y_mid - 8.0), ("hy_tri_R", y_mid + 8.0)] {
            self.add_channel(
                id,
                "hy_tri",
                "periph_merge",
                vec![(52.0, y_mid), (62.0, y), (92.0, y), (104.0, y_mid)],
                tri_side,
                req.hybrid_branch_length_m,
                req.channel_height_m,
                ChannelVisualRole::PeripheralBypass,
                TherapyZone::HealthyBypass,
                None,
                None,
            );
        }
        let tri_path = if let Some(spec) = center_serp {
            self.serpentine_path((52.0, y_mid), (64.0, y_mid), tri_ctr, spec, 1.3)
        } else {
            vec![(52.0, y_mid), (64.0, y_mid)]
        };
        self.add_channel(
            "hy_tri_center",
            "hy_tri",
            "hy_bi",
            tri_path,
            tri_ctr,
            req.hybrid_branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::CenterTreatment,
            TherapyZone::CancerTarget,
            center_serp.map(|spec| ChannelShape::Serpentine {
                segments: spec.segments,
                bend_radius_m: spec.bend_radius_m,
                wave_type: spec.wave_type,
            }),
            None,
        );
        let treat_w = (tri_ctr * bi_treat_frac).max(req.throat_width_m);
        let bypass_w = (tri_ctr * (1.0 - bi_treat_frac)).max(req.throat_width_m);
        self.add_channel(
            "hy_bi_bypass",
            "hy_bi",
            "periph_merge",
            vec![
                (64.0, y_mid),
                (76.0, y_mid + 10.0),
                (92.0, y_mid + 10.0),
                (104.0, y_mid),
            ],
            bypass_w,
            req.hybrid_branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::PeripheralBypass,
            TherapyZone::HealthyBypass,
            None,
            None,
        );
        let bi_path = if let Some(spec) = center_serp {
            self.serpentine_path((64.0, y_mid), (72.0, y_mid), treat_w, spec, 1.2)
        } else {
            vec![(64.0, y_mid), (72.0, y_mid)]
        };
        self.add_channel(
            "hy_bi_treat",
            "hy_bi",
            "throat_in",
            bi_path,
            treat_w,
            req.hybrid_branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::CenterTreatment,
            TherapyZone::CancerTarget,
            center_serp.map(|spec| ChannelShape::Serpentine {
                segments: spec.segments,
                bend_radius_m: spec.bend_radius_m,
                wave_type: spec.wave_type,
            }),
            None,
        );
        if venturi {
            self.add_channel(
                "throat_section",
                "throat_in",
                "throat_out",
                vec![(72.0, y_mid), (76.0, y_mid), (80.0, y_mid)],
                req.throat_width_m,
                req.throat_length_m,
                req.channel_height_m,
                ChannelVisualRole::VenturiThroat,
                TherapyZone::CancerTarget,
                None,
                Some(VenturiGeometryMetadata {
                    throat_width_m: req.throat_width_m,
                    throat_height_m: req.channel_height_m,
                    throat_length_m: req.throat_length_m,
                    inlet_width_m: treat_w,
                    outlet_width_m: treat_w,
                    convergent_half_angle_deg: 15.0,
                    divergent_half_angle_deg: 15.0,
                    throat_position: 0.5,
                }),
            );
        } else {
            let path = if let Some(spec) = center_serp {
                self.serpentine_path((72.0, y_mid), (80.0, y_mid), treat_w, spec, 1.2)
            } else {
                vec![(72.0, y_mid), (80.0, y_mid)]
            };
            self.add_channel(
                "treatment_section",
                "throat_in",
                "throat_out",
                path,
                treat_w,
                if center_serp.is_some() {
                    0.0
                } else {
                    req.hybrid_branch_length_m.max(req.throat_length_m)
                },
                req.channel_height_m,
                ChannelVisualRole::CenterTreatment,
                TherapyZone::CancerTarget,
                center_serp.map(|spec| ChannelShape::Serpentine {
                    segments: spec.segments,
                    bend_radius_m: spec.bend_radius_m,
                    wave_type: spec.wave_type,
                }),
                None,
            );
        }
        self.add_channel(
            "center_to_merge",
            "throat_out",
            "outlet_merge",
            vec![(80.0, y_mid), (116.0, y_mid)],
            treat_w,
            req.hybrid_branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::MergeCollector,
            TherapyZone::MixedFlow,
            None,
            None,
        );
        self.add_channel(
            "periph_to_merge",
            "periph_merge",
            "outlet_merge",
            vec![(104.0, y_mid), (116.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::MergeCollector,
            TherapyZone::HealthyBypass,
            None,
            None,
        );
        self.add_channel(
            "trunk_out",
            "outlet_merge",
            "outlet",
            vec![(116.0, y_mid), (self.box_dims.0, y_mid)],
            req.main_width_m,
            outlet_tail_length_m,
            req.channel_height_m,
            ChannelVisualRole::Trunk,
            TherapyZone::MixedFlow,
            None,
            None,
        );
    }

    fn build_tbt(
        &mut self,
        first_center_frac: f64,
        bi_treat_frac: f64,
        second_center_frac: f64,
        req: &SelectiveTreeRequest,
    ) {
        use crate::domain::therapy_metadata::TherapyZone;
        let y_mid = self.y_mid();
        for (name, x, kind, family) in [
            ("inlet", 0.0, NodeKind::Inlet, None),
            (
                "split_lv0",
                24.0,
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            ),
            (
                "split_bi",
                46.0,
                NodeKind::Junction,
                Some(JunctionFamily::Bifurcation),
            ),
            (
                "split_lv1",
                70.0,
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            ),
            ("throat_in", 86.0, NodeKind::Junction, None),
            ("throat_out", 92.0, NodeKind::Junction, None),
            (
                "periph_merge",
                104.0,
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            ),
            (
                "outlet_merge",
                116.0,
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            ),
            ("outlet", self.box_dims.0, NodeKind::Outlet, None),
        ] {
            self.add_node(name, (x, y_mid), kind, family);
        }
        self.add_channel(
            "inlet_section",
            "inlet",
            "split_lv0",
            vec![(0.0, y_mid), (24.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::Trunk,
            TherapyZone::MixedFlow,
            None,
            None,
        );
        let side1 = (req.main_width_m * (1.0 - first_center_frac) * 0.5).max(req.throat_width_m);
        let ctr1 = (req.main_width_m * first_center_frac).max(req.throat_width_m);
        for (id, y) in [("L_lv0", y_mid - 18.0), ("R_lv0", y_mid + 18.0)] {
            self.add_channel(
                id,
                "split_lv0",
                "periph_merge",
                vec![(24.0, y_mid), (40.0, y), (96.0, y), (104.0, y_mid)],
                side1,
                req.branch_length_m,
                req.channel_height_m,
                ChannelVisualRole::PeripheralBypass,
                TherapyZone::HealthyBypass,
                None,
                None,
            );
        }
        let bi_bypass_w = (ctr1 * (1.0 - bi_treat_frac)).max(req.throat_width_m);
        let bi_treat_w = (ctr1 * bi_treat_frac).max(req.throat_width_m);
        self.add_channel(
            "center_lv0",
            "split_lv0",
            "split_bi",
            vec![(24.0, y_mid), (46.0, y_mid)],
            ctr1,
            req.branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::CenterTreatment,
            TherapyZone::CancerTarget,
            None,
            None,
        );
        self.add_channel(
            "bi_bypass",
            "split_bi",
            "periph_merge",
            vec![
                (46.0, y_mid),
                (60.0, y_mid + 10.0),
                (96.0, y_mid + 10.0),
                (104.0, y_mid),
            ],
            bi_bypass_w,
            req.hybrid_branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::PeripheralBypass,
            TherapyZone::HealthyBypass,
            None,
            None,
        );
        self.add_channel(
            "center_bi",
            "split_bi",
            "split_lv1",
            vec![(46.0, y_mid), (70.0, y_mid)],
            bi_treat_w,
            req.hybrid_branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::CenterTreatment,
            TherapyZone::CancerTarget,
            None,
            None,
        );
        let side2 = (bi_treat_w * (1.0 - second_center_frac) * 0.5).max(req.throat_width_m);
        let ctr2 = (bi_treat_w * second_center_frac).max(req.throat_width_m);
        for (id, y) in [("L_lv1", y_mid - 10.0), ("R_lv1", y_mid + 10.0)] {
            self.add_channel(
                id,
                "split_lv1",
                "periph_merge",
                vec![(70.0, y_mid), (82.0, y), (96.0, y), (104.0, y_mid)],
                side2,
                req.hybrid_branch_length_m,
                req.channel_height_m,
                ChannelVisualRole::PeripheralBypass,
                TherapyZone::HealthyBypass,
                None,
                None,
            );
        }
        self.add_channel(
            "throat_section",
            "throat_in",
            "throat_out",
            vec![(86.0, y_mid), (89.0, y_mid), (92.0, y_mid)],
            req.throat_width_m,
            req.throat_length_m,
            req.channel_height_m,
            ChannelVisualRole::VenturiThroat,
            TherapyZone::CancerTarget,
            None,
            Some(VenturiGeometryMetadata {
                throat_width_m: req.throat_width_m,
                throat_height_m: req.channel_height_m,
                throat_length_m: req.throat_length_m,
                inlet_width_m: ctr2,
                outlet_width_m: ctr2,
                convergent_half_angle_deg: 15.0,
                divergent_half_angle_deg: 15.0,
                throat_position: 0.5,
            }),
        );
        self.add_channel(
            "center_lv1",
            "split_lv1",
            "throat_in",
            vec![(70.0, y_mid), (86.0, y_mid)],
            ctr2,
            req.hybrid_branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::CenterTreatment,
            TherapyZone::CancerTarget,
            None,
            None,
        );
        self.add_channel(
            "center_to_merge",
            "throat_out",
            "outlet_merge",
            vec![(92.0, y_mid), (116.0, y_mid)],
            ctr2,
            req.hybrid_branch_length_m,
            req.channel_height_m,
            ChannelVisualRole::MergeCollector,
            TherapyZone::MixedFlow,
            None,
            None,
        );
        self.add_channel(
            "periph_to_merge",
            "periph_merge",
            "outlet_merge",
            vec![(104.0, y_mid), (116.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::MergeCollector,
            TherapyZone::HealthyBypass,
            None,
            None,
        );
        self.add_channel(
            "trunk_out",
            "outlet_merge",
            "outlet",
            vec![(116.0, y_mid), (self.box_dims.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::Trunk,
            TherapyZone::MixedFlow,
            None,
            None,
        );
    }

    fn build_dtcv(
        &mut self,
        split1_center_frac: f64,
        split2_center_frac: f64,
        center_throat_count: u8,
        inter_spacing_m: f64,
        req: &SelectiveTreeRequest,
    ) {
        use crate::domain::therapy_metadata::TherapyZone;
        let y_mid = self.y_mid();
        for (name, x, y, kind, family) in [
            ("inlet", 0.0, y_mid, NodeKind::Inlet, None),
            (
                "split1",
                18.0,
                y_mid,
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            ),
            (
                "split2_upper",
                36.0,
                y_mid - 21.0,
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            ),
            (
                "split2",
                36.0,
                y_mid,
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            ),
            (
                "split2_lower",
                36.0,
                y_mid + 21.0,
                NodeKind::Junction,
                Some(JunctionFamily::Trifurcation),
            ),
            (
                "merge_center",
                100.0,
                y_mid,
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            ),
            (
                "periph_merge",
                104.0,
                y_mid,
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            ),
            (
                "merge_out",
                118.0,
                y_mid,
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            ),
            ("outlet", self.box_dims.0, y_mid, NodeKind::Outlet, None),
        ] {
            self.add_node(name, (x, y), kind, family);
        }
        self.add_channel(
            "inlet_section",
            "inlet",
            "split1",
            vec![(0.0, y_mid), (18.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::Trunk,
            TherapyZone::MixedFlow,
            None,
            None,
        );
        let side1 = (req.main_width_m * (1.0 - split1_center_frac) * 0.5).max(req.throat_width_m);
        let ctr1 = (req.main_width_m * split1_center_frac).max(req.throat_width_m);
        for (id, to, y) in [
            ("upper_trunk", "split2_upper", y_mid - 21.0),
            ("center_trunk", "split2", y_mid),
            ("lower_trunk", "split2_lower", y_mid + 21.0),
        ] {
            let zone = if id == "center_trunk" {
                TherapyZone::CancerTarget
            } else {
                TherapyZone::HealthyBypass
            };
            let width = if id == "center_trunk" { ctr1 } else { side1 };
            self.add_channel(
                id,
                "split1",
                to,
                vec![(18.0, y_mid), (36.0, y)],
                width,
                req.branch_length_m,
                req.channel_height_m,
                if id == "center_trunk" {
                    ChannelVisualRole::CenterTreatment
                } else {
                    ChannelVisualRole::PeripheralBypass
                },
                zone,
                None,
                None,
            );
        }
        let upper = self.leaf_triplet(y_mid - 21.0, 4.0);
        let center = self.leaf_triplet(y_mid, 7.0);
        let lower = self.leaf_triplet(y_mid + 21.0, 4.0);
        let outer_leaf_w = (side1 * (1.0 - split2_center_frac) * 0.5).max(req.throat_width_m);
        let outer_ctr_w = (side1 * split2_center_frac).max(req.throat_width_m);
        for ((prefix, split), ys) in [
            (("upper", "split2_upper"), upper),
            (("lower", "split2_lower"), lower),
        ] {
            let m1_name = format!("merge_{prefix}_1");
            let m2_name = format!("merge_{prefix}_2");
            let dir = if prefix == "upper" { -1.0 } else { 1.0 };
            let m1_y = ys[1] + 2.0 * dir;
            let m2_y = ys[1];

            self.add_node(&m1_name, (96.0, m1_y), NodeKind::Junction, Some(JunctionFamily::Merge));
            self.add_node(&m2_name, (104.0, m2_y), NodeKind::Junction, Some(JunctionFamily::Merge));

            let start_pt = self.nodes.iter().find(|n| n.id == self.node_ids[split]).unwrap().point;

            self.add_channel(
                &format!("{prefix}_L"), split, &m1_name,
                vec![start_pt, (54.0, ys[0]), (90.0, ys[0]), (96.0, m1_y)], outer_leaf_w,
                req.branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None,
            );

            self.add_channel(
                &format!("{prefix}_C"), split, &m1_name,
                vec![start_pt, (54.0, ys[1]), (90.0, ys[1]), (96.0, m1_y)], outer_ctr_w,
                req.branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None,
            );

            self.add_channel(
                &format!("{prefix}_R"), split, &m2_name,
                vec![start_pt, (54.0, ys[2]), (90.0, ys[2]), (104.0, m2_y)], outer_leaf_w,
                req.branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None,
            );

            let m1_to_m2_w = (outer_leaf_w + outer_ctr_w).min(req.throat_width_m * 4.0);
            self.add_channel(
                &format!("{prefix}_m1_to_m2"), &m1_name, &m2_name,
                vec![(96.0, m1_y), (104.0, m2_y)], m1_to_m2_w,
                req.branch_length_m * 0.1, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None,
            );

            self.add_channel(
                &format!("{prefix}_to_bypass"), &m2_name, "periph_merge",
                vec![(104.0, m2_y), (104.0, y_mid)], side1,
                req.branch_length_m * 0.2, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None,
            );
        }
        let center_side_w = (ctr1 * (1.0 - split2_center_frac) * 0.5).max(req.throat_width_m);
        let center_ctr_w = (ctr1 * split2_center_frac).max(req.throat_width_m);
        for (suffix, y, width) in [
            ("L", center[0], center_side_w),
            ("C", center[1], center_ctr_w),
            ("R", center[2], center_side_w),
        ] {
            let sub = format!("center_{suffix}");
            let sub_in = format!("{sub}_in");
            self.add_node(&sub_in, (54.0, y), NodeKind::Junction, None);
            self.add_channel(
                &format!("{sub}_approach"),
                "split2",
                &sub_in,
                vec![(36.0, y_mid), (54.0, y)],
                width,
                req.branch_length_m * 0.35,
                req.channel_height_m,
                ChannelVisualRole::CenterTreatment,
                TherapyZone::CancerTarget,
                None,
                None,
            );
            if center_throat_count == 0 {
                self.add_channel(
                    &format!("{sub}_treatment"),
                    &sub_in,
                    "merge_center",
                    vec![(54.0, y), (90.0, y), (100.0, y_mid)],
                    width,
                    req.branch_length_m,
                    req.channel_height_m,
                    ChannelVisualRole::CenterTreatment,
                    TherapyZone::CancerTarget,
                    None,
                    None,
                );
                continue;
            }
            let mut prev = sub_in.clone();
            for k in 0..center_throat_count {
                let tin = format!("{sub}_t{k}_in");
                let tout = format!("{sub}_t{k}_out");
                let x0 = 60.0 + f64::from(k) * 8.0;
                self.add_node(&tin, (x0, y), NodeKind::Junction, None);
                self.add_node(&tout, (x0 + 4.0, y), NodeKind::Junction, None);
                let start_pt = self
                    .nodes
                    .iter()
                    .find(|n| n.id == self.node_ids[&prev])
                    .unwrap()
                    .point;
                self.add_channel(
                    &format!("{sub}_conv{k}"),
                    &prev,
                    &tin,
                    vec![start_pt, (x0, y)],
                    width,
                    req.branch_length_m * 0.12,
                    req.channel_height_m,
                    ChannelVisualRole::CenterTreatment,
                    TherapyZone::CancerTarget,
                    None,
                    None,
                );
                self.add_channel(
                    &format!("{sub}_throat{k}"),
                    &tin,
                    &tout,
                    vec![(x0, y), (x0 + 2.0, y), (x0 + 4.0, y)],
                    req.throat_width_m,
                    req.throat_length_m,
                    req.channel_height_m,
                    ChannelVisualRole::VenturiThroat,
                    TherapyZone::CancerTarget,
                    None,
                    Some(VenturiGeometryMetadata {
                        throat_width_m: req.throat_width_m,
                        throat_height_m: req.channel_height_m,
                        throat_length_m: req.throat_length_m,
                        inlet_width_m: width,
                        outlet_width_m: width,
                        convergent_half_angle_deg: 15.0,
                        divergent_half_angle_deg: 15.0,
                        throat_position: 0.5,
                    }),
                );
                prev = tout;
                if k + 1 < center_throat_count {
                    let rdv = format!("{sub}_rdv{k}");
                    self.add_node(&rdv, (x0 + 6.0, y), NodeKind::Junction, None);
                    self.add_channel(
                        &format!("{sub}_diffuser{k}"),
                        &prev,
                        &rdv,
                        vec![(x0 + 4.0, y), (x0 + 6.0, y)],
                        width,
                        inter_spacing_m,
                        req.channel_height_m,
                        ChannelVisualRole::Diffuser,
                        TherapyZone::CancerTarget,
                        None,
                        None,
                    );
                    prev = rdv;
                }
            }
            let start_pt = self
                .nodes
                .iter()
                .find(|n| n.id == self.node_ids[&prev])
                .unwrap()
                .point;
            self.add_channel(
                &format!("{sub}_recovery"),
                &prev,
                "merge_center",
                vec![start_pt, (90.0, y), (100.0, y_mid)],
                width,
                req.branch_length_m * 0.18,
                req.channel_height_m,
                ChannelVisualRole::Diffuser,
                TherapyZone::CancerTarget,
                None,
                None,
            );
        }
        // Removed bypass_to_merge since upper and lower sweep straight to periph_merge
        self.add_channel(
            "periph_to_merge",
            "periph_merge",
            "merge_out",
            vec![(104.0, y_mid), (118.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::MergeCollector,
            TherapyZone::HealthyBypass,
            None,
            None,
        );
        self.add_channel(
            "center_to_merge",
            "merge_center",
            "merge_out",
            vec![(100.0, y_mid), (118.0, y_mid)],
            ctr1,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::MergeCollector,
            TherapyZone::MixedFlow,
            None,
            None,
        );
        self.add_channel(
            "outlet_section",
            "merge_out",
            "outlet",
            vec![(118.0, y_mid), (self.box_dims.0, y_mid)],
            req.main_width_m,
            req.trunk_length_m,
            req.channel_height_m,
            ChannelVisualRole::Trunk,
            TherapyZone::MixedFlow,
            None,
            None,
        );
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimitiveSelectiveSplitKind {
    Bi,
    Tri,
    Quad,
    Penta,
}

#[derive(Debug, Clone, PartialEq)]
pub struct PrimitiveSelectiveTreeRequest {
    pub name: String,
    pub box_dims_mm: (f64, f64),
    pub split_sequence: Vec<PrimitiveSelectiveSplitKind>,
    pub main_width_m: f64,
    pub throat_width_m: f64,
    pub throat_length_m: f64,
    pub channel_height_m: f64,
    pub first_trifurcation_center_frac: f64,
    pub later_trifurcation_center_frac: f64,
    pub bifurcation_treatment_frac: f64,
    pub treatment_branch_venturi_enabled: bool,
    pub treatment_branch_throat_count: u8,
    pub center_serpentine: Option<CenterSerpentinePathSpec>,
}

pub fn create_primitive_selective_tree_geometry(
    request: &PrimitiveSelectiveTreeRequest,
) -> NetworkBlueprint {
    // All primitive selective split sequences use the canonical full-tree
    // generator and are then annotated with treatment/bypass metadata.
    // This preserves the actual split topology instead of collapsing all-tri
    // trees into a centerline cascade surrogate.
    let splits: Vec<SplitType> = request
        .split_sequence
        .iter()
        .enumerate()
        .map(|(idx, split)| match split {
            PrimitiveSelectiveSplitKind::Bi => SplitType::Bifurcation,
            PrimitiveSelectiveSplitKind::Tri => SplitType::SymmetricTrifurcation {
                center_ratio: if idx == 0 {
                    request.first_trifurcation_center_frac
                } else {
                    request.later_trifurcation_center_frac
                },
            },
            PrimitiveSelectiveSplitKind::Quad => SplitType::Quadfurcation,
            PrimitiveSelectiveSplitKind::Penta => SplitType::Pentafurcation,
        })
        .collect();

    let geometry_config = GeometryConfig {
        wall_clearance: 1.0,
        channel_width: (request.main_width_m * 1.0e3).clamp(0.2, 12.0),
        channel_height: (request.channel_height_m * 1.0e3).clamp(0.2, 5.0),
        ..Default::default()
    };

    let split_depth = splits.len().max(1);
    let total_branches = splits
        .iter()
        .fold(1usize, |product, split| product.saturating_mul(split.branch_count()));
    let retry_scales = [1.0, 1.2, 1.45, 1.75, 2.1, 2.5];
    let mut last_blueprint = None;

    for retry_scale in retry_scales {
        let internal_dims = expanded_generation_box_dims(
            request.box_dims_mm,
            split_depth,
            total_branches,
            retry_scale,
        );
        let mut blueprint = super::create_geometry(
            internal_dims,
            &splits,
            &geometry_config,
            &ChannelTypeConfig::AllStraight,
        );
        if internal_dims != request.box_dims_mm {
            scale_blueprint_geometry(&mut blueprint, request.box_dims_mm);
        }
        annotate_primitive_tree(&mut blueprint, request);
        if blueprint.unresolved_channel_overlap_count() == 0 {
            return blueprint;
        }
        last_blueprint = Some(blueprint);
    }

    last_blueprint.expect("primitive selective tree generation should produce a blueprint")
}

fn expanded_generation_box_dims(
    target_dims: (f64, f64),
    split_depth: usize,
    total_branches: usize,
    retry_scale: f64,
) -> (f64, f64) {
    let depth_scale = (1.0 + 0.18 * split_depth.saturating_sub(1) as f64) * retry_scale;
    let branch_scale = ((total_branches as f64) / 9.0).sqrt().clamp(1.0, 3.2) * retry_scale;
    let width = target_dims.0 * depth_scale.max(1.0);
    let height = target_dims.1 * branch_scale.max(1.0);
    (width, height)
}

fn scale_blueprint_geometry(blueprint: &mut NetworkBlueprint, target_dims: (f64, f64)) {
    let source_dims = blueprint.box_dims;
    if source_dims.0 <= 0.0 || source_dims.1 <= 0.0 {
        return;
    }
    let scale_x = target_dims.0 / source_dims.0;
    let scale_y = target_dims.1 / source_dims.1;

    let scale_point = |(x, y): Point2D| (x * scale_x, y * scale_y);

    for node in &mut blueprint.nodes {
        node.point = scale_point(node.point);
    }
    for channel in &mut blueprint.channels {
        channel.path = channel.path.iter().copied().map(scale_point).collect();
    }
    blueprint.box_outline = blueprint
        .box_outline
        .iter()
        .map(|(start, end)| (scale_point(*start), scale_point(*end)))
        .collect();
    blueprint.box_dims = target_dims;
}

/// Build a primitive selective tree directly from a declarative topology spec
/// when the spec matches the canonical selective-routing contract.
///
/// Returns `Ok(None)` if the spec cannot be represented without losing
/// topology information.
pub fn create_primitive_selective_tree_geometry_from_spec(
    spec: &BlueprintTopologySpec,
) -> Result<Option<NetworkBlueprint>, String> {
    if spec.split_stages.is_empty() || spec.has_series_path() || spec.has_parallel_paths() {
        return Ok(None);
    }

    let mut split_sequence = Vec::with_capacity(spec.split_stages.len());
    let mut parent_width_m = spec.inlet_width_m;
    let mut first_trifurcation_center_frac = 0.5;
    let mut later_trifurcation_center_frac = 0.5;
    let mut bifurcation_treatment_frac = 0.5;
    let mut saw_later_trifurcation = false;
    let mut center_serpentine = None;

    for (stage_index, stage) in spec.split_stages.iter().enumerate() {
        let treatment_branches: Vec<_> = stage
            .branches
            .iter()
            .filter(|branch| branch.treatment_path)
            .collect();
        if treatment_branches.len() != 1 {
            return Ok(None);
        }
        let treatment_branch = treatment_branches[0];
        if parent_width_m <= 0.0 || treatment_branch.route.width_m <= 0.0 {
            return Ok(None);
        }

        let treatment_fraction = (treatment_branch.route.width_m / parent_width_m).clamp(0.0, 1.0);
        center_serpentine =
            center_serpentine.or(treatment_branch
                .route
                .serpentine
                .as_ref()
                .map(|serpentine| CenterSerpentinePathSpec {
                    segments: serpentine.segments,
                    bend_radius_m: serpentine.bend_radius_m,
                    wave_type: serpentine.wave_type,
                }));

        match stage.split_kind {
            SplitKind::NFurcation(2) => {
                split_sequence.push(PrimitiveSelectiveSplitKind::Bi);
                bifurcation_treatment_frac = treatment_fraction;
            }
            SplitKind::NFurcation(3) => {
                split_sequence.push(PrimitiveSelectiveSplitKind::Tri);
                if stage_index == 0 {
                    first_trifurcation_center_frac = treatment_fraction;
                } else {
                    later_trifurcation_center_frac = treatment_fraction;
                    saw_later_trifurcation = true;
                }
            }
            SplitKind::NFurcation(4) => {
                split_sequence.push(PrimitiveSelectiveSplitKind::Quad);
            }
            SplitKind::NFurcation(5) => {
                split_sequence.push(PrimitiveSelectiveSplitKind::Penta);
            }
            SplitKind::NFurcation(_) => return Ok(None),
        }

        parent_width_m = treatment_branch.route.width_m;
    }

    if !saw_later_trifurcation {
        later_trifurcation_center_frac = first_trifurcation_center_frac;
    }

    let representative_height_m = spec
        .split_stages
        .first()
        .and_then(|stage| stage.branches.first())
        .map_or(1.0e-3, |branch| branch.route.height_m);
    let strongest_venturi = spec
        .venturi_placements
        .iter()
        .max_by_key(|placement| placement.serial_throat_count);

    let request = PrimitiveSelectiveTreeRequest {
        name: spec.design_name.clone(),
        box_dims_mm: spec.box_dims_mm,
        split_sequence,
        main_width_m: spec.inlet_width_m,
        throat_width_m: strongest_venturi.map_or(parent_width_m, |placement| {
            placement.throat_geometry.throat_width_m
        }),
        throat_length_m: strongest_venturi.map_or(spec.trunk_length_m / 8.0, |placement| {
            placement.throat_geometry.throat_length_m
        }),
        channel_height_m: representative_height_m,
        first_trifurcation_center_frac,
        later_trifurcation_center_frac,
        bifurcation_treatment_frac,
        treatment_branch_venturi_enabled: spec.treatment_mode
            == TreatmentActuationMode::VenturiCavitation
            && !spec.venturi_placements.is_empty(),
        treatment_branch_throat_count: strongest_venturi
            .map_or(0, |placement| placement.serial_throat_count),
        center_serpentine,
    };

    let mut blueprint = create_primitive_selective_tree_geometry(&request);
    blueprint.name.clone_from(&spec.design_name);
    blueprint.topology = Some(spec.clone());
    blueprint.lineage = Some(BlueprintTopologyFactory::lineage_for_spec(spec));
    blueprint.render_hints = Some(BlueprintRenderHints {
        stage_sequence: spec.stage_sequence_label(),
        split_layers: spec.visible_split_layers(),
        throat_count_hint: spec.venturi_count(),
        treatment_label: if spec.treatment_mode == TreatmentActuationMode::VenturiCavitation {
            "venturi".to_string()
        } else {
            "ultrasound".to_string()
        },
        mirror_x: false,
        mirror_y: false,
    });
    Ok(Some(blueprint))
}

fn channel_length_from_points_or_endpoints(
    points: &[Point2D],
    start: Option<Point2D>,
    end: Option<Point2D>,
    fallback_m: f64,
) -> f64 {
    let path_length_m = if points.len() >= 2 {
        polyline_length_mm(points) * 1.0e-3
    } else {
        0.0
    };
    if path_length_m.is_finite() && path_length_m > 0.0 {
        return path_length_m;
    }

    if let (Some(start), Some(end)) = (start, end) {
        let endpoint_length_m = ((end.0 - start.0).hypot(end.1 - start.1)) * 1.0e-3;
        if endpoint_length_m.is_finite() && endpoint_length_m > 0.0 {
            return endpoint_length_m;
        }
    }

    if fallback_m.is_finite() && fallback_m > 0.0 {
        fallback_m
    } else {
        f64::EPSILON
    }
}

fn annotate_primitive_tree(system: &mut NetworkBlueprint, request: &PrimitiveSelectiveTreeRequest) {
    if system.nodes.is_empty() || system.channels.is_empty() {
        return;
    }

    let mid_x = system.box_dims.0 * 0.5;
    let mid_y = system.box_dims.1 * 0.5;
    let inlet_node = system
        .nodes
        .iter()
        .min_by(|a, b| a.point.0.total_cmp(&b.point.0))
        .map(|node| node.id.clone());
    let outlet_node = system
        .nodes
        .iter()
        .max_by(|a, b| a.point.0.total_cmp(&b.point.0))
        .map(|node| node.id.clone());
    let (Some(inlet_node), Some(outlet_node)) = (inlet_node, outlet_node) else {
        return;
    };

    let mut adjacency: HashMap<NodeId, Vec<(NodeId, usize)>> = HashMap::new();
    for (idx, channel) in system.channels.iter().enumerate() {
        adjacency
            .entry(channel.from.clone())
            .or_default()
            .push((channel.to.clone(), idx));
        adjacency
            .entry(channel.to.clone())
            .or_default()
            .push((channel.from.clone(), idx));
    }

    let node_pts: HashMap<NodeId, Point2D> = system
        .nodes
        .iter()
        .map(|n| (n.id.clone(), n.point))
        .collect();

    let leaf_nodes: Vec<NodeId> = system
        .nodes
        .iter()
        .filter(|node| (node.point.0 - mid_x).abs() < 1e-6)
        .map(|node| node.id.clone())
        .collect();
    let mut ordered_leaves = leaf_nodes;
    ordered_leaves.sort_by(|a, b| node_pts[a].1.total_cmp(&node_pts[b].1));

    let treatment_leaf_indices = primitive_treatment_leaf_indices(
        ordered_leaves.len(),
        request.split_sequence.first().copied(),
    );
    let treatment_leaves: HashSet<NodeId> = treatment_leaf_indices
        .into_iter()
        .filter_map(|idx| ordered_leaves.get(idx).cloned())
        .collect();

    let mut treatment_channels = HashSet::new();
    for leaf in &treatment_leaves {
        if let Some(path) =
            channel_path_between(&inlet_node, leaf, &node_pts, &adjacency, Some(mid_x), true)
        {
            treatment_channels.extend(path);
        }
        if let Some(path) = channel_path_between(
            leaf,
            &outlet_node,
            &node_pts,
            &adjacency,
            Some(mid_x),
            false,
        ) {
            treatment_channels.extend(path);
        }
    }

    let first_tri = matches!(
        request.split_sequence.first(),
        Some(PrimitiveSelectiveSplitKind::Tri)
    );
    let bypass_width = if first_tri {
        (request.main_width_m * (1.0 - request.first_trifurcation_center_frac) * 0.5)
            .max(request.throat_width_m)
    } else {
        (request.main_width_m * (1.0 - request.bifurcation_treatment_frac))
            .max(request.throat_width_m)
    };
    let treatment_width = if first_tri {
        (request.main_width_m * request.first_trifurcation_center_frac).max(request.throat_width_m)
    } else {
        (request.main_width_m * request.bifurcation_treatment_frac).max(request.throat_width_m)
    };


    for node in &mut system.nodes {
        let degree = adjacency.get(&node.id).map_or(0, Vec::len);
        node.kind = if node.id == inlet_node {
            NodeKind::Inlet
        } else if node.id == outlet_node {
            NodeKind::Outlet
        } else {
            NodeKind::Junction
        };

        if degree >= 3 {
            let junction_family = if degree >= 4 {
                JunctionFamily::Trifurcation
            } else {
                JunctionFamily::Bifurcation
            };
            node.junction_geometry = Some(JunctionGeometryMetadata {
                junction_family,
                branch_angles_deg: Vec::new(),
                merge_angles_deg: Vec::new(),
            });
        }
    }

    let mut pending_venturi_paths = Vec::new();

    for (idx, channel) in system.channels.iter_mut().enumerate() {
        let points = channel.path.clone();
        let start_point = node_pts.get(&channel.from).copied();
        let end_point = node_pts.get(&channel.to).copied();
        let (min_x, max_x) = if points.is_empty() {
            // Straight channels have empty path; derive x-range from endpoint node positions
            let x0 = start_point.map_or(f64::INFINITY, |(x, _)| x);
            let x1 = end_point.map_or(f64::NEG_INFINITY, |(x, _)| x);
            (x0.min(x1), x0.max(x1))
        } else {
            let min_x = points.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
            let max_x = points
                .iter()
                .map(|(x, _)| *x)
                .fold(f64::NEG_INFINITY, f64::max);
            (min_x, max_x)
        };
        let avg_y = if points.is_empty() {
            mid_y
        } else {
            points.iter().map(|(_, y)| *y).sum::<f64>() / points.len() as f64
        };

        let is_treatment = treatment_channels.contains(&idx);
        let touches_inlet = channel.from == inlet_node || channel.to == inlet_node;
        let touches_outlet = channel.from == outlet_node || channel.to == outlet_node;
        let starts_at_treatment_leaf = treatment_leaves.contains(&channel.from);
        let is_treatment_window_channel =
            starts_at_treatment_leaf && max_x > mid_x + 1e-6 && min_x >= mid_x - 1e-6;
        let _touches_treatment_leaf =
            starts_at_treatment_leaf || treatment_leaves.contains(&channel.to);
        let is_trunk = touches_inlet || touches_outlet;

        let physical_width = if touches_inlet || touches_outlet {
            request.main_width_m
        } else if is_treatment {
            treatment_width
        } else {
            bypass_width
        };

        channel.cross_section = crate::domain::model::CrossSectionSpec::Rectangular {
            width_m: physical_width,
            height_m: request.channel_height_m,
        };
        channel.length_m = channel_length_from_points_or_endpoints(
            &points,
            start_point,
            end_point,
            channel.length_m,
        );

        let mut role = if is_treatment {
            ChannelVisualRole::CenterTreatment
        } else {
            ChannelVisualRole::PeripheralBypass
        };
        if touches_inlet || touches_outlet {
            role = ChannelVisualRole::Trunk;
        } else if !is_treatment && max_x > mid_x {
            role = ChannelVisualRole::MergeCollector;
        }
        channel.visual_role = Some(role);
        channel.therapy_zone = Some(if is_trunk {
            TherapyZone::MixedFlow
        } else if is_treatment {
            TherapyZone::CancerTarget
        } else {
            TherapyZone::HealthyBypass
        });

        let should_overlay_serpentine = is_treatment
            && !is_trunk
            && request.center_serpentine.is_some()
            && (!request.treatment_branch_venturi_enabled || !is_treatment_window_channel);

        if should_overlay_serpentine {
            if let Some(spec) = request.center_serpentine {
                // When the generator produced a straight channel (empty path), seed
                // the serpentine overlay with the endpoint node positions so the
                // overlay has a valid bounding box to work with.
                let source_points: Vec<Point2D> = if points.is_empty() {
                    [start_point, end_point].into_iter().flatten().collect()
                } else {
                    points.clone()
                };
                let serpentine_path =
                    serpentine_overlay_path(
                        &source_points,
                        physical_width,
                        spec,
                        if request.treatment_branch_venturi_enabled {
                            1.1
                        } else {
                            1.5
                        },
                    );
                channel.path = serpentine_path;
                channel.length_m = channel_length_from_points_or_endpoints(
                    &channel.path,
                    start_point,
                    end_point,
                    channel.length_m,
                );
                if let Some((start, end)) = channel.path.first().zip(channel.path.last()) {
                    channel.channel_shape =
                        infer_serpentine_shape(&channel.path, *start, *end, physical_width * 1.0e3);
                }
            }
        }

        if is_treatment && is_treatment_window_channel && request.treatment_branch_venturi_enabled {
            pending_venturi_paths.push(PendingVenturiPath {
                channel_idx: idx,
                start: points.first().copied().or(start_point),
                end: points.last().copied().or(end_point),
                preferred_y: preferred_treatment_lane_y(&points, start_point, end_point, avg_y),
                fallback_length_m: channel.length_m,
            });
            channel.visual_role = Some(ChannelVisualRole::VenturiThroat);
            channel.venturi_geometry = Some(VenturiGeometryMetadata {
                throat_width_m: request.throat_width_m,
                throat_height_m: request.channel_height_m,
                throat_length_m: request.throat_length_m,
                inlet_width_m: physical_width,
                outlet_width_m: physical_width,
                convergent_half_angle_deg: 15.0,
                divergent_half_angle_deg: 15.0,
                throat_position: 0.5,
            });
            let metadata = channel.metadata.get_or_insert_with(MetadataContainer::new);
            metadata.insert(ChannelVenturiSpec {
                n_throats: request.treatment_branch_throat_count.max(1),
                is_ctc_stream: true,
                throat_width_m: request.throat_width_m,
                height_m: request.channel_height_m,
                inter_throat_spacing_m: request.throat_length_m * 2.0,
            });
            metadata.insert(channel.venturi_geometry.clone().unwrap());
        } else if !is_treatment {
            let metadata = channel.metadata.get_or_insert_with(MetadataContainer::new);
            metadata.insert(ChannelVenturiSpec {
                n_throats: 0,
                is_ctc_stream: false,
                throat_width_m: request.throat_width_m,
                height_m: request.channel_height_m,
                inter_throat_spacing_m: request.throat_length_m * 2.0,
            });
        }
    }

    route_pending_venturi_paths(&mut system.channels, &pending_venturi_paths, mid_y);
}

fn primitive_treatment_leaf_indices(
    leaf_count: usize,
    first_stage: Option<PrimitiveSelectiveSplitKind>,
) -> Vec<usize> {
    if leaf_count == 0 {
        return Vec::new();
    }
    match first_stage {
        Some(PrimitiveSelectiveSplitKind::Bi) => {
            // Leaves sorted by Y ascending: index 0 = lower (bypass), index 1 = upper (treatment)
            vec![leaf_count - 1]
        }
        Some(PrimitiveSelectiveSplitKind::Tri) => {
            if leaf_count == 3 {
                vec![1]
            } else {
                let block = leaf_count / 3;
                (block..(2 * block)).collect()
            }
        }
        Some(PrimitiveSelectiveSplitKind::Quad) => {
            // 4 leaves sorted by Y: indices 1,2 are the center pair (treatment)
            let quarter = leaf_count / 4;
            (quarter..(3 * quarter)).collect()
        }
        Some(PrimitiveSelectiveSplitKind::Penta) => {
            // 5 leaves sorted by Y: index 2 is the center (treatment)
            if leaf_count == 5 {
                vec![2]
            } else {
                let fifth = leaf_count / 5;
                (2 * fifth..(3 * fifth)).collect()
            }
        }
        None => Vec::new(),
    }
}

use crate::domain::model::NodeId;

fn channel_path_between(
    start: &NodeId,
    goal: &NodeId,
    node_pts: &HashMap<NodeId, Point2D>,
    adjacency: &HashMap<NodeId, Vec<(NodeId, usize)>>,
    mid_x: Option<f64>,
    left_half: bool,
) -> Option<Vec<usize>> {
    let mut visited = HashSet::new();
    let mut path = Vec::new();
    if dfs_channel_path(
        start,
        goal,
        node_pts,
        adjacency,
        mid_x,
        left_half,
        &mut visited,
        &mut path,
    ) {
        Some(path)
    } else {
        None
    }
}

fn dfs_channel_path(
    current: &NodeId,
    goal: &NodeId,
    node_pts: &HashMap<NodeId, Point2D>,
    adjacency: &HashMap<NodeId, Vec<(NodeId, usize)>>,
    mid_x: Option<f64>,
    left_half: bool,
    visited: &mut HashSet<NodeId>,
    path: &mut Vec<usize>,
) -> bool {
    if current == goal {
        return true;
    }
    visited.insert(current.clone());
    let Some(neighbors) = adjacency.get(current) else {
        return false;
    };
    for (next, channel_idx) in neighbors {
        if visited.contains(next) {
            continue;
        }
        if let Some(mid_x) = mid_x {
            let from = node_pts[current].0;
            let to = node_pts[next].0;
            let keep = if left_half {
                from <= mid_x + 1e-6 && to <= mid_x + 1e-6
            } else {
                from >= mid_x - 1e-6 && to >= mid_x - 1e-6
            };
            if !keep {
                continue;
            }
        }
        path.push(*channel_idx);
        if dfs_channel_path(
            next, goal, node_pts, adjacency, mid_x, left_half, visited, path,
        ) {
            return true;
        }
        path.pop();
    }
    visited.remove(current);
    false
}

#[cfg(test)]
mod tests {
    use super::{
        create_primitive_selective_tree_geometry, path_intersects_any,
        route_monotone_treatment_path, CenterSerpentinePathSpec, PrimitiveSelectiveSplitKind,
        PrimitiveSelectiveTreeRequest,
    };

    #[test]
    fn primitive_selective_tree_annotation_preserves_positive_channel_lengths() {
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "primitive-selective-lengths".to_string(),
            box_dims_mm: (127.76, 85.47),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: 8.0e-3,
            throat_width_m: 55.0e-6,
            throat_length_m: 110.0e-6,
            channel_height_m: 1.0e-3,
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: false,
            treatment_branch_throat_count: 1,
            center_serpentine: None,
        });

        let invalid_lengths: Vec<String> = blueprint
            .channels
            .iter()
            .filter(|channel| !(channel.length_m.is_finite() && channel.length_m > 0.0))
            .map(|channel| format!("{}={}", channel.id.as_str(), channel.length_m))
            .collect();

        assert!(
            invalid_lengths.is_empty(),
            "primitive selective annotation should preserve positive channel lengths: {}",
            invalid_lengths.join(", ")
        );
    }

    #[test]
    fn primitive_selective_venturi_paths_have_no_unresolved_crossings() {
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "primitive-selective-no-crossings".to_string(),
            box_dims_mm: (127.76, 85.47),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: 8.0e-3,
            throat_width_m: 55.0e-6,
            throat_length_m: 110.0e-6,
            channel_height_m: 1.0e-3,
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: true,
            treatment_branch_throat_count: 1,
            center_serpentine: None,
        });

        assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
        blueprint
            .validate()
            .expect("venturi treatment paths must remain planar");

        for channel in blueprint.venturi_channels() {
            if let (Some(start), Some(end)) = (channel.path.first(), channel.path.last()) {
                if (start.1 - end.1).abs() < 1e-9 {
                    assert!(
                        channel
                            .path
                            .iter()
                            .all(|point| (point.1 - start.1).abs() < 1e-9),
                        "equal-y treatment channel {} must stay on its own lane",
                        channel.id.as_str()
                    );
                }
            }
        }
    }

    #[test]
    fn primitive_selective_penta_quad_tri_is_geometry_authorable() {
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "primitive-selective-penta-quad-tri".to_string(),
            box_dims_mm: (127.76, 85.47),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Penta,
                PrimitiveSelectiveSplitKind::Quad,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: 4.0e-3,
            throat_width_m: 35.0e-6,
            throat_length_m: 180.0e-6,
            channel_height_m: 1.0e-3,
            first_trifurcation_center_frac: 0.45,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: true,
            treatment_branch_throat_count: 1,
            center_serpentine: None,
        });

        assert_eq!(blueprint.box_dims, (127.76, 85.47));
        assert_eq!(blueprint.unresolved_channel_overlap_count(), 0);
        blueprint
            .validate()
            .expect("Penta->Quad->Tri primitive selective tree should remain planar");
    }

    #[test]
    fn primitive_selective_tri_tri_retains_multiple_treatment_window_lanes() {
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "primitive-selective-tritri-lanes".to_string(),
            box_dims_mm: (127.76, 85.47),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: 8.0e-3,
            throat_width_m: 55.0e-6,
            throat_length_m: 110.0e-6,
            channel_height_m: 1.0e-3,
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: false,
            treatment_branch_throat_count: 1,
            center_serpentine: None,
        });

        let mid_x = blueprint.box_dims.0 * 0.5;
        let mut lane_keys = std::collections::BTreeSet::new();
        for channel in &blueprint.channels {
            if channel.therapy_zone
                != Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)
            {
                continue;
            }
            let points = if channel.path.is_empty() {
                let start = blueprint
                    .nodes
                    .iter()
                    .find(|node| node.id == channel.from)
                    .map(|node| node.point);
                let end = blueprint
                    .nodes
                    .iter()
                    .find(|node| node.id == channel.to)
                    .map(|node| node.point);
                [start, end].into_iter().flatten().collect::<Vec<_>>()
            } else {
                channel.path.clone()
            };
            if points.is_empty() {
                continue;
            }
            let min_x = points.iter().map(|(x, _)| *x).fold(f64::INFINITY, f64::min);
            let max_x = points
                .iter()
                .map(|(x, _)| *x)
                .fold(f64::NEG_INFINITY, f64::max);
            if max_x <= mid_x + 1.0e-6 || min_x < mid_x - 1.0e-6 {
                continue;
            }
            let mean_y = points.iter().map(|(_, y)| *y).sum::<f64>() / points.len() as f64;
            lane_keys.insert((mean_y * 100.0).round() as i64);
        }

        assert!(
            lane_keys.len() >= 3,
            "Tri->Tri selective trees must preserve multiple treatment-window lanes instead of collapsing to a centerline surrogate, got {:?}",
            lane_keys
        );
    }

    #[test]
    fn monotone_treatment_routing_preserves_equal_y_lane() {
        let routed =
            route_monotone_treatment_path(Some((10.0, 24.0)), Some((40.0, 24.0)), 24.0, 42.0, &[]);
        assert_eq!(routed, vec![(10.0, 24.0), (40.0, 24.0)]);
    }

    #[test]
    fn monotone_treatment_routing_doglegs_around_existing_branch() {
        let existing_paths = vec![vec![(3.0, 0.0), (7.0, 4.0)]];
        let routed = route_monotone_treatment_path(
            Some((0.0, 2.0)),
            Some((10.0, 2.0)),
            0.0,
            2.0,
            &existing_paths,
        );

        assert!(
            routed.len() > 2,
            "crossing direct path should reroute with a dogleg"
        );
        assert!(
            !path_intersects_any(&routed, &existing_paths),
            "rerouted path must not cross the existing treatment lane"
        );
    }

    #[test]
    fn primitive_selective_serpentine_mirrors_on_both_sides_of_midline() {
        use crate::domain::model::ChannelShape;
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "serp-mirror-check".to_string(),
            box_dims_mm: (127.76, 85.47),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: 4.0e-3,
            throat_width_m: 55.0e-6,
            throat_length_m: 110.0e-6,
            channel_height_m: 1.0e-3,
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: false,
            treatment_branch_throat_count: 1,
            center_serpentine: Some(CenterSerpentinePathSpec {
                segments: 5,
                bend_radius_m: 3.0e-3,
                wave_type: crate::SerpentineWaveType::default(),
            }),
        });

        let mid_x = blueprint.box_dims.0 * 0.5;
        let serpentine_channels: Vec<_> = blueprint
            .channels
            .iter()
            .filter(|ch| matches!(ch.channel_shape, ChannelShape::Serpentine { .. }))
            .collect();

        assert!(
            serpentine_channels.len() >= 2,
            "at least one serpentine channel on each side expected, got {}",
            serpentine_channels.len()
        );

        let node_pts: std::collections::HashMap<_, _> = blueprint
            .nodes
            .iter()
            .map(|n| (n.id.clone(), n.point))
            .collect();
        let left_count = serpentine_channels
            .iter()
            .filter(|ch| {
                let max_x = [
                    node_pts.get(&ch.from).map(|p| p.0),
                    node_pts.get(&ch.to).map(|p| p.0),
                ]
                .into_iter()
                .flatten()
                .fold(f64::NEG_INFINITY, f64::max);
                max_x <= mid_x + 1e-6
            })
            .count();
        let right_count = serpentine_channels
            .iter()
            .filter(|ch| {
                let min_x = [
                    node_pts.get(&ch.from).map(|p| p.0),
                    node_pts.get(&ch.to).map(|p| p.0),
                ]
                .into_iter()
                .flatten()
                .fold(f64::INFINITY, f64::min);
                min_x >= mid_x - 1e-6
            })
            .count();

        assert!(
            left_count >= 1,
            "split-side (left) treatment channels must have serpentine overlays, got {left_count}"
        );
        assert!(
            right_count >= 1,
            "merge-side (right) treatment channels must have serpentine overlays, got {right_count}"
        );
    }

    #[test]
    fn treatment_channels_retain_center_treatment_role_on_merge_side() {
        use crate::geometry::metadata::ChannelVisualRole;
        let blueprint = create_primitive_selective_tree_geometry(&PrimitiveSelectiveTreeRequest {
            name: "role-mirror-check".to_string(),
            box_dims_mm: (127.76, 85.47),
            split_sequence: vec![
                PrimitiveSelectiveSplitKind::Tri,
                PrimitiveSelectiveSplitKind::Tri,
            ],
            main_width_m: 4.0e-3,
            throat_width_m: 55.0e-6,
            throat_length_m: 110.0e-6,
            channel_height_m: 1.0e-3,
            first_trifurcation_center_frac: 0.55,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_branch_venturi_enabled: false,
            treatment_branch_throat_count: 1,
            center_serpentine: Some(CenterSerpentinePathSpec {
                segments: 5,
                bend_radius_m: 3.0e-3,
                wave_type: crate::SerpentineWaveType::default(),
            }),
        });

        let mid_x = blueprint.box_dims.0 * 0.5;
        let node_pts: std::collections::HashMap<_, _> = blueprint
            .nodes
            .iter()
            .map(|n| (n.id.clone(), n.point))
            .collect();

        let merge_side_treatment: Vec<_> = blueprint
            .channels
            .iter()
            .filter(|ch| {
                ch.therapy_zone
                    == Some(crate::domain::therapy_metadata::TherapyZone::CancerTarget)
            })
            .filter(|ch| {
                let max_x = [
                    node_pts.get(&ch.from).map(|p| p.0),
                    node_pts.get(&ch.to).map(|p| p.0),
                ]
                .into_iter()
                .flatten()
                .fold(f64::NEG_INFINITY, f64::max);
                max_x > mid_x + 1e-6
            })
            .collect();

        assert!(
            !merge_side_treatment.is_empty(),
            "should have treatment channels on the merge side"
        );
        for ch in &merge_side_treatment {
            assert_eq!(
                ch.visual_role,
                Some(ChannelVisualRole::CenterTreatment),
                "merge-side treatment channel {} must retain CenterTreatment role",
                ch.id.as_str()
            );
        }
    }
}

fn serpentine_overlay_path(
    points: &[Point2D],
    width_m: f64,
    spec: CenterSerpentinePathSpec,
    amplitude_mm: f64,
) -> Vec<Point2D> {
    let Some(start) = points.first().copied() else {
        return Vec::new();
    };
    let Some(end) = points.last().copied() else {
        return vec![start];
    };
    build_serpentine_lobe_path(start, end, width_m, spec, amplitude_mm)
}

fn preferred_treatment_lane_y(
    points: &[Point2D],
    start: Option<Point2D>,
    end: Option<Point2D>,
    y_fallback: f64,
) -> f64 {
    match (start, end) {
        (Some(a), Some(b)) => f64::midpoint(a.1, b.1),
        (Some(a), None) => a.1,
        (None, Some(b)) => b.1,
        (None, None) => match (points.first(), points.last()) {
            (Some(first), Some(last)) => f64::midpoint(first.1, last.1),
            (Some(point), None) | (None, Some(point)) => point.1,
            (None, None) => y_fallback,
        },
    }
}

fn route_pending_venturi_paths(
    channels: &mut [ChannelSpec],
    pending_paths: &[PendingVenturiPath],
    mid_y: f64,
) {
    let mut ordered_paths = pending_paths.to_vec();
    ordered_paths.sort_by(|left, right| {
        right
            .preferred_y
            .partial_cmp(&left.preferred_y)
            .unwrap_or(Ordering::Equal)
            .then(left.channel_idx.cmp(&right.channel_idx))
    });

    let mut assigned_paths = Vec::with_capacity(ordered_paths.len());
    for pending in ordered_paths {
        let routed_path = route_monotone_treatment_path(
            pending.start,
            pending.end,
            pending.preferred_y,
            mid_y,
            &assigned_paths,
        );
        let channel = &mut channels[pending.channel_idx];
        channel.path = routed_path;
        channel.length_m = channel_length_from_points_or_endpoints(
            &channel.path,
            pending.start,
            pending.end,
            pending.fallback_length_m,
        );
        assigned_paths.push(channel.path.clone());
    }
}

fn route_monotone_treatment_path(
    start: Option<Point2D>,
    end: Option<Point2D>,
    preferred_y: f64,
    mid_y: f64,
    existing_paths: &[Vec<Point2D>],
) -> Vec<Point2D> {
    match (start, end) {
        (Some(a), Some(b)) => {
            let direct = simplify_polyline_points(vec![a, b]);
            if !path_intersects_any(&direct, existing_paths) {
                return direct;
            }

            let direction = if preferred_y >= mid_y { 1.0 } else { -1.0 };
            let offset_step = ((a.1 - b.1).abs() * 0.5).max(2.0);
            for attempt in 0..8 {
                let lane_y = preferred_y + direction * f64::from(attempt) * offset_step;
                let dogleg = simplify_polyline_points(monotone_dogleg_path(a, b, lane_y));
                if !path_intersects_any(&dogleg, existing_paths) {
                    return dogleg;
                }
            }

            direct
        }
        (Some(a), None) => vec![a],
        (None, Some(b)) => vec![b],
        (None, None) => Vec::new(),
    }
}

fn monotone_dogleg_path(start: Point2D, end: Point2D, lane_y: f64) -> Vec<Point2D> {
    let span_x = end.0 - start.0;
    let lead_frac = if span_x.abs() > 1e-9 { 0.25 } else { 0.0 };
    let x1 = start.0 + span_x * lead_frac;
    let x2 = end.0 - span_x * lead_frac;
    vec![
        start,
        (x1, start.1),
        (x1, lane_y),
        (x2, lane_y),
        (x2, end.1),
        end,
    ]
}

fn simplify_polyline_points(points: Vec<Point2D>) -> Vec<Point2D> {
    let mut simplified = Vec::with_capacity(points.len());
    for point in points {
        if simplified.last().is_some_and(|last: &Point2D| {
            (last.0 - point.0).abs() < 1e-9 && (last.1 - point.1).abs() < 1e-9
        }) {
            continue;
        }
        simplified.push(point);
    }
    if simplified.len() < 3 {
        return simplified;
    }

    let mut cleaned = Vec::with_capacity(simplified.len());
    for point in simplified {
        cleaned.push(point);
        while cleaned.len() >= 3 {
            let len = cleaned.len();
            if points_are_colinear(cleaned[len - 3], cleaned[len - 2], cleaned[len - 1]) {
                cleaned.remove(len - 2);
            } else {
                break;
            }
        }
    }
    cleaned
}

fn points_are_colinear(a: Point2D, b: Point2D, c: Point2D) -> bool {
    let cross = (b.0 - a.0) * (c.1 - b.1) - (b.1 - a.1) * (c.0 - b.0);
    cross.abs() < 1e-9
}

fn segment_intersection_strict(p1: Point2D, p2: Point2D, p3: Point2D, p4: Point2D) -> bool {
    let dx1 = p2.0 - p1.0;
    let dy1 = p2.1 - p1.1;
    let dx2 = p4.0 - p3.0;
    let dy2 = p4.1 - p3.1;
    let denom = dx1 * dy2 - dy1 * dx2;
    if denom.abs() < 1e-12 {
        return false;
    }

    let t = ((p3.0 - p1.0) * dy2 - (p3.1 - p1.1) * dx2) / denom;
    let u = ((p3.0 - p1.0) * dy1 - (p3.1 - p1.1) * dx1) / denom;
    let eps = 1e-6;
    t > eps && t < (1.0 - eps) && u > eps && u < (1.0 - eps)
}

fn path_intersects_any(candidate: &[Point2D], existing_paths: &[Vec<Point2D>]) -> bool {
    candidate.windows(2).any(|segment_a| {
        existing_paths.iter().any(|path| {
            path.windows(2).any(|segment_b| {
                segment_intersection_strict(segment_a[0], segment_a[1], segment_b[0], segment_b[1])
            })
        })
    })
}

fn polyline_length_mm(points: &[Point2D]) -> f64 {
    points
        .windows(2)
        .map(|segment| {
            let dx = segment[1].0 - segment[0].0;
            let dy = segment[1].1 - segment[0].1;
            dx.hypot(dy)
        })
        .sum()
}

fn build_serpentine_lobe_path(
    start: Point2D,
    end: Point2D,
    width_m: f64,
    spec: CenterSerpentinePathSpec,
    amplitude_mm: f64,
) -> Vec<Point2D> {
    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let chord_length_mm = dx.hypot(dy);
    if chord_length_mm <= 1e-9 || spec.segments < 2 {
        return vec![start, end];
    }

    let channel_diameter_mm = (width_m * 1.0e3).max(1.0e-3);
    let bend_radius_mm = (spec.bend_radius_m * 1.0e3).max(channel_diameter_mm * 0.5);
    let requested_turns = spec.segments.saturating_sub(1).max(2);
    // Always materialize at least a mirrored 2-lobe rotated-S when a treatment
    // serpentine is requested. The physical bend radius still modulates the
    // shoulder and amplitude, but the authored path never degrades to a single apex.
    let safe_turns = requested_turns;
    let step_mm = chord_length_mm / safe_turns as f64;
    let shoulder_mm = (step_mm * 0.28)
        .max(channel_diameter_mm * 0.75)
        .min(bend_radius_mm.max(channel_diameter_mm * 1.5));
    let lobe_amplitude_mm = amplitude_mm
        .max(channel_diameter_mm * 1.5)
        .min(step_mm * 0.40);
    let tx = dx / chord_length_mm;
    let ty = dy / chord_length_mm;
    let nx = -ty;
    let ny = tx;

    let mut path = Vec::with_capacity(2 + safe_turns * 5);
    path.push(start);
    let mut last_along_mm = 0.0;
    for turn_idx in 0..safe_turns {
        let center_mm = step_mm * (turn_idx as f64 + 0.5);
        let direction = if turn_idx % 2 == 0 { 1.0 } else { -1.0 };
        for (along_mm, offset_scale) in [
            (center_mm - shoulder_mm, 0.35),
            (center_mm - shoulder_mm * 0.45, 0.72),
            (center_mm, 1.0),
            (center_mm + shoulder_mm * 0.45, 0.72),
            (center_mm + shoulder_mm, 0.35),
        ] {
            if along_mm <= last_along_mm + 1.0e-6 || along_mm >= chord_length_mm - 1.0e-6 {
                continue;
            }
            let offset_mm = direction * lobe_amplitude_mm * offset_scale;
            path.push((
                start.0 + tx * along_mm + nx * offset_mm,
                start.1 + ty * along_mm + ny * offset_mm,
            ));
            last_along_mm = along_mm;
        }
    }
    path.push(end);
    path
}

fn estimate_local_bend_radius(path: &[Point2D], idx: usize) -> Option<f64> {
    if idx == 0 || idx + 1 >= path.len() {
        return None;
    }

    let a = path[idx - 1];
    let b = path[idx];
    let c = path[idx + 1];
    let ab = ((b.0 - a.0).powi(2) + (b.1 - a.1).powi(2)).sqrt();
    let bc = ((c.0 - b.0).powi(2) + (c.1 - b.1).powi(2)).sqrt();
    let ac = ((c.0 - a.0).powi(2) + (c.1 - a.1).powi(2)).sqrt();
    if ab <= 1e-9 || bc <= 1e-9 || ac <= 1e-9 {
        return None;
    }

    let twice_area = ((b.0 - a.0) * (c.1 - a.1) - (b.1 - a.1) * (c.0 - a.0)).abs();
    if twice_area <= 1e-12 {
        return Some(f64::INFINITY);
    }

    Some((ab * bc * ac) / (2.0 * twice_area))
}

fn infer_serpentine_shape(
    path: &[Point2D],
    start: Point2D,
    end: Point2D,
    channel_width_mm: f64,
) -> ChannelShape {
    if path.len() < 3 {
        return ChannelShape::Serpentine {
            segments: 2,
            bend_radius_m: (channel_width_mm * 0.5) * 1.0e-3,
            wave_type: crate::topology::SerpentineWaveType::Sine,
        };
    }

    let dx = end.0 - start.0;
    let dy = end.1 - start.1;
    let length = dx.hypot(dy);
    if length <= 1e-9 {
        return ChannelShape::Serpentine {
            segments: 2,
            bend_radius_m: (channel_width_mm * 0.5) * 1.0e-3,
            wave_type: crate::topology::SerpentineWaveType::Sine,
        };
    }

    let nx = -dy / length;
    let ny = dx / length;
    let offsets: Vec<f64> = path
        .iter()
        .map(|point| ((point.0 - start.0) * nx) + ((point.1 - start.1) * ny))
        .collect();
    let extrema_threshold = channel_width_mm * 0.15;
    let mut turns = 0usize;
    for idx in 1..offsets.len() - 1 {
        let prev_delta = offsets[idx] - offsets[idx - 1];
        let next_delta = offsets[idx + 1] - offsets[idx];
        if prev_delta.abs() <= 1e-9 || next_delta.abs() <= 1e-9 {
            continue;
        }
        if prev_delta.signum() != next_delta.signum() && offsets[idx].abs() > extrema_threshold {
            turns += 1;
        }
    }

    let min_radius_mm = path
        .iter()
        .enumerate()
        .filter_map(|(idx, _)| estimate_local_bend_radius(path, idx))
        .filter(|radius| radius.is_finite())
        .fold(f64::INFINITY, f64::min);
    let bend_radius_mm = if min_radius_mm.is_finite() {
        min_radius_mm.max(channel_width_mm * 0.5)
    } else {
        channel_width_mm * 0.5
    };

    ChannelShape::Serpentine {
        segments: turns.saturating_add(1).max(2),
        bend_radius_m: bend_radius_mm * 1.0e-3,
        wave_type: crate::topology::SerpentineWaveType::Sine,
    }
}
