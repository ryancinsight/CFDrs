use super::super::types::{centerline_for_channel, Channel, ChannelSystem, ChannelType, Node, Point2D, SplitType};
use crate::config::{ChannelTypeConfig, GeometryConfig, TaperProfile};
use crate::domain::model::{ChannelShape, NodeKind};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::{
    ChannelVenturiSpec, ChannelVisualRole, JunctionFamily, JunctionGeometryMetadata,
    MetadataContainer, VenturiGeometryMetadata,
};
use std::collections::{HashMap, HashSet};

#[derive(Debug, Clone, Copy, PartialEq)]
pub struct CenterSerpentinePathSpec {
    pub segments: usize,
    pub bend_radius_m: f64,
}

#[derive(Debug, Clone, PartialEq)]
pub enum SelectiveTreeTopology {
    CascadeCenterTrifurcation {
        n_levels: usize,
        center_frac: f64,
        venturi_treatment_enabled: bool,
        center_serpentine: Option<CenterSerpentinePathSpec>,
    },
    IncrementalFiltrationTriBi {
        n_pretri: usize,
        pretri_center_frac: f64,
        terminal_tri_center_frac: f64,
        bi_treat_frac: f64,
        venturi_treatment_enabled: bool,
        center_serpentine: Option<CenterSerpentinePathSpec>,
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

pub fn create_selective_tree_geometry(request: &SelectiveTreeRequest) -> ChannelSystem {
    let mut builder = SelectiveTreeBuilder::new(request.box_dims_mm);
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
        } => builder.build_tbt(*first_center_frac, *bi_treat_frac, *second_center_frac, request),
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
    box_dims: (f64, f64),
    nodes: Vec<Node>,
    channels: Vec<Channel>,
    node_ids: HashMap<String, usize>,
}

impl SelectiveTreeBuilder {
    fn new(box_dims: (f64, f64)) -> Self {
        Self { box_dims, nodes: Vec::new(), channels: Vec::new(), node_ids: HashMap::new() }
    }

    fn finish(self) -> ChannelSystem {
        let (w, h) = self.box_dims;
        ChannelSystem {
            box_dims: self.box_dims,
            nodes: self.nodes,
            channels: self.channels,
            box_outline: vec![((0.0, 0.0), (w, 0.0)), ((w, 0.0), (w, h)), ((w, h), (0.0, h)), ((0.0, h), (0.0, 0.0))],
        }
    }

    fn add_node(&mut self, name: &str, point: Point2D, kind: NodeKind, family: Option<JunctionFamily>) {
        let junction_geometry = family.map(|junction_family| JunctionGeometryMetadata {
            junction_family,
            branch_angles_deg: Vec::new(),
            merge_angles_deg: Vec::new(),
        });
        let id = self.nodes.len();
        self.nodes.push(Node {
            id,
            name: Some(name.to_string()),
            point,
            kind: Some(kind),
            junction_geometry,
            metadata: None,
        });
        self.node_ids.insert(name.to_string(), id);
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
        zone: TherapyZone,
        shape: Option<ChannelShape>,
        venturi: Option<VenturiGeometryMetadata>,
    ) {
        let render_width_mm = (width_m * 1e3).clamp(0.20, 3.0);
        let path_len = path.len();
        let inlet_render = (venturi.as_ref().map_or(width_m, |v| v.inlet_width_m) * 1e3).clamp(0.20, 3.0);
        let throat_render = (venturi.as_ref().map_or(width_m, |v| v.throat_width_m) * 1e3).clamp(0.12, 2.0);
        let outlet_render = (venturi.as_ref().map_or(width_m, |v| v.outlet_width_m) * 1e3).clamp(0.20, 3.0);
        let channel_type = if let Some(venturi_geometry) = &venturi {
            let throat_idx = path_len / 2;
            let widths = (0..path_len)
                .map(|idx| if idx < throat_idx { inlet_render } else if idx == throat_idx { throat_render } else { outlet_render })
                .collect();
            ChannelType::Frustum {
                path,
                widths,
                inlet_width: inlet_render,
                throat_width: throat_render,
                outlet_width: outlet_render,
                taper_profile: crate::config::TaperProfile::Smooth,
                throat_position: 0.5,
                has_venturi_throat: venturi_geometry.throat_width_m > 0.0,
            }
        } else if matches!(shape, Some(ChannelShape::Serpentine { .. })) {
            ChannelType::Serpentine { path }
        } else if path_len > 2 {
            ChannelType::SmoothStraight { path }
        } else {
            ChannelType::Straight
        };
        let mut metadata = MetadataContainer::new();
        if let Some(venturi_geometry) = &venturi {
            metadata.insert(venturi_geometry.clone());
        }
        self.channels.push(Channel {
            id: self.channels.len(),
            name: Some(name.to_string()),
            from_node: self.node_ids[from],
            to_node: self.node_ids[to],
            width: render_width_mm,
            height: (height_m * 1e3).clamp(0.20, 3.0),
            channel_type,
            visual_role: Some(role),
            physical_length_m: Some(length_m),
            physical_width_m: Some(width_m),
            physical_height_m: Some(height_m),
            physical_shape: shape,
            therapy_zone: Some(zone),
            venturi_geometry: venturi,
            metadata: if metadata.is_empty() { None } else { Some(metadata) },
        });
    }

    fn y_mid(&self) -> f64 { self.box_dims.1 * 0.5 }
    fn leaf_triplet(&self, center_y: f64, gap: f64) -> [f64; 3] { [center_y - gap, center_y, center_y + gap] }
    fn serpentine_path(&self, start: Point2D, end: Point2D, segments: usize, amplitude_mm: f64) -> Vec<Point2D> {
        let mut path = vec![start];
        for idx in 1..segments {
            let t = idx as f64 / segments as f64;
            let x = start.0 + (end.0 - start.0) * t;
            let y = start.1 + if idx % 2 == 0 { amplitude_mm } else { -amplitude_mm };
            path.push((x, y));
        }
        path.push(end);
        path
    }

    fn build_cct(&mut self, n_levels: usize, center_frac: f64, venturi: bool, center_serp: Option<CenterSerpentinePathSpec>, req: &SelectiveTreeRequest) {
        let y_mid = self.y_mid();
        self.add_node("inlet", (0.0, y_mid), NodeKind::Inlet, None);
        self.add_node("periph_merge", (104.0, y_mid), NodeKind::Junction, Some(JunctionFamily::Merge));
        self.add_node("outlet_merge", (116.0, y_mid), NodeKind::Junction, Some(JunctionFamily::Merge));
        self.add_node("throat_in", (68.0, y_mid), NodeKind::Junction, None);
        self.add_node("throat_out", (76.0, y_mid), NodeKind::Junction, None);
        self.add_node("outlet", (self.box_dims.0, y_mid), NodeKind::Outlet, None);
        for lv in 0..n_levels { self.add_node(&format!("split_lv{lv}"), (18.0 + lv as f64 * 12.0, y_mid), NodeKind::Junction, Some(JunctionFamily::Trifurcation)); }
        self.add_channel("inlet_section", "inlet", "split_lv0", vec![(0.0, y_mid), (18.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::Trunk, TherapyZone::MixedFlow, None, None);
        let mut center_width = req.main_width_m;
        let side_frac = (1.0 - center_frac) * 0.5;
        for lv in 0..n_levels {
            let split = format!("split_lv{lv}");
            let next = if lv + 1 < n_levels { format!("split_lv{}", lv + 1) } else { "throat_in".to_string() };
            let off = 12.0 + 5.0 * (n_levels.saturating_sub(lv + 1)) as f64;
            let side_w = (center_width * side_frac).max(req.throat_width_m);
            let ctr_w = (center_width * center_frac).max(req.throat_width_m);
            for (id, y) in [("L", y_mid - off), ("R", y_mid + off)] {
                self.add_channel(&format!("{id}_lv{lv}"), &split, "periph_merge", vec![self.nodes[self.node_ids[&split]].point, (34.0 + lv as f64 * 10.0, y), (90.0, y), (104.0, y_mid)], side_w, req.branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None);
            }
            let start = self.nodes[self.node_ids[&split]].point;
            let end = self.nodes[self.node_ids[&next]].point;
            let path = if center_serp.is_some() { self.serpentine_path(start, end, center_serp.unwrap().segments.max(3), 1.6) } else { vec![start, end] };
            self.add_channel(&format!("center_lv{lv}"), &split, &next, path, ctr_w, req.branch_length_m, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, center_serp.map(|s| ChannelShape::Serpentine { segments: s.segments, bend_radius_m: s.bend_radius_m }), None);
            center_width = ctr_w;
        }
        let zone = if venturi { TherapyZone::CancerTarget } else { TherapyZone::CancerTarget };
        if venturi {
            self.add_channel("throat_section", "throat_in", "throat_out", vec![(68.0, y_mid), (72.0, y_mid), (76.0, y_mid)], req.throat_width_m, req.throat_length_m, req.channel_height_m, ChannelVisualRole::VenturiThroat, zone, None, Some(VenturiGeometryMetadata { throat_width_m: req.throat_width_m, throat_height_m: req.channel_height_m, throat_length_m: req.throat_length_m, inlet_width_m: center_width, outlet_width_m: center_width, convergent_half_angle_deg: 15.0, divergent_half_angle_deg: 15.0 }));
        } else {
            let path = if let Some(spec) = center_serp { self.serpentine_path((68.0, y_mid), (76.0, y_mid), spec.segments.max(3), 1.5) } else { vec![(68.0, y_mid), (76.0, y_mid)] };
            self.add_channel("treatment_section", "throat_in", "throat_out", path, center_width, req.branch_length_m.max(req.throat_length_m), req.channel_height_m, ChannelVisualRole::CenterTreatment, zone, center_serp.map(|s| ChannelShape::Serpentine { segments: s.segments, bend_radius_m: s.bend_radius_m }), None);
        }
        self.add_channel("center_to_merge", "throat_out", "outlet_merge", vec![(76.0, y_mid), (116.0, y_mid)], center_width, req.branch_length_m, req.channel_height_m, ChannelVisualRole::MergeCollector, TherapyZone::MixedFlow, None, None);
        self.add_channel("periph_to_merge", "periph_merge", "outlet_merge", vec![(104.0, y_mid), (116.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::MergeCollector, TherapyZone::HealthyBypass, None, None);
        self.add_channel("trunk_out", "outlet_merge", "outlet", vec![(116.0, y_mid), (self.box_dims.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::Trunk, TherapyZone::MixedFlow, None, None);
    }

    fn build_cif(&mut self, n_pretri: usize, pretri_center_frac: f64, terminal_tri_center_frac: f64, bi_treat_frac: f64, venturi: bool, center_serp: Option<CenterSerpentinePathSpec>, outlet_tail_length_m: f64, req: &SelectiveTreeRequest) {
        let y_mid = self.y_mid();
        for (name, x, kind, family) in [("inlet", 0.0, NodeKind::Inlet, None), ("periph_merge", 104.0, NodeKind::Junction, Some(JunctionFamily::Merge)), ("outlet_merge", 116.0, NodeKind::Junction, Some(JunctionFamily::Merge)), ("hy_tri", 52.0, NodeKind::Junction, Some(JunctionFamily::Trifurcation)), ("hy_bi", 64.0, NodeKind::Junction, Some(JunctionFamily::Bifurcation)), ("throat_in", 72.0, NodeKind::Junction, None), ("throat_out", 80.0, NodeKind::Junction, None), ("outlet", self.box_dims.0, NodeKind::Outlet, None)] { self.add_node(name, (x, y_mid), kind, family); }
        for lv in 0..n_pretri { self.add_node(&format!("split_lv{lv}"), (18.0 + lv as f64 * 11.0, y_mid), NodeKind::Junction, Some(JunctionFamily::Trifurcation)); }
        self.add_channel("inlet_section", "inlet", "split_lv0", vec![(0.0, y_mid), (18.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::Trunk, TherapyZone::MixedFlow, None, None);
        let mut center_width = req.main_width_m;
        for lv in 0..n_pretri {
            let split = format!("split_lv{lv}");
            let next = if lv + 1 < n_pretri { format!("split_lv{}", lv + 1) } else { "hy_tri".to_string() };
            let side_off = 10.0 + 4.0 * (n_pretri.saturating_sub(lv + 1)) as f64;
            let side_w = (center_width * (1.0 - pretri_center_frac) * 0.5).max(req.throat_width_m);
            let ctr_w = (center_width * pretri_center_frac).max(req.throat_width_m);
            for (id, y) in [("L", y_mid - side_off), ("R", y_mid + side_off)] {
                self.add_channel(&format!("{id}_lv{lv}"), &split, "periph_merge", vec![self.nodes[self.node_ids[&split]].point, (34.0 + lv as f64 * 9.0, y), (50.0 + lv as f64 * 8.0, y), (92.0, y), (104.0, y_mid)], side_w, req.branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None);
            }
            let start = self.nodes[self.node_ids[&split]].point;
            let end = self.nodes[self.node_ids[&next]].point;
            let path = if let Some(spec) = center_serp { self.serpentine_path(start, end, spec.segments.max(3), 1.4) } else { vec![start, end] };
            self.add_channel(&format!("center_lv{lv}"), &split, &next, path, ctr_w, req.branch_length_m, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, center_serp.map(|s| ChannelShape::Serpentine { segments: s.segments, bend_radius_m: s.bend_radius_m }), None);
            center_width = ctr_w;
        }
        let tri_side = (center_width * (1.0 - terminal_tri_center_frac) * 0.5).max(req.throat_width_m);
        let tri_ctr = (center_width * terminal_tri_center_frac).max(req.throat_width_m);
        for (id, y) in [("hy_tri_L", y_mid - 8.0), ("hy_tri_R", y_mid + 8.0)] {
            self.add_channel(id, "hy_tri", "periph_merge", vec![(52.0, y_mid), (62.0, y), (92.0, y), (104.0, y_mid)], tri_side, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None);
        }
        let tri_path = if let Some(spec) = center_serp { self.serpentine_path((52.0, y_mid), (64.0, y_mid), spec.segments.max(3), 1.3) } else { vec![(52.0, y_mid), (64.0, y_mid)] };
        self.add_channel("hy_tri_center", "hy_tri", "hy_bi", tri_path, tri_ctr, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, center_serp.map(|s| ChannelShape::Serpentine { segments: s.segments, bend_radius_m: s.bend_radius_m }), None);
        let treat_w = (tri_ctr * bi_treat_frac).max(req.throat_width_m);
        let bypass_w = (tri_ctr * (1.0 - bi_treat_frac)).max(req.throat_width_m);
        self.add_channel("hy_bi_bypass", "hy_bi", "periph_merge", vec![(64.0, y_mid), (76.0, y_mid + 10.0), (92.0, y_mid + 10.0), (104.0, y_mid)], bypass_w, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None);
        let bi_path = if let Some(spec) = center_serp { self.serpentine_path((64.0, y_mid), (72.0, y_mid), spec.segments.max(3), 1.2) } else { vec![(64.0, y_mid), (72.0, y_mid)] };
        self.add_channel("hy_bi_treat", "hy_bi", "throat_in", bi_path, treat_w, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, center_serp.map(|s| ChannelShape::Serpentine { segments: s.segments, bend_radius_m: s.bend_radius_m }), None);
        if venturi {
            self.add_channel("throat_section", "throat_in", "throat_out", vec![(72.0, y_mid), (76.0, y_mid), (80.0, y_mid)], req.throat_width_m, req.throat_length_m, req.channel_height_m, ChannelVisualRole::VenturiThroat, TherapyZone::CancerTarget, None, Some(VenturiGeometryMetadata { throat_width_m: req.throat_width_m, throat_height_m: req.channel_height_m, throat_length_m: req.throat_length_m, inlet_width_m: treat_w, outlet_width_m: treat_w, convergent_half_angle_deg: 15.0, divergent_half_angle_deg: 15.0 }));
        } else {
            let path = if let Some(spec) = center_serp { self.serpentine_path((72.0, y_mid), (80.0, y_mid), spec.segments.max(3), 1.2) } else { vec![(72.0, y_mid), (80.0, y_mid)] };
            self.add_channel("treatment_section", "throat_in", "throat_out", path, treat_w, req.hybrid_branch_length_m.max(req.throat_length_m), req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, center_serp.map(|s| ChannelShape::Serpentine { segments: s.segments, bend_radius_m: s.bend_radius_m }), None);
        }
        self.add_channel("center_to_merge", "throat_out", "outlet_merge", vec![(80.0, y_mid), (116.0, y_mid)], treat_w, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::MergeCollector, TherapyZone::MixedFlow, None, None);
        self.add_channel("periph_to_merge", "periph_merge", "outlet_merge", vec![(104.0, y_mid), (116.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::MergeCollector, TherapyZone::HealthyBypass, None, None);
        self.add_channel("trunk_out", "outlet_merge", "outlet", vec![(116.0, y_mid), (self.box_dims.0, y_mid)], req.main_width_m, outlet_tail_length_m, req.channel_height_m, ChannelVisualRole::Trunk, TherapyZone::MixedFlow, None, None);
    }

    fn build_tbt(&mut self, first_center_frac: f64, bi_treat_frac: f64, second_center_frac: f64, req: &SelectiveTreeRequest) {
        let y_mid = self.y_mid();
        for (name, x, kind, family) in [("inlet", 0.0, NodeKind::Inlet, None), ("split_lv0", 24.0, NodeKind::Junction, Some(JunctionFamily::Trifurcation)), ("split_bi", 46.0, NodeKind::Junction, Some(JunctionFamily::Bifurcation)), ("split_lv1", 70.0, NodeKind::Junction, Some(JunctionFamily::Trifurcation)), ("throat_in", 86.0, NodeKind::Junction, None), ("throat_out", 92.0, NodeKind::Junction, None), ("periph_merge", 104.0, NodeKind::Junction, Some(JunctionFamily::Merge)), ("outlet_merge", 116.0, NodeKind::Junction, Some(JunctionFamily::Merge)), ("outlet", self.box_dims.0, NodeKind::Outlet, None)] { self.add_node(name, (x, y_mid), kind, family); }
        self.add_channel("inlet_section", "inlet", "split_lv0", vec![(0.0, y_mid), (24.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::Trunk, TherapyZone::MixedFlow, None, None);
        let side1 = (req.main_width_m * (1.0 - first_center_frac) * 0.5).max(req.throat_width_m);
        let ctr1 = (req.main_width_m * first_center_frac).max(req.throat_width_m);
        for (id, y) in [("L_lv0", y_mid - 18.0), ("R_lv0", y_mid + 18.0)] { self.add_channel(id, "split_lv0", "periph_merge", vec![(24.0, y_mid), (40.0, y), (96.0, y), (104.0, y_mid)], side1, req.branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None); }
        let bi_bypass_w = (ctr1 * (1.0 - bi_treat_frac)).max(req.throat_width_m);
        let bi_treat_w = (ctr1 * bi_treat_frac).max(req.throat_width_m);
        self.add_channel("center_lv0", "split_lv0", "split_bi", vec![(24.0, y_mid), (46.0, y_mid)], ctr1, req.branch_length_m, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, None, None);
        self.add_channel("bi_bypass", "split_bi", "periph_merge", vec![(46.0, y_mid), (60.0, y_mid + 10.0), (96.0, y_mid + 10.0), (104.0, y_mid)], bi_bypass_w, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None);
        self.add_channel("center_bi", "split_bi", "split_lv1", vec![(46.0, y_mid), (70.0, y_mid)], bi_treat_w, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, None, None);
        let side2 = (bi_treat_w * (1.0 - second_center_frac) * 0.5).max(req.throat_width_m);
        let ctr2 = (bi_treat_w * second_center_frac).max(req.throat_width_m);
        for (id, y) in [("L_lv1", y_mid - 10.0), ("R_lv1", y_mid + 10.0)] { self.add_channel(id, "split_lv1", "periph_merge", vec![(70.0, y_mid), (82.0, y), (96.0, y), (104.0, y_mid)], side2, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None); }
        self.add_channel("throat_section", "throat_in", "throat_out", vec![(86.0, y_mid), (89.0, y_mid), (92.0, y_mid)], req.throat_width_m, req.throat_length_m, req.channel_height_m, ChannelVisualRole::VenturiThroat, TherapyZone::CancerTarget, None, Some(VenturiGeometryMetadata { throat_width_m: req.throat_width_m, throat_height_m: req.channel_height_m, throat_length_m: req.throat_length_m, inlet_width_m: ctr2, outlet_width_m: ctr2, convergent_half_angle_deg: 15.0, divergent_half_angle_deg: 15.0 }));
        self.add_channel("center_lv1", "split_lv1", "throat_in", vec![(70.0, y_mid), (86.0, y_mid)], ctr2, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, None, None);
        self.add_channel("center_to_merge", "throat_out", "outlet_merge", vec![(92.0, y_mid), (116.0, y_mid)], ctr2, req.hybrid_branch_length_m, req.channel_height_m, ChannelVisualRole::MergeCollector, TherapyZone::MixedFlow, None, None);
        self.add_channel("periph_to_merge", "periph_merge", "outlet_merge", vec![(104.0, y_mid), (116.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::MergeCollector, TherapyZone::HealthyBypass, None, None);
        self.add_channel("trunk_out", "outlet_merge", "outlet", vec![(116.0, y_mid), (self.box_dims.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::Trunk, TherapyZone::MixedFlow, None, None);
    }

fn build_dtcv(&mut self, split1_center_frac: f64, split2_center_frac: f64, center_throat_count: u8, inter_spacing_m: f64, req: &SelectiveTreeRequest) {
        let y_mid = self.y_mid();
        for (name, x, y, kind, family) in [("inlet", 0.0, y_mid, NodeKind::Inlet, None), ("split1", 18.0, y_mid, NodeKind::Junction, Some(JunctionFamily::Trifurcation)), ("split2_upper", 36.0, y_mid - 21.0, NodeKind::Junction, Some(JunctionFamily::Trifurcation)), ("split2", 36.0, y_mid, NodeKind::Junction, Some(JunctionFamily::Trifurcation)), ("split2_lower", 36.0, y_mid + 21.0, NodeKind::Junction, Some(JunctionFamily::Trifurcation)), ("merge_bypass", 100.0, y_mid, NodeKind::Junction, Some(JunctionFamily::Merge)), ("merge_center", 100.0, y_mid, NodeKind::Junction, Some(JunctionFamily::Merge)), ("merge_out", 118.0, y_mid, NodeKind::Junction, Some(JunctionFamily::Merge)), ("outlet", self.box_dims.0, y_mid, NodeKind::Outlet, None)] { self.add_node(name, (x, y), kind, family); }
        self.add_channel("inlet_section", "inlet", "split1", vec![(0.0, y_mid), (18.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::Trunk, TherapyZone::MixedFlow, None, None);
        let side1 = (req.main_width_m * (1.0 - split1_center_frac) * 0.5).max(req.throat_width_m);
        let ctr1 = (req.main_width_m * split1_center_frac).max(req.throat_width_m);
        for (id, to, y) in [("upper_trunk", "split2_upper", y_mid - 21.0), ("center_trunk", "split2", y_mid), ("lower_trunk", "split2_lower", y_mid + 21.0)] {
            let zone = if id == "center_trunk" { TherapyZone::CancerTarget } else { TherapyZone::HealthyBypass };
            let width = if id == "center_trunk" { ctr1 } else { side1 };
            self.add_channel(id, "split1", to, vec![(18.0, y_mid), (36.0, y)], width, req.branch_length_m, req.channel_height_m, if id == "center_trunk" { ChannelVisualRole::CenterTreatment } else { ChannelVisualRole::PeripheralBypass }, zone, None, None);
        }
        let upper = self.leaf_triplet(y_mid - 21.0, 4.0);
        let center = self.leaf_triplet(y_mid, 7.0);
        let lower = self.leaf_triplet(y_mid + 21.0, 4.0);
        let outer_leaf_w = (side1 * (1.0 - split2_center_frac) * 0.5).max(req.throat_width_m);
        let outer_ctr_w = (side1 * split2_center_frac).max(req.throat_width_m);
        for ((prefix, split), ys) in [(("upper", "split2_upper"), upper), (("lower", "split2_lower"), lower)] {
            for (suffix, y, width) in [("L", ys[0], outer_leaf_w), ("C", ys[1], outer_ctr_w), ("R", ys[2], outer_leaf_w)] {
                self.add_channel(&format!("{prefix}_{suffix}"), split, "merge_bypass", vec![self.nodes[self.node_ids[split]].point, (54.0, y), (90.0, y), (100.0, y_mid)], width, req.branch_length_m, req.channel_height_m, ChannelVisualRole::PeripheralBypass, TherapyZone::HealthyBypass, None, None);
            }
        }
        let center_side_w = (ctr1 * (1.0 - split2_center_frac) * 0.5).max(req.throat_width_m);
        let center_ctr_w = (ctr1 * split2_center_frac).max(req.throat_width_m);
        for (suffix, y, width) in [("L", center[0], center_side_w), ("C", center[1], center_ctr_w), ("R", center[2], center_side_w)] {
            let sub = format!("center_{suffix}");
            let sub_in = format!("{sub}_in");
            self.add_node(&sub_in, (54.0, y), NodeKind::Junction, None);
            self.add_channel(&format!("{sub}_approach"), "split2", &sub_in, vec![(36.0, y_mid), (54.0, y)], width, req.branch_length_m * 0.35, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, None, None);
            if center_throat_count == 0 {
                self.add_channel(&format!("{sub}_treatment"), &sub_in, "merge_center", vec![(54.0, y), (90.0, y), (100.0, y_mid)], width, req.branch_length_m, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, None, None);
                continue;
            }
            let mut prev = sub_in.clone();
            for k in 0..center_throat_count {
                let tin = format!("{sub}_t{k}_in");
                let tout = format!("{sub}_t{k}_out");
                let x0 = 60.0 + f64::from(k) * 8.0;
                self.add_node(&tin, (x0, y), NodeKind::Junction, None);
                self.add_node(&tout, (x0 + 4.0, y), NodeKind::Junction, None);
                self.add_channel(&format!("{sub}_conv{k}"), &prev, &tin, vec![self.nodes[self.node_ids[&prev]].point, (x0, y)], width, req.branch_length_m * 0.12, req.channel_height_m, ChannelVisualRole::CenterTreatment, TherapyZone::CancerTarget, None, None);
                self.add_channel(&format!("{sub}_throat{k}"), &tin, &tout, vec![(x0, y), (x0 + 2.0, y), (x0 + 4.0, y)], req.throat_width_m, req.throat_length_m, req.channel_height_m, ChannelVisualRole::VenturiThroat, TherapyZone::CancerTarget, None, Some(VenturiGeometryMetadata { throat_width_m: req.throat_width_m, throat_height_m: req.channel_height_m, throat_length_m: req.throat_length_m, inlet_width_m: width, outlet_width_m: width, convergent_half_angle_deg: 15.0, divergent_half_angle_deg: 15.0 }));
                prev = tout;
                if k + 1 < center_throat_count {
                    let rdv = format!("{sub}_rdv{k}");
                    self.add_node(&rdv, (x0 + 6.0, y), NodeKind::Junction, None);
                    self.add_channel(&format!("{sub}_diffuser{k}"), &prev, &rdv, vec![(x0 + 4.0, y), (x0 + 6.0, y)], width, inter_spacing_m, req.channel_height_m, ChannelVisualRole::Diffuser, TherapyZone::CancerTarget, None, None);
                    prev = rdv;
                }
            }
            self.add_channel(&format!("{sub}_recovery"), &prev, "merge_center", vec![self.nodes[self.node_ids[&prev]].point, (90.0, y), (100.0, y_mid)], width, req.branch_length_m * 0.18, req.channel_height_m, ChannelVisualRole::Diffuser, TherapyZone::CancerTarget, None, None);
        }
        self.add_channel("bypass_to_merge", "merge_bypass", "merge_out", vec![(100.0, y_mid), (118.0, y_mid)], req.main_width_m * (1.0 - split1_center_frac), req.trunk_length_m, req.channel_height_m, ChannelVisualRole::MergeCollector, TherapyZone::HealthyBypass, None, None);
        self.add_channel("center_to_merge", "merge_center", "merge_out", vec![(100.0, y_mid), (118.0, y_mid)], ctr1, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::MergeCollector, TherapyZone::MixedFlow, None, None);
        self.add_channel("outlet_section", "merge_out", "outlet", vec![(118.0, y_mid), (self.box_dims.0, y_mid)], req.main_width_m, req.trunk_length_m, req.channel_height_m, ChannelVisualRole::Trunk, TherapyZone::MixedFlow, None, None);
    }
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum PrimitiveSelectiveSplitKind {
    Bi,
    Tri,
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
) -> ChannelSystem {
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
        })
        .collect();

    let geometry_config = GeometryConfig {
        wall_clearance: 1.0,
        channel_width: (request.main_width_m * 1.0e3).clamp(0.2, 12.0),
        channel_height: (request.channel_height_m * 1.0e3).clamp(0.2, 5.0),
        ..Default::default()
    };

    let mut system = super::create_geometry(
        request.box_dims_mm,
        &splits,
        &geometry_config,
        &ChannelTypeConfig::AllStraight,
    );
    annotate_primitive_tree(&mut system, request);
    system
}

fn annotate_primitive_tree(system: &mut ChannelSystem, request: &PrimitiveSelectiveTreeRequest) {
    if system.nodes.is_empty() || system.channels.is_empty() {
        return;
    }

    let mid_x = system.box_dims.0 * 0.5;
    let mid_y = system.box_dims.1 * 0.5;
    let inlet_node = system
        .nodes
        .iter()
        .min_by(|a, b| a.point.0.total_cmp(&b.point.0))
        .map(|node| node.id);
    let outlet_node = system
        .nodes
        .iter()
        .max_by(|a, b| a.point.0.total_cmp(&b.point.0))
        .map(|node| node.id);
    let (Some(inlet_node), Some(outlet_node)) = (inlet_node, outlet_node) else {
        return;
    };

    let mut adjacency: HashMap<usize, Vec<(usize, usize)>> = HashMap::new();
    for (idx, channel) in system.channels.iter().enumerate() {
        adjacency
            .entry(channel.from_node)
            .or_default()
            .push((channel.to_node, idx));
        adjacency
            .entry(channel.to_node)
            .or_default()
            .push((channel.from_node, idx));
    }

    let leaf_nodes: Vec<usize> = system
        .nodes
        .iter()
        .filter(|node| (node.point.0 - mid_x).abs() < 1e-6)
        .map(|node| node.id)
        .collect();
    let mut ordered_leaves = leaf_nodes;
    ordered_leaves.sort_by(|a, b| system.nodes[*a].point.1.total_cmp(&system.nodes[*b].point.1));

    let treatment_leaf_indices = primitive_treatment_leaf_indices(
        ordered_leaves.len(),
        request.split_sequence.first().copied(),
    );
    let treatment_leaves: HashSet<usize> = treatment_leaf_indices
        .into_iter()
        .filter_map(|idx| ordered_leaves.get(idx).copied())
        .collect();

    let mut treatment_channels = HashSet::new();
    for &leaf in &treatment_leaves {
        if let Some(path) = channel_path_between(
            inlet_node,
            leaf,
            system,
            &adjacency,
            Some(mid_x),
            true,
        ) {
            treatment_channels.extend(path);
        }
        if let Some(path) = channel_path_between(
            leaf,
            outlet_node,
            system,
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
        request.main_width_m
    };

    for node in &mut system.nodes {
        let degree = adjacency.get(&node.id).map_or(0, Vec::len);
        node.kind = Some(if node.id == inlet_node {
            NodeKind::Inlet
        } else if node.id == outlet_node {
            NodeKind::Outlet
        } else if degree > 1 {
            NodeKind::Junction
        } else {
            NodeKind::Junction
        });
        if node.id == inlet_node {
            node.name = Some("inlet".to_string());
        } else if node.id == outlet_node {
            node.name = Some("outlet".to_string());
        } else {
            node.name = Some(format!("jn{}", node.id));
        }
    }

    let mut inlet_named = false;
    let mut outlet_named = false;
    let mut venturi_named = 0usize;
    for (idx, channel) in system.channels.iter_mut().enumerate() {
        let points = centerline_for_channel(channel, &system.nodes);
        let min_x = points
            .iter()
            .map(|(x, _)| *x)
            .fold(f64::INFINITY, f64::min);
        let max_x = points
            .iter()
            .map(|(x, _)| *x)
            .fold(f64::NEG_INFINITY, f64::max);
        let avg_y = if points.is_empty() {
            mid_y
        } else {
            points.iter().map(|(_, y)| *y).sum::<f64>() / points.len() as f64
        };

        let is_treatment = treatment_channels.contains(&idx);
        let touches_inlet = channel.from_node == inlet_node || channel.to_node == inlet_node;
        let touches_outlet = channel.from_node == outlet_node || channel.to_node == outlet_node;
        let starts_at_treatment_leaf = treatment_leaves.contains(&channel.from_node);
        let is_treatment_window_channel =
            starts_at_treatment_leaf && max_x > mid_x + 1e-6 && min_x >= mid_x - 1e-6;
        let is_trunk = touches_inlet || touches_outlet || (avg_y - mid_y).abs() < 1e-6;

        let physical_width = if touches_inlet || touches_outlet {
            request.main_width_m
        } else if is_treatment {
            treatment_width
        } else {
            bypass_width
        };
        channel.width = (physical_width * 1.0e3).clamp(0.15, 3.5);
        channel.height = (request.channel_height_m * 1.0e3).clamp(0.15, 3.0);
        channel.physical_width_m = Some(physical_width);
        channel.physical_height_m = Some(request.channel_height_m);
        channel.physical_length_m = Some(polyline_length_mm(&points) * 1.0e-3);

        let mut role = if is_treatment {
            ChannelVisualRole::CenterTreatment
        } else {
            ChannelVisualRole::PeripheralBypass
        };
        if touches_inlet || touches_outlet {
            role = ChannelVisualRole::Trunk;
        } else if max_x > mid_x {
            role = ChannelVisualRole::MergeCollector;
        }
        channel.visual_role = Some(role);
        channel.therapy_zone = Some(if is_treatment {
            TherapyZone::CancerTarget
        } else if is_trunk {
            TherapyZone::MixedFlow
        } else {
            TherapyZone::HealthyBypass
        });

        if touches_inlet && !inlet_named {
            channel.name = Some("inlet_section".to_string());
            inlet_named = true;
        } else if touches_outlet && !outlet_named {
            channel.name = Some("trunk_out".to_string());
            outlet_named = true;
        } else {
            channel.name = Some(format!("seg_{}", idx));
        }

        if is_treatment
            && is_treatment_window_channel
            && request.center_serpentine.is_some()
            && !request.treatment_branch_venturi_enabled
        {
            if let Some(spec) = request.center_serpentine {
                channel.channel_type = ChannelType::Serpentine {
                    path: serpentine_overlay_path(&points, spec.segments.max(3), 1.5),
                };
                channel.physical_shape = Some(ChannelShape::Serpentine {
                    segments: spec.segments,
                    bend_radius_m: spec.bend_radius_m,
                });
            }
        }

        if is_treatment && is_treatment_window_channel && request.treatment_branch_venturi_enabled {
            venturi_named += 1;
            let path = simple_path(points.first().copied(), points.last().copied(), avg_y);
            channel.channel_type = ChannelType::Frustum {
                path,
                widths: vec![
                    (physical_width * 1.0e3).clamp(0.15, 3.5),
                    (request.throat_width_m * 1.0e3).clamp(0.12, 2.0),
                    (physical_width * 1.0e3).clamp(0.15, 3.5),
                ],
                inlet_width: (physical_width * 1.0e3).clamp(0.15, 3.5),
                throat_width: (request.throat_width_m * 1.0e3).clamp(0.12, 2.0),
                outlet_width: (physical_width * 1.0e3).clamp(0.15, 3.5),
                taper_profile: TaperProfile::Smooth,
                throat_position: 0.5,
                has_venturi_throat: true,
            };
            channel.visual_role = Some(ChannelVisualRole::VenturiThroat);
            channel.venturi_geometry = Some(VenturiGeometryMetadata {
                throat_width_m: request.throat_width_m,
                throat_height_m: request.channel_height_m,
                throat_length_m: request.throat_length_m,
                inlet_width_m: physical_width,
                outlet_width_m: physical_width,
                convergent_half_angle_deg: 15.0,
                divergent_half_angle_deg: 15.0,
            });
            let metadata = channel
                .metadata
                .get_or_insert_with(MetadataContainer::new);
            metadata.insert(ChannelVenturiSpec {
                n_throats: request.treatment_branch_throat_count.max(1),
                is_ctc_stream: true,
                throat_width_m: request.throat_width_m,
                height_m: request.channel_height_m,
                inter_throat_spacing_m: request.throat_length_m * 2.0,
            });
            channel.name = Some(format!("throat_section_{venturi_named}"));
        } else if !is_treatment {
            let metadata = channel
                .metadata
                .get_or_insert_with(MetadataContainer::new);
            metadata.insert(ChannelVenturiSpec {
                n_throats: 0,
                is_ctc_stream: false,
                throat_width_m: request.throat_width_m,
                height_m: request.channel_height_m,
                inter_throat_spacing_m: request.throat_length_m * 2.0,
            });
        }
    }
}

fn primitive_treatment_leaf_indices(
    leaf_count: usize,
    first_stage: Option<PrimitiveSelectiveSplitKind>,
) -> Vec<usize> {
    if leaf_count == 0 {
        return Vec::new();
    }
    match first_stage {
        Some(PrimitiveSelectiveSplitKind::Tri) => {
            if leaf_count == 3 {
                vec![1]
            } else {
                let block = leaf_count / 3;
                (block..(2 * block)).collect()
            }
        }
        _ => (0..leaf_count).collect(),
    }
}

fn channel_path_between(
    start: usize,
    goal: usize,
    system: &ChannelSystem,
    adjacency: &HashMap<usize, Vec<(usize, usize)>>,
    mid_x: Option<f64>,
    left_half: bool,
) -> Option<Vec<usize>> {
    let mut visited = HashSet::new();
    let mut path = Vec::new();
    if dfs_channel_path(
        start,
        goal,
        system,
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
    current: usize,
    goal: usize,
    system: &ChannelSystem,
    adjacency: &HashMap<usize, Vec<(usize, usize)>>,
    mid_x: Option<f64>,
    left_half: bool,
    visited: &mut HashSet<usize>,
    path: &mut Vec<usize>,
) -> bool {
    if current == goal {
        return true;
    }
    visited.insert(current);
    let Some(neighbors) = adjacency.get(&current) else {
        return false;
    };
    for &(next, channel_idx) in neighbors {
        if visited.contains(&next) {
            continue;
        }
        if let Some(mid_x) = mid_x {
            let from = system.nodes[current].point.0;
            let to = system.nodes[next].point.0;
            let keep = if left_half {
                from <= mid_x + 1e-6 && to <= mid_x + 1e-6
            } else {
                from >= mid_x - 1e-6 && to >= mid_x - 1e-6
            };
            if !keep {
                continue;
            }
        }
        path.push(channel_idx);
        if dfs_channel_path(next, goal, system, adjacency, mid_x, left_half, visited, path) {
            return true;
        }
        path.pop();
    }
    visited.remove(&current);
    false
}

fn serpentine_overlay_path(points: &[Point2D], segments: usize, amplitude_mm: f64) -> Vec<Point2D> {
    let Some(start) = points.first().copied() else {
        return Vec::new();
    };
    let Some(end) = points.last().copied() else {
        return vec![start];
    };
    let mut path = vec![start];
    for idx in 1..segments {
        let t = idx as f64 / segments as f64;
        let x = start.0 + (end.0 - start.0) * t;
        let y = start.1 + if idx % 2 == 0 { amplitude_mm } else { -amplitude_mm };
        path.push((x, y));
    }
    path.push(end);
    path
}

fn simple_path(start: Option<Point2D>, end: Option<Point2D>, y_fallback: f64) -> Vec<Point2D> {
    match (start, end) {
        (Some(a), Some(b)) => {
            let mid_x = (a.0 + b.0) * 0.5;
            vec![a, (mid_x, y_fallback), b]
        }
        (Some(a), None) => vec![a],
        (None, Some(b)) => vec![b],
        (None, None) => Vec::new(),
    }
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
