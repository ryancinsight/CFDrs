use super::{
    CenterSerpentinePathSpec, ChannelShape, ChannelVisualRole, JunctionFamily, NodeKind,
    SelectiveTreeBuilder, SelectiveTreeRequest, VenturiGeometryMetadata,
};
use crate::domain::therapy_metadata::TherapyZone;

impl SelectiveTreeBuilder {
    pub(in super::super) fn build_cct(
        &mut self,
        n_levels: usize,
        center_frac: f64,
        venturi: bool,
        center_serp: Option<CenterSerpentinePathSpec>,
        req: &SelectiveTreeRequest,
    ) {
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
                    .expect("split node should exist")
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
                .expect("split node should exist")
                .point;
            let end = self
                .nodes
                .iter()
                .find(|n| n.id == self.node_ids[&next])
                .expect("next node should exist")
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
}
