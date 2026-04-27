use super::{
    CenterSerpentinePathSpec, ChannelShape, ChannelVisualRole, JunctionFamily, NodeKind,
    SelectiveTreeBuilder, SelectiveTreeRequest, VenturiGeometryMetadata,
};
use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::builders::ChannelExt;
use crate::geometry::metadata::IncrementalFiltrationParams;

impl SelectiveTreeBuilder {
    pub(in super::super) fn build_cif(
        &mut self,
        n_pretri: usize,
        pretri_center_frac: f64,
        terminal_tri_center_frac: f64,
        bi_treat_frac: f64,
        venturi: bool,
        center_serp: Option<CenterSerpentinePathSpec>,
        outlet_tail_length_m: f64,
        req: &SelectiveTreeRequest,
    ) {
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

        if let Some(channel) = self.channels.last_mut() {
            channel.add_metadata(IncrementalFiltrationParams {
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
                    .expect("split node should exist")
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
                .expect("split node should exist")
                .point;
            let end = self
                .nodes
                .iter()
                .find(|n| n.id == self.node_ids[&next])
                .expect("next node should exist")
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
}
