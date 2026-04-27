use super::{
    ChannelVisualRole, JunctionFamily, NodeKind, SelectiveTreeBuilder, SelectiveTreeRequest,
    VenturiGeometryMetadata,
};
use crate::domain::therapy_metadata::TherapyZone;

impl SelectiveTreeBuilder {
    pub(in super::super) fn build_dtcv(
        &mut self,
        split1_center_frac: f64,
        split2_center_frac: f64,
        center_throat_count: u8,
        inter_spacing_m: f64,
        req: &SelectiveTreeRequest,
    ) {
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

            self.add_node(
                &m1_name,
                (96.0, m1_y),
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            );
            self.add_node(
                &m2_name,
                (104.0, m2_y),
                NodeKind::Junction,
                Some(JunctionFamily::Merge),
            );

            let start_pt = self
                .nodes
                .iter()
                .find(|n| n.id == self.node_ids[split])
                .expect("split node should exist")
                .point;

            self.add_channel(
                &format!("{prefix}_L"),
                split,
                &m1_name,
                vec![start_pt, (54.0, ys[0]), (90.0, ys[0]), (96.0, m1_y)],
                outer_leaf_w,
                req.branch_length_m,
                req.channel_height_m,
                ChannelVisualRole::PeripheralBypass,
                TherapyZone::HealthyBypass,
                None,
                None,
            );

            self.add_channel(
                &format!("{prefix}_C"),
                split,
                &m1_name,
                vec![start_pt, (54.0, ys[1]), (90.0, ys[1]), (96.0, m1_y)],
                outer_ctr_w,
                req.branch_length_m,
                req.channel_height_m,
                ChannelVisualRole::PeripheralBypass,
                TherapyZone::HealthyBypass,
                None,
                None,
            );

            self.add_channel(
                &format!("{prefix}_R"),
                split,
                &m2_name,
                vec![start_pt, (54.0, ys[2]), (90.0, ys[2]), (104.0, m2_y)],
                outer_leaf_w,
                req.branch_length_m,
                req.channel_height_m,
                ChannelVisualRole::PeripheralBypass,
                TherapyZone::HealthyBypass,
                None,
                None,
            );

            let m1_to_m2_w = (outer_leaf_w + outer_ctr_w).min(req.throat_width_m * 4.0);
            self.add_channel(
                &format!("{prefix}_m1_to_m2"),
                &m1_name,
                &m2_name,
                vec![(96.0, m1_y), (104.0, m2_y)],
                m1_to_m2_w,
                req.branch_length_m * 0.1,
                req.channel_height_m,
                ChannelVisualRole::PeripheralBypass,
                TherapyZone::HealthyBypass,
                None,
                None,
            );

            self.add_channel(
                &format!("{prefix}_to_bypass"),
                &m2_name,
                "periph_merge",
                vec![(104.0, m2_y), (104.0, y_mid)],
                side1,
                req.branch_length_m * 0.2,
                req.channel_height_m,
                ChannelVisualRole::PeripheralBypass,
                TherapyZone::HealthyBypass,
                None,
                None,
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
                    .expect("upstream node should exist")
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
                .expect("recovery node should exist")
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
