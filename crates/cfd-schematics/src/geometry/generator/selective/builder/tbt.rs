use super::{SelectiveTreeBuilder, SelectiveTreeRequest, NodeKind, JunctionFamily, ChannelVisualRole, VenturiGeometryMetadata};
use crate::domain::therapy_metadata::TherapyZone;

impl SelectiveTreeBuilder {
    pub(in super::super) fn build_tbt(
        &mut self,
        first_center_frac: f64,
        bi_treat_frac: f64,
        second_center_frac: f64,
        req: &SelectiveTreeRequest,
    ) {
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
}