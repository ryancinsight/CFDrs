//! Specialized composite presets: cell separation, asymmetric, constriction,
//! spiral, parallel microchannel array, and selective tree topologies.

mod basic;
mod filtration;
mod parallel_lane;
mod trifurcation;

pub use basic::{
    asymmetric_bifurcation_serpentine_rect, cell_separation_rect,
    constriction_expansion_array_rect, parallel_microchannel_array_rect,
    primitive_selective_split_tree_rect, spiral_channel_rect,
};
pub use filtration::{
    cascade_center_trifurcation_rect, incremental_filtration_tri_bi_rect,
    incremental_filtration_tri_bi_rect_staged, incremental_filtration_tri_bi_rect_staged_remerge,
};
pub use parallel_lane::CenterSerpentineSpec;
pub use trifurcation::{
    asymmetric_trifurcation_venturi_rect, cascade_tri_bi_tri_selective_rect,
    double_trifurcation_cif_venturi_rect,
};

// Sub-modules provide their own imports

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::model::{ChannelShape, CrossSectionSpec, NetworkBlueprint};
    use crate::geometry::metadata::IncrementalFiltrationParams;

    fn rect_width(bp: &NetworkBlueprint, id: &str) -> f64 {
        let ch = bp
            .channels
            .iter()
            .find(|c| c.id.as_str() == id)
            .expect("channel must exist");
        match ch.cross_section {
            CrossSectionSpec::Rectangular { width_m, .. } => width_m,
            CrossSectionSpec::Circular { .. } => panic!("expected rectangular cross-section"),
        }
    }

    fn channel_length(bp: &NetworkBlueprint, id: &str) -> f64 {
        bp.channels
            .iter()
            .find(|c| c.id.as_str() == id)
            .map_or_else(|| panic!("channel {id} must exist"), |c| c.length_m)
    }

    #[test]
    fn staged_cif_terminal_fraction_controls_hybrid_center_width() {
        let bp = incremental_filtration_tri_bi_rect_staged(
            "cif-staged",
            12e-3,
            8e-3,
            6e-3,
            2,
            2.0e-3,
            0.33,
            0.55,
            0.68,
            100e-6,
            300e-6,
            1.0e-3,
            true,
            None,
        );

        let center_lv1 = rect_width(&bp, "center_lv1");
        let hy_tri_center = rect_width(&bp, "hy_tri_center");
        let expected = center_lv1 * 0.55;

        assert!(
            (hy_tri_center - expected).abs() < 1e-12,
            "terminal staged fraction must set hy_tri_center width"
        );
    }

    #[test]
    fn staged_cif_metadata_encodes_all_stage_parameters() {
        let bp = incremental_filtration_tri_bi_rect_staged(
            "cif-meta", 12e-3, 8e-3, 6e-3, 3, 2.0e-3, 0.45, 0.55, 0.76, 100e-6, 300e-6, 1.0e-3,
            true, None,
        );
        let inlet = bp
            .channels
            .iter()
            .find(|c| c.id.as_str() == "inlet_section")
            .expect("inlet_section must exist");
        let params = inlet
            .metadata
            .as_ref()
            .and_then(
                crate::geometry::metadata::MetadataContainer::get::<IncrementalFiltrationParams>,
            )
            .expect("IncrementalFiltrationParams metadata must exist");

        assert_eq!(params.n_pretri, 3);
        assert!((params.pretri_center_frac - 0.45).abs() < 1e-12);
        assert!((params.terminal_tri_center_frac - 0.55).abs() < 1e-12);
        assert!((params.bi_treat_frac - 0.76).abs() < 1e-12);
        assert!((params.outlet_tail_length_m - 12e-3).abs() < 1e-12);
    }

    #[test]
    fn staged_cif_remerge_tail_controls_trunk_out_length() {
        let bp = incremental_filtration_tri_bi_rect_staged_remerge(
            "cif-remerge",
            12e-3,
            8e-3,
            6e-3,
            2,
            2.0e-3,
            0.45,
            0.55,
            0.68,
            100e-6,
            300e-6,
            1.5e-3,
            1.0e-3,
            true,
            None,
        );

        assert!((channel_length(&bp, "periph_to_merge") - 12e-3).abs() < 1e-12);
        assert!((channel_length(&bp, "trunk_out") - 1.5e-3).abs() < 1e-12);
    }

    #[test]
    fn staged_cif_pretri_levels_ramp_center_bias() {
        let bp = incremental_filtration_tri_bi_rect_staged(
            "cif-stage-ramp",
            12e-3,
            8e-3,
            6e-3,
            3,
            2.0e-3,
            0.45,
            0.60,
            0.68,
            100e-6,
            300e-6,
            1.0e-3,
            true,
            None,
        );

        let c0 = rect_width(&bp, "center_lv0");
        let c1 = rect_width(&bp, "center_lv1");
        let c2 = rect_width(&bp, "center_lv2");
        let l0 = rect_width(&bp, "L_lv0");
        let l1 = rect_width(&bp, "L_lv1");
        let l2 = rect_width(&bp, "L_lv2");

        let f0 = c0 / (c0 + 2.0 * l0);
        let f1 = c1 / (c1 + 2.0 * l1);
        let f2 = c2 / (c2 + 2.0 * l2);

        assert!(f1 >= f0 - 1e-12);
        assert!(f2 >= f1 - 1e-12);
        assert!(f2 > f0 + 1e-6);
    }

    #[test]
    fn staged_cif_acoustic_serpentine_omits_throat_and_keeps_bypass_straight() {
        let serpentine = CenterSerpentineSpec {
            segments: 5,
            bend_radius_m: 3.5e-3,
            segment_length_m: 8.0e-3,
        };
        let bp = incremental_filtration_tri_bi_rect_staged(
            "cif-acoustic-serp",
            12e-3,
            8e-3,
            6e-3,
            2,
            2.0e-3,
            0.45,
            0.55,
            0.68,
            100e-6,
            300e-6,
            1.0e-3,
            false,
            Some(serpentine),
        );

        assert!(
            bp.channels
                .iter()
                .all(|c| c.id.as_str() != "throat_section"),
            "acoustic-only selective-routing layout must not emit a venturi throat section"
        );
        assert!(
            bp.channels
                .iter()
                .any(|c| c.id.as_str() == "treatment_section"),
            "acoustic-only selective-routing layout must emit a center treatment section"
        );

        let center_channel = bp
            .channels
            .iter()
            .find(|c| c.id.as_str() == "center_lv0")
            .expect("center_lv0 must exist");
        assert!(matches!(
            center_channel.channel_shape,
            ChannelShape::Serpentine { .. }
        ));

        let bypass_channel = bp
            .channels
            .iter()
            .find(|c| c.id.as_str() == "L_lv0")
            .expect("L_lv0 must exist");
        assert_eq!(bypass_channel.channel_shape, ChannelShape::Straight);
    }

    #[test]
    fn constriction_expansion_array_embeds_series_topology_metadata() {
        let bp =
            constriction_expansion_array_rect("ce-meta", 4, 3.0e-3, 1.5e-3, 200e-6, 80e-6, 60e-6);
        let topology = bp
            .topology_spec()
            .expect("constriction_expansion_array_rect should attach topology metadata");

        assert!(bp.is_geometry_authored());
        assert!(bp.channels.iter().all(|channel| channel.path.len() >= 2));
        assert_eq!(bp.unresolved_channel_overlap_count(), 0);
        assert_eq!(topology.treatment_channel_ids().len(), 4);
        assert!(topology.is_leukapheresis_topology());
        assert_eq!(topology.parallel_venturi_count(), 0);
    }

    #[test]
    fn spiral_channel_embeds_series_topology_metadata() {
        let bp = spiral_channel_rect("spiral-meta", 5, 6.0e-3, 200e-6, 60e-6);
        let topology = bp
            .topology_spec()
            .expect("spiral_channel_rect should attach topology metadata");

        assert!(bp.is_geometry_authored());
        assert!(bp.channels.iter().all(|channel| channel.path.len() >= 2));
        assert_eq!(bp.unresolved_channel_overlap_count(), 0);
        assert_eq!(topology.treatment_channel_ids().len(), 5);
        assert!(topology.has_serpentine());
        assert!(topology.is_leukapheresis_topology());
    }

    #[test]
    fn parallel_microchannel_array_embeds_parallel_topology_metadata() {
        let bp = parallel_microchannel_array_rect("parallel-meta", 8, 12.0e-3, 140e-6, 55e-6);
        let topology = bp
            .topology_spec()
            .expect("parallel_microchannel_array_rect should attach topology metadata");

        assert!(bp.is_geometry_authored());
        assert!(bp.channels.iter().all(|channel| channel.path.len() >= 2));
        assert_eq!(bp.unresolved_channel_overlap_count(), 0);
        assert!(topology.has_parallel_paths());
        assert_eq!(topology.parallel_channels.len(), 8);
        assert_eq!(topology.treatment_channel_ids().len(), 8);
        assert!(topology.is_leukapheresis_topology());
        assert_eq!(topology.parallel_venturi_count(), 0);
    }
}
