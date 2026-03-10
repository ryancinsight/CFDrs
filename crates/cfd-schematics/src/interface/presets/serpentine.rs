use super::finalize_preset_blueprint;
use crate::topology::presets::{serpentine_bend_venturi_series_spec, serpentine_series_spec};
use crate::BlueprintTopologyFactory;

#[must_use]
pub fn serpentine_chain(
    name: impl Into<String>,
    segments: usize,
    segment_length_m: f64,
    diameter_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name = name.into();
    let topology =
        serpentine_series_spec(&name, segments, segment_length_m, diameter_m, diameter_m);
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology).expect("valid serpentine_chain topology spec"),
    )
}

#[must_use]
pub fn serpentine_rect(
    name: impl Into<String>,
    segments: usize,
    segment_length_m: f64,
    width_m: f64,
    height_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name_str = name.into();
    let topology = serpentine_series_spec(&name_str, segments, segment_length_m, width_m, height_m);
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology).expect("valid serpentine topology spec"),
    )
}

#[must_use]
pub fn serpentine_venturi_rect(
    name: impl Into<String>,
    segments: usize,
    segment_length_m: f64,
    width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    bend_radius_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name_str = name.into();
    let topology = serpentine_bend_venturi_series_spec(
        &name_str,
        segments,
        segment_length_m,
        width_m,
        throat_width_m,
        height_m,
        throat_length_m,
        bend_radius_m,
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology).expect("valid venturi topology spec"),
    )
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::model::NodeKind;

    #[test]
    fn serpentine_chain_has_correct_segment_count() {
        let bp = serpentine_chain("test", 4, 0.01, 0.001);
        let segs = bp
            .channels
            .iter()
            .filter(|c| c.id.as_str().starts_with("segment_"))
            .count();
        assert_eq!(segs, 4);
        assert!(bp.is_geometry_authored());
    }

    #[test]
    fn serpentine_rect_has_correct_endpoints() {
        let bp = serpentine_rect("test", 3, 0.01, 0.002, 0.001);
        let inlets = bp
            .nodes
            .iter()
            .filter(|n| n.kind == NodeKind::Inlet)
            .count();
        let outlets = bp
            .nodes
            .iter()
            .filter(|n| n.kind == NodeKind::Outlet)
            .count();
        assert_eq!(inlets, 1);
        assert_eq!(outlets, 1);
    }

    #[test]
    fn serpentine_rect_embeds_canonical_topology_metadata() {
        let bp = serpentine_rect("meta-serp", 4, 0.01, 0.002, 0.001);
        let topology = bp
            .topology_spec()
            .expect("serpentine_rect should attach topology metadata");

        assert_eq!(topology.design_name, "meta-serp");
        assert_eq!(topology.treatment_channel_ids().len(), 4);
        assert!(topology.has_serpentine());
        assert_eq!(topology.parallel_venturi_count(), 0);
    }

    #[test]
    fn serpentine_venturi_rect_has_venturi_at_each_bend() {
        let bp = serpentine_venturi_rect("test_sv", 5, 0.01, 0.002, 0.001, 0.001, 0.003, 0.001);

        let venturi_channels = bp.venturi_channels();
        assert_eq!(venturi_channels.len(), 4);

        let segments = bp
            .channels
            .iter()
            .filter(|c| c.id.as_str().starts_with("segment_"))
            .count();
        assert_eq!(segments, 5);
    }

    #[test]
    fn serpentine_venturi_rect_minimum_segments() {
        let bp = serpentine_venturi_rect("test_min", 1, 0.01, 0.002, 0.001, 0.001, 0.003, 0.001);
        let venturi_channels = bp.venturi_channels();
        assert_eq!(venturi_channels.len(), 1);

        let segments = bp
            .channels
            .iter()
            .filter(|c| c.id.as_str().starts_with("segment_"))
            .count();
        assert_eq!(segments, 2);
    }
}
