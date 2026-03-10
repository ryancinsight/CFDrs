use super::finalize_preset_blueprint;
use crate::topology::presets::{serial_double_venturi_series_spec, venturi_serpentine_series_spec};
use crate::BlueprintTopologyFactory;

#[must_use]
pub fn venturi_serpentine_rect(
    name: impl Into<String>,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    segments: usize,
    segment_length_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name = name.into();
    let topology = venturi_serpentine_series_spec(
        &name,
        main_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
        segments,
        segment_length_m,
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology)
            .expect("valid venturi_serpentine_rect topology spec"),
    )
}

#[must_use]
pub fn serial_double_venturi_rect(
    name: impl Into<String>,
    main_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
    inter_length_m: f64,
) -> crate::domain::model::NetworkBlueprint {
    let name = name.into();
    let topology = serial_double_venturi_series_spec(
        &name,
        main_width_m,
        throat_width_m,
        height_m,
        throat_length_m,
        inter_length_m,
    );
    finalize_preset_blueprint(
        BlueprintTopologyFactory::build(&topology)
            .expect("valid serial_double_venturi_rect topology spec"),
    )
}

#[cfg(test)]
mod tests {
    use super::{serial_double_venturi_rect, venturi_serpentine_rect};

    #[test]
    fn venturi_serpentine_rect_embeds_canonical_topology_metadata() {
        let blueprint =
            venturi_serpentine_rect("vs-meta", 2.0e-3, 0.7e-3, 1.0e-3, 1.8e-3, 4, 8.0e-3);
        let topology = blueprint
            .topology_spec()
            .expect("venturi_serpentine_rect should attach topology metadata");

        assert_eq!(topology.parallel_venturi_count(), 1);
        assert_eq!(topology.serial_venturi_stages(), 1);
        assert!(topology.has_serpentine());
        assert!(topology
            .treatment_channel_ids()
            .iter()
            .any(|channel_id| channel_id == "throat_section"));
    }

    #[test]
    fn serial_double_venturi_rect_embeds_serial_path_metadata() {
        let blueprint =
            serial_double_venturi_rect("sdv-meta", 2.0e-3, 0.7e-3, 1.0e-3, 1.8e-3, 4.0e-3);
        let topology = blueprint
            .topology_spec()
            .expect("serial_double_venturi_rect should attach topology metadata");

        assert_eq!(topology.parallel_venturi_count(), 1);
        assert_eq!(topology.serial_venturi_stages(), 2);
        assert_eq!(topology.treatment_channel_ids().len(), 2);
    }
}
