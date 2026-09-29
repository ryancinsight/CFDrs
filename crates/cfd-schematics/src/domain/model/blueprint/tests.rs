use crate::domain::therapy_metadata::TherapyZone;
use crate::geometry::metadata::ChannelVenturiSpec;
use crate::topology::presets::parallel_path_spec;
use crate::topology::{ChannelRouteSpec, ParallelChannelSpec, VenturiConfig, VenturiPlacementMode};
use crate::{BlueprintTopologyFactory, TreatmentActuationMode};
use aequitas::systems::si::quantities::{Angle, Length};

#[test]
fn add_venturi_attaches_metadata_to_existing_parallel_channel() {
    let topology = parallel_path_spec(
        "venturi-blueprint",
        2.0e-3,
        2.0e-3,
        12.0e-3,
        12.0e-3,
        vec![ParallelChannelSpec {
            channel_id: "treatment_lane".to_string(),
            route: ChannelRouteSpec {
                length_m: Length::from_base(10.0e-3),
                width_m: Length::from_base(1.6e-3),
                height_m: Length::from_base(1.0e-3),
                serpentine: None,
                therapy_zone: TherapyZone::CancerTarget,
            },
        }],
        TreatmentActuationMode::UltrasoundOnly,
    );
    let mut blueprint =
        BlueprintTopologyFactory::build(&topology).expect("parallel topology should build");

    blueprint
        .add_venturi(&VenturiConfig {
            target_channel_ids: vec!["treatment_lane".to_string()],
            serial_throat_count: 2,
            throat_geometry: crate::topology::ThroatGeometrySpec {
                throat_width_m: Length::from_base(80.0e-6),
                throat_height_m: Length::from_base(1.0e-3),
                throat_length_m: Length::from_base(300.0e-6),
                inlet_width_m: Length::from_base(0.0),
                outlet_width_m: Length::from_base(0.0),
                convergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
                divergent_half_angle: Angle::from_base(7.0_f64.to_radians()),
            },
            placement_mode: VenturiPlacementMode::StraightSegment,
        })
        .expect("venturi should attach");

    let channel = blueprint
        .channels
        .iter()
        .find(|channel| channel.id.as_str() == "treatment_lane")
        .expect("target channel must exist");
    let venturi = channel
        .venturi_geometry
        .as_ref()
        .expect("the treatment lane must carry venturi geometry");
    // A venturi constricts: the throat is narrower than its inlet.
    assert!(
        venturi.throat_width_m < venturi.inlet_width_m,
        "throat {:?} must be narrower than inlet {:?}",
        venturi.throat_width_m,
        venturi.inlet_width_m
    );
    assert_eq!(
        channel
            .metadata
            .as_ref()
            .and_then(|metadata| metadata.get::<ChannelVenturiSpec>())
            .expect("ChannelVenturiSpec must be inserted")
            .n_throats,
        2
    );
    assert_eq!(
        blueprint
            .topology
            .as_ref()
            .expect("topology metadata must be preserved")
            .venturi_placements
            .len(),
        1
    );
}
