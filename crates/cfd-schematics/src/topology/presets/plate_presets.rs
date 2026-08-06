//! Shared helper functions and constants for preset construction.

use crate::domain::therapy_metadata::TherapyZone;
use aequitas::systems::si::quantities::{Angle, Length};

use super::super::model::{
    ChannelRouteSpec, ParallelChannelSpec, SeriesChannelSpec, SerpentineSpec, ThroatGeometrySpec,
};

/// Well plate dimensions \[mm] per ANSI/SLAS 1-2004.
pub(super) const PLATE_WIDTH_MM: f64 = 127.76;
pub(super) const PLATE_HEIGHT_MM: f64 = 85.47;

/// Default convergent/divergent half-angle for venturi nozzle \[degrees].
pub(super) const VENTURI_HALF_ANGLE_DEG: f64 = 7.0;

pub(super) fn default_venturi_half_angle() -> Angle<f64> {
    Angle::from_base(VENTURI_HALF_ANGLE_DEG.to_radians())
}

pub(super) fn series_channel(
    channel_id: impl Into<String>,
    length_m: f64,
    width_m: f64,
    height_m: f64,
    therapy_zone: TherapyZone,
    serpentine: Option<SerpentineSpec>,
) -> SeriesChannelSpec {
    SeriesChannelSpec {
        channel_id: channel_id.into(),
        route: ChannelRouteSpec {
            length_m: Length::from_base(length_m),
            width_m: Length::from_base(width_m),
            height_m: Length::from_base(height_m),
            serpentine,
            therapy_zone,
        },
    }
}

pub(super) fn parallel_channel(
    channel_id: impl Into<String>,
    length_m: f64,
    width_m: f64,
    height_m: f64,
    therapy_zone: TherapyZone,
    serpentine: Option<SerpentineSpec>,
) -> ParallelChannelSpec {
    ParallelChannelSpec {
        channel_id: channel_id.into(),
        route: ChannelRouteSpec {
            length_m: Length::from_base(length_m),
            width_m: Length::from_base(width_m),
            height_m: Length::from_base(height_m),
            serpentine,
            therapy_zone,
        },
    }
}

pub(super) fn throat_geometry(
    inlet_width_m: f64,
    throat_width_m: f64,
    height_m: f64,
    throat_length_m: f64,
) -> ThroatGeometrySpec {
    ThroatGeometrySpec {
        throat_width_m: Length::from_base(throat_width_m),
        throat_height_m: Length::from_base(height_m),
        throat_length_m: Length::from_base(throat_length_m),
        inlet_width_m: Length::from_base(inlet_width_m),
        outlet_width_m: Length::from_base(inlet_width_m),
        convergent_half_angle: default_venturi_half_angle(),
        divergent_half_angle: default_venturi_half_angle(),
    }
}
