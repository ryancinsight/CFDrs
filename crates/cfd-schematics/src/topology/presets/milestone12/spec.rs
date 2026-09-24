use super::layout::Milestone12StageLayout;
use crate::topology::model::{
    SerpentineSpec, SplitKind, TreatmentActuationMode, VenturiPlacementMode,
};
use crate::topology::presets::plate_presets::{PLATE_HEIGHT_MM, PLATE_WIDTH_MM};
use aequitas::systems::si::quantities::Length;

/// Canonical declarative Milestone 12 topology request.
#[derive(Debug, Clone)]
pub struct Milestone12PrimitiveSelectiveSpec {
    /// Stable topology identifier.
    pub topology_id: String,
    /// Human-readable design name used for catalog lookup.
    pub design_name: String,
    /// Mirror the geometry about the vertical axis.
    pub mirror_x: bool,
    /// Mirror the geometry about the horizontal axis.
    pub mirror_y: bool,
    /// Authored plate envelope in meters.
    pub box_dims_m: (Length<f64>, Length<f64>),
    /// Ordered split kinds per stage.
    pub split_kinds: Vec<SplitKind>,
    /// Inlet channel width.
    pub inlet_width_m: Length<f64>,
    /// Uniform channel height.
    pub channel_height_m: Length<f64>,
    /// Authored branch segment length.
    pub branch_length_m: Length<f64>,
    /// Length of the outlet tail segments.
    pub outlet_tail_length_m: Length<f64>,
    /// Explicit stage layouts, overriding derived defaults when non-empty.
    pub stage_layouts: Vec<Milestone12StageLayout>,
    /// Center-line fraction for the first trifurcation stage.
    pub first_trifurcation_center_frac: f64,
    /// Center-line fraction for later trifurcation stages.
    pub later_trifurcation_center_frac: f64,
    /// Fraction of the bifurcation consumed by the treatment path.
    pub bifurcation_treatment_frac: f64,
    /// Treatment actuation mode.
    pub treatment_mode: TreatmentActuationMode,
    /// Number of serial venturi throats.
    pub venturi_throat_count: u8,
    /// Venturi throat width.
    pub venturi_throat_width_m: Length<f64>,
    /// Venturi throat length.
    pub venturi_throat_length_m: Length<f64>,
    /// Serpentine applied to the center branch, if any.
    pub center_serpentine: Option<SerpentineSpec>,
    /// Placement mode for venturi throats.
    pub venturi_placement_mode: VenturiPlacementMode,
    /// Explicit venturi target channel identifiers.
    pub venturi_target_channel_ids: Vec<String>,
}

impl Milestone12PrimitiveSelectiveSpec {
    /// Returns the authored plate envelope in millimetres for layout geometry.
    #[must_use]
    pub fn box_dims_mm(&self) -> (f64, f64) {
        (
            self.box_dims_m.0.into_base() * 1.0e3,
            self.box_dims_m.1.into_base() * 1.0e3,
        )
    }

    /// Construct a primitive selective spec from its structural inputs,
    /// applying canonical Milestone 12 defaults for the remaining fields.
    #[must_use]
    pub fn new(
        topology_id: impl Into<String>,
        design_name: impl Into<String>,
        split_kinds: Vec<SplitKind>,
        inlet_width_m: Length<f64>,
        channel_height_m: Length<f64>,
        branch_length_m: Length<f64>,
        outlet_tail_length_m: Length<f64>,
    ) -> Self {
        Self {
            topology_id: topology_id.into(),
            design_name: design_name.into(),
            mirror_x: false,
            mirror_y: false,
            box_dims_m: (
                Length::from_base(PLATE_WIDTH_MM * 1.0e-3),
                Length::from_base(PLATE_HEIGHT_MM * 1.0e-3),
            ),
            split_kinds,
            inlet_width_m,
            channel_height_m,
            branch_length_m,
            outlet_tail_length_m,
            stage_layouts: Vec::new(),
            first_trifurcation_center_frac: 0.45,
            later_trifurcation_center_frac: 0.45,
            bifurcation_treatment_frac: 0.68,
            treatment_mode: TreatmentActuationMode::UltrasoundOnly,
            venturi_throat_count: 0,
            venturi_throat_width_m: inlet_width_m,
            venturi_throat_length_m: Length::from_base(branch_length_m.into_base() / 8.0),
            center_serpentine: None,
            venturi_placement_mode: VenturiPlacementMode::StraightSegment,
            venturi_target_channel_ids: Vec::new(),
        }
    }
}

/// Alias of the canonical Milestone 12 topology request shape.
///
/// This is the single-source request shape consumed by
/// [`build_milestone12_blueprint`](super::build_milestone12_blueprint) and
/// [`build_milestone12_topology_spec`](super::build_milestone12_topology_spec).
pub type Milestone12TopologyRequest = Milestone12PrimitiveSelectiveSpec;
