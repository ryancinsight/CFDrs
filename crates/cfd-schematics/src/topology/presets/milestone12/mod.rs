use super::super::model::{
    BranchRole, SerpentineSpec, SplitKind, TreatmentActuationMode, VenturiPlacementMode,
};
use super::plate_presets::{PLATE_HEIGHT_MM, PLATE_WIDTH_MM};
use aequitas::systems::si::quantities::Length;

mod build;
mod catalog;
mod support;

pub use build::{
    build_milestone12_blueprint, build_milestone12_topology_spec,
    milestone12_primitive_selective_tree_spec, promote_milestone12_option1_to_option2,
};
pub use catalog::enumerate_milestone12_topologies;
pub use support::milestone12_default_stage_layouts;

/// Branch specification within a Milestone 12 stage layout.
#[derive(Debug, Clone, PartialEq)]
pub struct Milestone12StageBranchSpec {
    /// Branch label used to derive channel identifiers.
    pub label: String,
    /// Branch role within the treatment network.
    pub role: BranchRole,
    /// Whether this branch carries the treatment path.
    pub treatment_path: bool,
    /// Authored branch width.
    pub width_m: Length<f64>,
}

/// Stage layout within a Milestone 12 topology.
#[derive(Debug, Clone, PartialEq)]
pub struct Milestone12StageLayout {
    /// Split kind of the stage.
    pub split_kind: SplitKind,
    /// Branch specifications for the stage.
    pub branches: Vec<Milestone12StageBranchSpec>,
}

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
/// [`build_milestone12_blueprint`] and
/// [`build_milestone12_topology_spec`].
pub type Milestone12TopologyRequest = Milestone12PrimitiveSelectiveSpec;

#[cfg(test)]
mod tests {
    use super::super::sequence::MILESTONE12_SWEEP_SEQUENCES;
    use super::catalog::MILESTONE12_MIRROR_VARIANTS;
    use super::*;
    use crate::topology::BlueprintTopologySpec;
    use crate::topology::{
        BlueprintTopologyFactory, BlueprintTopologyMutation, TopologyOptimizationStage,
    };

    #[test]
    fn default_stage_layouts_cover_bi_tri_quad_penta_roots() {
        let request = Milestone12PrimitiveSelectiveSpec::new(
            "m12",
            "m12",
            vec![
                SplitKind::NFurcation(2),
                SplitKind::NFurcation(3),
                SplitKind::NFurcation(4),
                SplitKind::NFurcation(5),
            ],
            Length::from_base(6.0e-3),
            Length::from_base(1.0e-3),
            Length::from_base(8.0e-3),
            Length::from_base(8.0e-3),
        );
        let layouts = milestone12_default_stage_layouts(&request);
        assert_eq!(layouts.len(), 4);
        assert_eq!(layouts[0].branches.len(), 2);
        assert_eq!(layouts[1].branches.len(), 3);
        assert_eq!(layouts[2].branches.len(), 4);
        assert_eq!(layouts[3].branches.len(), 5);
    }

    #[test]
    fn explicit_even_way_treatment_lane_sets_are_preserved() {
        let mut request = Milestone12PrimitiveSelectiveSpec::new(
            "quad",
            "quad",
            vec![SplitKind::NFurcation(4)],
            Length::from_base(6.0e-3),
            Length::from_base(1.0e-3),
            Length::from_base(8.0e-3),
            Length::from_base(8.0e-3),
        );
        request.stage_layouts = vec![Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(4),
            branches: vec![
                Milestone12StageBranchSpec {
                    label: "arm_0".to_string(),
                    role: BranchRole::Neutral,
                    treatment_path: false,
                    width_m: Length::from_base(1.0e-3),
                },
                Milestone12StageBranchSpec {
                    label: "arm_1".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: Length::from_base(2.0e-3),
                },
                Milestone12StageBranchSpec {
                    label: "arm_2".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: Length::from_base(2.0e-3),
                },
                Milestone12StageBranchSpec {
                    label: "arm_3".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: Length::from_base(1.0e-3),
                },
            ],
        }];
        request.treatment_mode = TreatmentActuationMode::VenturiCavitation;
        request.venturi_throat_count = 1;
        request.venturi_target_channel_ids = vec![
            BlueprintTopologySpec::branch_channel_id("stage_0", "arm_1"),
            BlueprintTopologySpec::branch_channel_id("stage_0", "arm_2"),
        ];

        let spec = milestone12_primitive_selective_tree_spec(&request);
        assert_eq!(spec.treatment_channel_ids().len(), 2);
        assert_eq!(spec.venturi_placements.len(), 2);
    }

    #[test]
    fn asymmetric_branch_widths_are_conserved() {
        let mut request = Milestone12PrimitiveSelectiveSpec::new(
            "tri",
            "tri",
            vec![SplitKind::NFurcation(3)],
            Length::from_base(6.0e-3),
            Length::from_base(1.0e-3),
            Length::from_base(8.0e-3),
            Length::from_base(8.0e-3),
        );
        request.stage_layouts = vec![Milestone12StageLayout {
            split_kind: SplitKind::NFurcation(3),
            branches: vec![
                Milestone12StageBranchSpec {
                    label: "left".to_string(),
                    role: BranchRole::Neutral,
                    treatment_path: false,
                    width_m: Length::from_base(1.0e-3),
                },
                Milestone12StageBranchSpec {
                    label: "center".to_string(),
                    role: BranchRole::Treatment,
                    treatment_path: true,
                    width_m: Length::from_base(2.0e-3),
                },
                Milestone12StageBranchSpec {
                    label: "right".to_string(),
                    role: BranchRole::RbcBypass,
                    treatment_path: false,
                    width_m: Length::from_base(3.0e-3),
                },
            ],
        }];

        let spec = milestone12_primitive_selective_tree_spec(&request);
        let total: f64 = spec.split_stages[0]
            .branches
            .iter()
            .map(|branch| branch.route.width_m.into_base())
            .sum();
        assert!((total - request.inlet_width_m.into_base()).abs() < 1.0e-12);
    }

    #[test]
    fn enumerate_milestone12_topologies_covers_mirrored_root_catalog() {
        let catalog = enumerate_milestone12_topologies();
        assert_eq!(
            catalog.len(),
            MILESTONE12_SWEEP_SEQUENCES.len() * MILESTONE12_MIRROR_VARIANTS.len(),
        );
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "Bi-BASE"));
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "Quad-Y"));
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "Penta-XY"));
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "QuadTriBi-BASE"));
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "PentaQuadBi-BASE"));
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "PentaQuadTri-XY"));
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "PentaTriBi-XY"));
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "QuadBi-X"));
        assert!(catalog
            .iter()
            .any(|request| request.design_name == "PentaTri-Y"));
    }

    #[test]
    fn build_milestone12_blueprint_uses_canonical_geometry_authoring() {
        let request = enumerate_milestone12_topologies()
            .into_iter()
            .find(|request| request.design_name == "Tri-XY")
            .expect("mirrored Tri scaffold should exist");

        let blueprint =
            build_milestone12_blueprint(&request).expect("Milestone 12 blueprint should build");

        assert!(blueprint.is_geometry_authored());
        assert!(blueprint
            .render_hints()
            .is_some_and(|hints| hints.mirror_x && hints.mirror_y));
        blueprint
            .validate()
            .expect("Milestone 12 blueprint should validate");
    }

    #[test]
    fn promote_option1_to_option2_rebuilds_geometry_authored_blueprint() {
        let mut request = enumerate_milestone12_topologies()
            .into_iter()
            .find(|request| request.design_name == "Quad-Y")
            .expect("mirrored Quad scaffold should exist");
        request.treatment_mode = TreatmentActuationMode::UltrasoundOnly;
        request.venturi_throat_count = 0;
        let option1 =
            build_milestone12_blueprint(&request).expect("Milestone 12 Option 1 should build");

        let promoted = promote_milestone12_option1_to_option2(
            &option1,
            1,
            VenturiPlacementMode::StraightSegment,
        )
        .expect("promotion should succeed");

        assert!(promoted.is_geometry_authored());
        assert!(promoted
            .topology_spec()
            .is_some_and(BlueprintTopologySpec::has_venturi));
        assert!(promoted
            .render_hints()
            .is_some_and(|hints| !hints.mirror_x && hints.mirror_y));
        assert_eq!(
            promoted
                .lineage()
                .and_then(|lineage| lineage.option1_source_blueprint.as_deref()),
            Some(option1.name.as_str())
        );
    }

    #[test]
    fn selective_mutation_rejects_non_geometry_authored_blueprint() {
        let request = enumerate_milestone12_topologies()
            .into_iter()
            .find(|request| request.design_name == "Tri-BASE")
            .expect("Tri scaffold should exist");
        let mut blueprint =
            build_milestone12_blueprint(&request).expect("Milestone 12 blueprint should build");
        blueprint.metadata = None;
        blueprint.geometry_authored = false;

        let error = BlueprintTopologyFactory::mutate(
            &blueprint,
            BlueprintTopologyMutation::UpdateBranchWidth {
                stage_id: "stage_0".to_string(),
                branch_label: "center".to_string(),
                new_width_m: 1.8e-3,
            },
            TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        )
        .expect_err("non-canonical selective mutation must fail");
        assert!(error.contains("create_geometry-authored provenance"));
    }
}
