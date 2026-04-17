use crate::application::objectives::evaluate_goal;
use crate::domain::{BlueprintCandidate, OptimizationGoal};
use crate::error::OptimError;
use cfd_schematics::{
    domain::model::NetworkBlueprint, promote_milestone12_option1_to_option2,
    BlueprintTopologyFactory, BlueprintTopologyMutation, SerpentineSpec, SplitKind,
    ThroatGeometrySpec, TopologyLineageEvent, TopologyLineageMetadata, TopologyOptimizationStage,
    VenturiPlacementMode,
};

fn route_serpentine(route: &cfd_schematics::ChannelRouteSpec) -> SerpentineSpec {
    route.serpentine.clone().unwrap_or(SerpentineSpec {
        segments: 4,
        bend_radius_m: (route.width_m * 2.5).max(route.height_m),
        segment_length_m: (route.length_m / 4.0).max(route.width_m),
        wave_type: cfd_schematics::SerpentineWaveType::Sine,
    })
}

fn route_throat_geometry(route: &cfd_schematics::ChannelRouteSpec) -> ThroatGeometrySpec {
    let inlet_width_m = route.width_m;
    let throat_width_m = (route.width_m * 0.22).clamp(35.0e-6, route.width_m * 0.85);
    ThroatGeometrySpec {
        throat_width_m,
        throat_height_m: route.height_m,
        throat_length_m: (route.length_m / 8.0).clamp(180.0e-6, route.length_m.max(220.0e-6)),
        inlet_width_m,
        outlet_width_m: inlet_width_m,
        convergent_half_angle_deg: 7.0,
        divergent_half_angle_deg: 7.0,
    }
}

fn stamp_lineage_metadata(
    mut blueprint: NetworkBlueprint,
    metadata: impl Into<String>,
) -> NetworkBlueprint {
    if let Some(lineage) = blueprint.lineage.as_mut() {
        if let Some(last) = lineage.mutations.last_mut() {
            last.mutation = metadata.into();
        }
    }
    blueprint
}

fn append_lineage_event(
    mut blueprint: NetworkBlueprint,
    stage: TopologyOptimizationStage,
    metadata: impl Into<String>,
    source_blueprint: impl Into<String>,
) -> NetworkBlueprint {
    let topology = blueprint.topology.clone();
    let lineage = blueprint.lineage.get_or_insert_with(|| {
        topology.as_ref().map_or_else(
            TopologyLineageMetadata::default,
            BlueprintTopologyFactory::lineage_for_spec,
        )
    });
    lineage.current_stage = stage;
    lineage.mutations.push(TopologyLineageEvent {
        stage,
        mutation: metadata.into(),
        source_blueprint: Some(source_blueprint.into()),
    });
    blueprint
}

fn merge_parent_lineages(
    mut child_blueprint: NetworkBlueprint,
    donor_blueprint: &NetworkBlueprint,
) -> NetworkBlueprint {
    match (
        child_blueprint.lineage.as_mut(),
        donor_blueprint.lineage.as_ref(),
    ) {
        (Some(child_lineage), Some(donor_lineage)) => {
            for donor_event in &donor_lineage.mutations {
                if !child_lineage
                    .mutations
                    .iter()
                    .any(|event| event == donor_event)
                {
                    child_lineage.mutations.push(donor_event.clone());
                }
            }
        }
        (None, Some(donor_lineage)) => {
            child_blueprint.lineage = Some(donor_lineage.clone());
        }
        _ => {}
    }
    child_blueprint
}

fn append_operating_point_lineage(
    blueprint: NetworkBlueprint,
    source_candidate_id: &str,
    family: &str,
    operator: &str,
) -> NetworkBlueprint {
    append_lineage_event(
        blueprint,
        TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
        format!("family={family};lane=global;operator={operator}"),
        source_candidate_id,
    )
}

fn apply_labeled_mutation(
    blueprint: &NetworkBlueprint,
    mutation: BlueprintTopologyMutation,
    stage: TopologyOptimizationStage,
    metadata: impl Into<String>,
) -> Result<NetworkBlueprint, OptimError> {
    let mutated = BlueprintTopologyFactory::mutate(blueprint, mutation, stage)
        .map_err(OptimError::InvalidParameter)?;
    Ok(stamp_lineage_metadata(mutated, metadata))
}

fn classify_serpentine_family(
    route: &cfd_schematics::ChannelRouteSpec,
    serpentine: &SerpentineSpec,
) -> &'static str {
    let neutral_segments = 4.0;
    let neutral_bend_radius_m = (route.width_m * 2.5).max(route.height_m);
    let neutral_segment_length_m = (route.length_m / 4.0).max(route.width_m);

    let segments_factor = serpentine.segments as f64 / neutral_segments;
    let bend_factor = serpentine.bend_radius_m / neutral_bend_radius_m.max(1.0e-12);
    let length_factor = serpentine.segment_length_m / neutral_segment_length_m.max(1.0e-12);

    let candidates = [
        (
            "serpentine_dense",
            positive_ratio((segments_factor - 1.75) / 0.75)
                + positive_ratio((0.85 - length_factor) / 0.15)
                + positive_ratio((0.85 - bend_factor) / 0.15),
        ),
        (
            "serpentine_compact",
            positive_ratio((1.0 - bend_factor) / 0.20)
                + positive_ratio((segments_factor - 1.25) / 0.50)
                + positive_ratio((1.0 - length_factor) / 0.20),
        ),
        (
            "serpentine_smooth",
            positive_ratio((bend_factor - 1.30) / 0.30)
                + positive_ratio((length_factor - 1.10) / 0.25),
        ),
        (
            "serpentine_long",
            positive_ratio((length_factor - 1.30) / 0.30)
                + positive_ratio((segments_factor - 1.25) / 0.75),
        ),
    ];

    candidates
        .into_iter()
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .map_or("serpentine_neutral", |(label, _)| label)
}

fn positive_ratio(value: f64) -> f64 {
    value.max(0.0)
}

fn serpentine_variants(route: &cfd_schematics::ChannelRouteSpec) -> Vec<(String, SerpentineSpec)> {
    let base = route_serpentine(route);
    use cfd_schematics::SerpentineWaveType;

    let candidates = [
        // Sine variants: compact, spread, dense, long
        (
            "sc",
            SerpentineSpec {
                segments: (base.segments + 1).max(3),
                bend_radius_m: (base.bend_radius_m * 0.8).max(route.height_m),
                segment_length_m: (base.segment_length_m * 0.9).max(route.width_m),
                wave_type: SerpentineWaveType::Sine,
            },
        ),
        (
            "ss",
            SerpentineSpec {
                segments: base.segments.saturating_sub(1).max(2),
                bend_radius_m: (base.bend_radius_m * 1.3).max(route.height_m),
                segment_length_m: (base.segment_length_m * 1.15).max(route.width_m),
                wave_type: SerpentineWaveType::Sine,
            },
        ),
        (
            "sd",
            SerpentineSpec {
                segments: (base.segments + 3).max(4),
                bend_radius_m: (base.bend_radius_m * 0.65).max(route.height_m),
                segment_length_m: (base.segment_length_m * 0.75).max(route.width_m),
                wave_type: SerpentineWaveType::Sine,
            },
        ),
        (
            "sl",
            SerpentineSpec {
                segments: (base.segments + 1).max(3),
                bend_radius_m: (base.bend_radius_m * 1.05).max(route.height_m),
                segment_length_m: (base.segment_length_m * 1.35).max(route.width_m),
                wave_type: SerpentineWaveType::Sine,
            },
        ),
        // Square wave variant: sharp U-turns with higher K-factor losses
        // but most compact footprint.  The BendType::Sharp in the 1D
        // solver gives K = 2.2 + 250/Re per bend vs ~0.3-0.9 for smooth.
        (
            "sq",
            SerpentineSpec {
                segments: base.segments.max(3),
                bend_radius_m: (base.bend_radius_m * 0.6).max(route.height_m),
                segment_length_m: (base.segment_length_m * 0.85).max(route.width_m),
                wave_type: SerpentineWaveType::Square,
            },
        ),
        // Triangular wave variant: linear ramps with sharp apices.
        // Zero curvature along ramps, localized Dean peaks at apices.
        // Lower path length than sine (2/pi ratio) for same amplitude.
        (
            "tr",
            SerpentineSpec {
                segments: (base.segments + 1).max(3),
                bend_radius_m: (base.bend_radius_m * 0.9).max(route.height_m),
                segment_length_m: (base.segment_length_m * 0.95).max(route.width_m),
                wave_type: SerpentineWaveType::Triangular,
            },
        ),
    ];
    let mut dedup = std::collections::HashSet::with_capacity(candidates.len());
    let mut variants = Vec::with_capacity(candidates.len());
    for (label, spec) in candidates {
        let wave_tag = match spec.wave_type {
            SerpentineWaveType::Sine => 0_i64,
            SerpentineWaveType::Square => 1,
            SerpentineWaveType::Triangular => 2,
        };
        let key = (
            spec.segments,
            (spec.bend_radius_m * 1.0e9).round() as i64,
            (spec.segment_length_m * 1.0e9).round() as i64,
            wave_tag,
        );
        if dedup.insert(key) {
            variants.push((label.to_string(), spec));
        }
    }
    variants
}

fn strongest_venturi_template(
    candidate: &BlueprintCandidate,
) -> Result<Option<(u8, ThroatGeometrySpec, VenturiPlacementMode)>, OptimError> {
    Ok(candidate
        .topology_spec()?
        .venturi_placements
        .iter()
        .max_by_key(|placement| placement.serial_throat_count)
        .map(|placement| {
            (
                placement.serial_throat_count.max(1),
                placement.throat_geometry.clone(),
                placement.placement_mode,
            )
        }))
}

fn compatible_treatment_pairs(
    left: &BlueprintCandidate,
    right: &BlueprintCandidate,
) -> Result<Vec<(String, String)>, OptimError> {
    let left_topology = left.topology_spec()?;
    let right_topology = right.topology_spec()?;
    if left_topology.split_stages.len() != right_topology.split_stages.len() {
        return Ok(Vec::new());
    }

    let left_channels = left.treatment_channel_ids();
    let right_channels = right.treatment_channel_ids();
    if left_channels.len() != right_channels.len() || left_channels.is_empty() {
        return Ok(Vec::new());
    }

    Ok(left_channels.into_iter().zip(right_channels).collect())
}

fn crossover_child(
    base: &BlueprintCandidate,
    donor: &BlueprintCandidate,
    label: &str,
    operating_point: crate::domain::OperatingPoint,
) -> Result<Option<BlueprintCandidate>, OptimError> {
    let treatment_pairs = compatible_treatment_pairs(base, donor)?;
    if treatment_pairs.is_empty() {
        return Ok(None);
    }

    let donor_topology = donor.topology_spec()?;
    let donor_venturi = strongest_venturi_template(donor)?;
    let mut child_blueprint = merge_parent_lineages(base.blueprint.clone(), &donor.blueprint);

    for (base_channel_id, donor_channel_id) in treatment_pairs {
        if let Some(donor_route) = donor_topology.channel_route(&donor_channel_id) {
            if donor_route.serpentine.is_some() {
                let family = donor_route
                    .serpentine
                    .as_ref()
                    .map_or("serpentine_neutral", |serpentine| {
                        classify_serpentine_family(donor_route, serpentine)
                    });
                child_blueprint = apply_labeled_mutation(
                    &child_blueprint,
                    BlueprintTopologyMutation::SetTreatmentChannelSerpentine {
                        target_channel_id: base_channel_id.clone(),
                        serpentine: donor_route.serpentine.clone(),
                    },
                    TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                    format!("family={family};lane={base_channel_id};operator=crossover_serpentine_transfer"),
                )
                ?;
            }
        }
        if let Some((serial_throat_count, throat_geometry, placement_mode)) = &donor_venturi {
            child_blueprint = apply_labeled_mutation(
                &child_blueprint,
                BlueprintTopologyMutation::SetTreatmentChannelVenturi {
                    target_channel_id: base_channel_id,
                    serial_throat_count: *serial_throat_count,
                    throat_geometry: throat_geometry.clone(),
                    placement_mode: *placement_mode,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                format!(
                    "family=venturi_transfer;lane={donor_channel_id};operator=crossover_venturi_transfer"
                ),
            )
            ?;
        }
    }

    child_blueprint = append_operating_point_lineage(
        child_blueprint,
        &base.id,
        &format!("operating_point_{label}"),
        "operating_point_crossover",
    );

    Ok(Some(BlueprintCandidate::new(
        format!("{}-x-{}", base.id, label),
        child_blueprint,
        operating_point,
    )))
}

pub fn generate_ga_crossover_children(
    left: &BlueprintCandidate,
    right: &BlueprintCandidate,
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let mut children = Vec::with_capacity(3);
    if let Some(child) = crossover_child(
        left,
        right,
        "qrp",
        crate::domain::OperatingPoint {
            flow_rate_m3_s: left.operating_point.flow_rate_m3_s,
            inlet_gauge_pa: right.operating_point.inlet_gauge_pa,
            feed_hematocrit: 0.5
                * (left.operating_point.feed_hematocrit + right.operating_point.feed_hematocrit),
            patient_context: left.operating_point.patient_context.clone(),
        },
    )? {
        children.push(child);
    }
    if let Some(child) = crossover_child(
        right,
        left,
        "qpr",
        crate::domain::OperatingPoint {
            flow_rate_m3_s: right.operating_point.flow_rate_m3_s,
            inlet_gauge_pa: left.operating_point.inlet_gauge_pa,
            feed_hematocrit: 0.5
                * (left.operating_point.feed_hematocrit + right.operating_point.feed_hematocrit),
            patient_context: right.operating_point.patient_context.clone(),
        },
    )? {
        children.push(child);
    }
    if let Some(child) = crossover_child(
        left,
        right,
        "blend",
        crate::domain::OperatingPoint {
            flow_rate_m3_s: 0.5
                * (left.operating_point.flow_rate_m3_s + right.operating_point.flow_rate_m3_s),
            inlet_gauge_pa: 0.5
                * (left.operating_point.inlet_gauge_pa + right.operating_point.inlet_gauge_pa),
            feed_hematocrit: 0.5
                * (left.operating_point.feed_hematocrit + right.operating_point.feed_hematocrit),
            patient_context: left.operating_point.patient_context.clone(),
        },
    )? {
        children.push(child);
    }

    Ok(children)
}

pub fn seed_option2_candidates(candidates: &[BlueprintCandidate]) -> Vec<&BlueprintCandidate> {
    candidates
        .iter()
        .filter(|candidate| {
            candidate.inferred_goal()
                == OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity
                || candidate
                    .blueprint
                    .topology_spec()
                    .is_some_and(|spec| !spec.venturi_placements.is_empty())
        })
        .collect()
}

pub fn generate_ga_mutations(
    seed: &BlueprintCandidate,
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let topology = seed.topology_spec()?;
    if topology.venturi_placements.is_empty() {
        return Err(OptimError::InvalidParameter(format!(
            "GA refinement requires canonical venturi-capable seeds, but candidate '{}' has no venturi placements",
            seed.id
        )));
    }
    let treatment_ids = seed.treatment_channel_ids();
    // Each treatment channel produces ~20 mutations (4 venturi + 1 serpentine
    // + 12 split-merge variants) plus 8 operating-point perturbations.
    let mut mutated = Vec::with_capacity(treatment_ids.len() * 20 + 8);

    for target_channel_id in treatment_ids {
        let route = topology
            .channel_route(&target_channel_id)
            .ok_or_else(|| {
                OptimError::InvalidParameter(format!(
                    "treatment channel '{}' not found in seed '{}'",
                    target_channel_id, seed.id
                ))
            })?
            .clone();
        let venturi_geometry = route_throat_geometry(&route);
        let serpentine_variants = serpentine_variants(&route);
        for (variant_label, variant_serpentine) in &serpentine_variants {
            let family = classify_serpentine_family(&route, variant_serpentine);
            let serpentine_mutation = apply_labeled_mutation(
                &seed.blueprint,
                BlueprintTopologyMutation::SetTreatmentChannelSerpentine {
                    target_channel_id: target_channel_id.clone(),
                    serpentine: Some(variant_serpentine.clone()),
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                format!("family={family};lane={target_channel_id};operator=serpentine_variant"),
            )?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-s{}-{}", seed.id, variant_label, target_channel_id),
                serpentine_mutation.clone(),
                seed.operating_point.clone(),
            ));

            let serpentine_dean_venturi = apply_labeled_mutation(
                &serpentine_mutation,
                BlueprintTopologyMutation::SetTreatmentChannelVenturi {
                    target_channel_id: target_channel_id.clone(),
                    serial_throat_count: 4,
                    throat_geometry: venturi_geometry.clone(),
                    placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                format!(
                    "family=venturi_dean;lane={target_channel_id};operator=serpentine_venturi_pair"
                ),
            )?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-svd{}-{}", seed.id, variant_label, target_channel_id),
                serpentine_dean_venturi,
                seed.operating_point.clone(),
            ));
        }

        // Venturi mutations with varied throat counts and widths.
        for &(vt_count, width_factor, label) in &[
            (2_u8, 1.0, "v2"),
            (3, 1.0, "v3"),
            (2, 0.7, "v2n"), // narrower throat → lower σ
            (2, 1.3, "v2w"), // wider throat → higher σ, less hemolysis
        ] {
            let mut varied_geom = venturi_geometry.clone();
            varied_geom.throat_width_m *= width_factor;
            let venturi_mutation = apply_labeled_mutation(
                &seed.blueprint,
                BlueprintTopologyMutation::SetTreatmentChannelVenturi {
                    target_channel_id: target_channel_id.clone(),
                    serial_throat_count: vt_count,
                    throat_geometry: varied_geom,
                    placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                format!("family=venturi_{label};lane={target_channel_id};operator=venturi_variant"),
            )?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-{}-{}", seed.id, label, target_channel_id),
                venturi_mutation,
                seed.operating_point.clone(),
            ));
        }

        for split_kind in [
            SplitKind::NFurcation(2),
            SplitKind::NFurcation(3),
            SplitKind::NFurcation(4),
            SplitKind::NFurcation(5),
        ] {
            let split_merge = apply_labeled_mutation(
                &seed.blueprint,
                BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                    target_channel_id: target_channel_id.clone(),
                    split_kind,
                    treatment_serpentine: None,
                    venturi_serial_throat_count: None,
                    venturi_throat_geometry: None,
                    venturi_placement_mode: VenturiPlacementMode::StraightSegment,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                format!("family=split_merge_{split_kind:?};lane={target_channel_id};operator=split_merge"),
            )
            ?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-sm-{}-{:?}", seed.id, target_channel_id, split_kind),
                split_merge,
                seed.operating_point.clone(),
            ));

            let split_merge_serpentine = apply_labeled_mutation(
                &seed.blueprint,
                BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                    target_channel_id: target_channel_id.clone(),
                    split_kind,
                    treatment_serpentine: serpentine_variants.first().map(|(_, spec)| spec.clone()),
                    venturi_serial_throat_count: None,
                    venturi_throat_geometry: None,
                    venturi_placement_mode: VenturiPlacementMode::StraightSegment,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                format!("family=split_merge_serpentine_{split_kind:?};lane={target_channel_id};operator=split_merge_serpentine"),
            )
            ?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-sms-{}-{:?}", seed.id, target_channel_id, split_kind),
                split_merge_serpentine,
                seed.operating_point.clone(),
            ));

            let split_merge_venturi = apply_labeled_mutation(
                &seed.blueprint,
                BlueprintTopologyMutation::InsertTreatmentSplitMerge {
                    target_channel_id: target_channel_id.clone(),
                    split_kind,
                    treatment_serpentine: serpentine_variants.get(1).map(|(_, spec)| spec.clone()).or_else(|| serpentine_variants.first().map(|(_, spec)| spec.clone())),
                    venturi_serial_throat_count: Some(2),
                    venturi_throat_geometry: Some(venturi_geometry.clone()),
                    venturi_placement_mode: VenturiPlacementMode::CurvaturePeakDeanNumber,
                },
                TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                format!("family=split_merge_venturi_{split_kind:?};lane={target_channel_id};operator=split_merge_venturi"),
            )
            ?;
            mutated.push(BlueprintCandidate::new(
                format!("{}-ga-smv-{}-{:?}", seed.id, target_channel_id, split_kind),
                split_merge_venturi,
                seed.operating_point.clone(),
            ));
        }
    }

    // Operating-point perturbation mutations: vary flow rate and pressure
    // to explore the physical design space around the seed topology.
    // These produce genuinely different scores because they change the
    // 1D CFD solve inputs (Re, ΔP, σ), unlike stacked serpentine
    // mutations which are nearly idempotent on the same channel.
    let base_q = seed.operating_point.flow_rate_m3_s;
    let base_p = seed.operating_point.inlet_gauge_pa;
    let base_ht = seed.operating_point.feed_hematocrit;
    for &q_factor in &[0.7, 0.85, 1.15, 1.3] {
        let blueprint = append_operating_point_lineage(
            seed.blueprint.clone(),
            &seed.id,
            &format!("operating_point_qf{:.0}", q_factor * 100.0),
            "operating_point_scan",
        );
        mutated.push(BlueprintCandidate::new(
            format!("{}-ga-qf{:.0}", seed.id, q_factor * 100.0),
            blueprint,
            crate::domain::OperatingPoint {
                flow_rate_m3_s: base_q * q_factor,
                inlet_gauge_pa: base_p,
                feed_hematocrit: base_ht,
                patient_context: seed.operating_point.patient_context.clone(),
            },
        ));
    }
    for &p_factor in &[0.7, 0.85, 1.15, 1.3] {
        let blueprint = append_operating_point_lineage(
            seed.blueprint.clone(),
            &seed.id,
            &format!("operating_point_pf{:.0}", p_factor * 100.0),
            "operating_point_scan",
        );
        mutated.push(BlueprintCandidate::new(
            format!("{}-ga-pf{:.0}", seed.id, p_factor * 100.0),
            blueprint,
            crate::domain::OperatingPoint {
                flow_rate_m3_s: base_q,
                inlet_gauge_pa: base_p * p_factor,
                feed_hematocrit: base_ht,
                patient_context: seed.operating_point.patient_context.clone(),
            },
        ));
    }

    // Branch width scaling mutations: adjust the Zweifach-Fung routing
    // by changing the relative widths at each split junction.  This is
    // the ONLY way to change cell routing metrics (cancer_center_fraction,
    // rbc_venturi_exposure) since these depend on the flow-rate split
    // which is determined by the rectangular branch-conductance weighting.
    //
    // Wider treatment branches → more flow to center → more CTCs and RBCs
    // Narrower treatment branches → less flow → fewer CTCs but also fewer RBCs
    // The optimizer must find the width that maximizes CTC enrichment while
    // minimizing RBC exposure.
    for stage in &topology.split_stages {
        for branch in &stage.branches {
            if !branch.treatment_path {
                continue; // Only scale treatment branches for now
            }
            let base_w = branch.route.width_m;
            for &(factor, label) in &[
                (1.10, "w110"), // 10% wider → more flow → more cells (both CTC and RBC)
                (1.20, "w120"), // 20% wider
                (0.90, "w90"),  // 10% narrower → less flow → fewer cells
                (0.80, "w80"),  // 20% narrower
            ] {
                let new_w = (base_w * factor).clamp(50e-6, 10e-3); // clamp to physical range
                let width_mutation = apply_labeled_mutation(
                    &seed.blueprint,
                    BlueprintTopologyMutation::UpdateBranchWidth {
                        stage_id: stage.stage_id.clone(),
                        branch_label: branch.label.clone(),
                        new_width_m: new_w,
                    },
                    TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
                    format!(
                        "family=width_{label};lane={}/{};operator=branch_width_scale",
                        stage.stage_id, branch.label
                    ),
                )?;
                mutated.push(BlueprintCandidate::new(
                    format!(
                        "{}-ga-{}-{}-{}",
                        seed.id, label, stage.stage_id, branch.label
                    ),
                    width_mutation,
                    seed.operating_point.clone(),
                ));
            }
        }
    }

    Ok(mutated)
}

#[cfg(test)]
mod tests {
    use crate::domain::fixtures::{canonical_option2_candidate, operating_point};

    use super::{generate_ga_crossover_children, generate_ga_mutations};

    #[test]
    fn ga_crossover_children_remain_valid() {
        let seed = canonical_option2_candidate("seed", operating_point(2.0e-6, 30_000.0, 0.18));
        let donor_topology = generate_ga_mutations(&seed)
            .expect("mutations")
            .into_iter()
            .find(|candidate| candidate.id.contains("-ga-s"))
            .expect("serpentine mutation");
        let donor = crate::domain::BlueprintCandidate::new(
            "donor".to_string(),
            donor_topology.blueprint.clone(),
            crate::domain::OperatingPoint {
                flow_rate_m3_s: seed.operating_point.flow_rate_m3_s * 1.1,
                inlet_gauge_pa: seed.operating_point.inlet_gauge_pa * 0.9,
                feed_hematocrit: seed.operating_point.feed_hematocrit,
                patient_context: seed.operating_point.patient_context.clone(),
            },
        );

        let children = generate_ga_crossover_children(&seed, &donor).expect("crossovers");
        assert!(!children.is_empty());
        for child in children {
            let topology = child.topology_spec().expect("topology metadata");
            cfd_schematics::BlueprintTopologyFactory::validate_spec(topology)
                .expect("crossover topology must remain valid");
        }
    }
}

pub fn promote_option1_candidate_to_ga_seed(
    seed: &BlueprintCandidate,
) -> Result<BlueprintCandidate, OptimError> {
    let promoted = promote_milestone12_option1_to_option2(
        &seed.blueprint,
        1,
        VenturiPlacementMode::CurvaturePeakDeanNumber,
    )
    .map_err(OptimError::InvalidParameter)?;

    let promoted_seed = BlueprintCandidate::new(
        format!("{}-ga-promoted", seed.id),
        promoted,
        seed.operating_point.clone(),
    );
    let topology = promoted_seed.topology_spec()?;
    let target_channel_id = promoted_seed
        .treatment_channel_ids()
        .into_iter()
        .next()
        .ok_or_else(|| {
            OptimError::InvalidParameter(format!(
                "promoted GA seed '{}' has no treatment channels",
                promoted_seed.id
            ))
        })?;
    let route = topology.channel_route(&target_channel_id).ok_or_else(|| {
        OptimError::InvalidParameter(format!(
            "promoted GA seed '{}' could not resolve treatment channel '{}'",
            promoted_seed.id, target_channel_id
        ))
    })?;
    let serpentine_seed = BlueprintTopologyFactory::mutate(
        promoted_seed.blueprint(),
        BlueprintTopologyMutation::SetTreatmentChannelSerpentine {
            target_channel_id: target_channel_id.clone(),
            serpentine: Some(route_serpentine(route)),
        },
        TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
    )
    .map_err(OptimError::InvalidParameter)?;

    Ok(BlueprintCandidate::new(
        format!("{}-ga-promoted", seed.id),
        serpentine_seed,
        seed.operating_point.clone(),
    ))
}

pub fn promote_option2_candidate_to_ga_seed(
    seed: &BlueprintCandidate,
) -> Result<BlueprintCandidate, OptimError> {
    let topology = seed.topology_spec()?;
    if topology.venturi_placements.is_empty() {
        return Err(OptimError::InvalidParameter(format!(
            "Option 2 GA seed '{}' has no venturi placements",
            seed.id
        )));
    }

    let mut promoted = topology.clone();
    for placement in &mut promoted.venturi_placements {
        placement.placement_mode = VenturiPlacementMode::CurvaturePeakDeanNumber;
    }
    let promoted_blueprint =
        BlueprintTopologyFactory::build(&promoted).map_err(OptimError::InvalidParameter)?;
    let promoted_seed = BlueprintCandidate::new(
        format!("{}-ga-dean-seed", seed.id),
        promoted_blueprint,
        seed.operating_point.clone(),
    );
    let target_channel_id = promoted_seed
        .treatment_channel_ids()
        .into_iter()
        .next()
        .ok_or_else(|| {
            OptimError::InvalidParameter(format!(
                "promoted Option 2 GA seed '{}' has no treatment channels",
                promoted_seed.id
            ))
        })?;
    let route = promoted_seed
        .topology_spec()?
        .channel_route(&target_channel_id)
        .ok_or_else(|| {
            OptimError::InvalidParameter(format!(
                "promoted Option 2 GA seed '{}' could not resolve treatment channel '{}'",
                promoted_seed.id, target_channel_id
            ))
        })?;
    let serpentine_blueprint = BlueprintTopologyFactory::mutate(
        promoted_seed.blueprint(),
        BlueprintTopologyMutation::SetTreatmentChannelSerpentine {
            target_channel_id: target_channel_id.clone(),
            serpentine: Some(route_serpentine(route)),
        },
        TopologyOptimizationStage::InPlaceDeanSerpentineRefinement,
    )
    .map_err(OptimError::InvalidParameter)?;

    Ok(BlueprintCandidate::new(
        format!("{}-ga-dean-seed", seed.id),
        serpentine_blueprint,
        seed.operating_point.clone(),
    ))
}

pub fn build_milestone12_ga_seed_pair(
    candidates: &[BlueprintCandidate],
) -> Result<Vec<BlueprintCandidate>, OptimError> {
    let option1 = candidates
        .iter()
        .filter(|candidate| {
            candidate
                .topology_spec()
                .is_ok_and(|topology| !topology.has_venturi() && !topology.split_stages.is_empty())
        })
        .map(|candidate| {
            evaluate_goal(
                candidate,
                OptimizationGoal::AsymmetricSplitResidenceSeparation,
            )
            .map(|evaluation| (candidate, evaluation.score.unwrap_or(0.0)))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .map(|(candidate, _)| candidate.clone())
        .ok_or_else(|| {
            OptimError::InvalidParameter(
                "Milestone 12 GA seed selection found no Option 1 baseline candidate".to_string(),
            )
        })?;
    let option2 = candidates
        .iter()
        .filter(|candidate| {
            candidate
                .topology_spec()
                .is_ok_and(|topology| !topology.venturi_placements.is_empty())
        })
        .map(|candidate| {
            evaluate_goal(
                candidate,
                OptimizationGoal::AsymmetricSplitVenturiCavitationSelectivity,
            )
            .map(|evaluation| (candidate, evaluation.score.unwrap_or(0.0)))
        })
        .collect::<Result<Vec<_>, _>>()?
        .into_iter()
        .max_by(|left, right| left.1.total_cmp(&right.1))
        .map(|(candidate, _)| candidate.clone())
        .ok_or_else(|| {
            OptimError::InvalidParameter(
                "Milestone 12 GA seed selection found no Option 2 baseline candidate".to_string(),
            )
        })?;

    Ok(vec![
        promote_option1_candidate_to_ga_seed(&option1)?,
        promote_option2_candidate_to_ga_seed(&option2)?,
    ])
}
