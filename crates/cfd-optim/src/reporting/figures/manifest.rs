//! Dynamic figure generation and manifest assembly for Milestone 12 narrative.

use std::path::Path;

use cfd_schematics::domain::model::ChannelShape;
use cfd_schematics::geometry::metadata::ChannelVisualRole;
use serde::{Deserialize, Serialize};

use super::process::{write_creation_optimization_process_figure, write_placeholder};
use super::svg::{
    write_cavitation_distribution_figure, write_cross_mode_figure,
    write_dean_venturi_placement_figure, write_ga_convergence_figure,
    write_head_to_head_figure, write_pareto_figure, write_pediatric_ecv_figure,
    DeanVenturiPoint,
};
use super::treatment_lane::write_treatment_lane_zoom_figure;
use crate::constraints::{
    BLOOD_DENSITY_KG_M3, BLOOD_VAPOR_PRESSURE_PA, BLOOD_VISCOSITY_PA_S, DIFFUSER_DISCHARGE_COEFF,
    PEDIATRIC_BLOOD_VOLUME_ML_PER_KG, PEDIATRIC_REFERENCE_WEIGHT_KG, VENTURI_CC,
};
use crate::reporting::{rank_ga_hydrosdt_report_designs, Milestone12ReportDesign, ParetoPoint};

/// Figure metadata used for dynamic table-of-contents and section rendering.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct NarrativeFigureSpec {
    pub number: usize,
    pub title: String,
    pub path: String,
    pub caption: String,
    pub alt: String,
    /// True when this figure was generated as a placeholder because the source
    /// schematic or data was not available in the current run.
    #[serde(default)]
    pub is_placeholder: bool,
}

/// Inputs needed to generate dynamic report figures.
pub struct FigureGenerationInput<'a> {
    pub option1_ranked: &'a [Milestone12ReportDesign],
    pub option2_ranked: &'a [Milestone12ReportDesign],
    pub ga_top: &'a [Milestone12ReportDesign],
    pub option2_pool_all: &'a [ParetoPoint],
    pub ga_pool_all: &'a [ParetoPoint],
    pub ga_best_per_gen: &'a [f64],
    pub fast_mode: bool,
}

/// Generate report figure assets and return ordered figure specs.
///
/// # Errors
/// Returns an error if any figure cannot be written.
pub fn generate_m12_report_figures(
    figures_dir: &Path,
    figure_path_prefix: &str,
    input: &FigureGenerationInput<'_>,
) -> Result<Vec<NarrativeFigureSpec>, Box<dyn std::error::Error>> {
    let ga_ranked_for_figures: Vec<Milestone12ReportDesign> = input
        .option2_ranked
        .first()
        .map(|baseline_option2| rank_ga_hydrosdt_report_designs(input.ga_top, baseline_option2))
        .filter(|filtered| !filtered.is_empty())
        .unwrap_or_else(|| input.ga_top.to_vec());

    std::fs::create_dir_all(figures_dir)?;

    write_creation_optimization_process_figure(
        &figures_dir.join("milestone12_creation_optimization_process.svg"),
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("treatment_zone_plate.svg"),
        "Treatment Zone Layout - Bifurcation Concept",
    )?;
    ensure_existing_or_placeholder(
        &figures_dir.join("treatment_zone_plate_trifurcation.svg"),
        "Treatment Zone Layout - Trifurcation Concept",
    )?;
    let option1_is_placeholder = if input.option1_ranked.is_empty() {
        write_placeholder(
            &figures_dir.join("selected_option1_schematic.svg"),
            "Option 1 Unavailable Under Current Physics",
            "No selective acoustic design satisfied strict eligibility under the current physics regime.",
        )?;
        true
    } else {
        ensure_existing_or_placeholder(
            &figures_dir.join("selected_option1_schematic.svg"),
            "Selected Option 1 Design - Selective Acoustic",
        )?
    };
    let option2_is_placeholder = ensure_existing_or_placeholder(
        &figures_dir.join("selected_option2_combined_schematic.svg"),
        "Selected Option 2 Design - Combined Selective Venturi",
    )?;
    let ga_is_placeholder = ensure_existing_or_placeholder(
        &figures_dir.join("top_hydrosdt_schematic.svg"),
        "Best HydroSDT GA-Optimized Design",
    )?;

    write_cross_mode_figure(
        &figures_dir.join("m12_cross_mode_scoring.svg"),
        input.option1_ranked,
        input.option2_ranked,
        &ga_ranked_for_figures,
    )?;
    write_head_to_head_figure(
        &figures_dir.join("m12_head_to_head_scores.svg"),
        input.option2_ranked,
    )?;
    write_cavitation_distribution_figure(
        &figures_dir.join("m12_cavitation_distribution.svg"),
        input.option2_ranked,
        &ga_ranked_for_figures,
    )?;
    write_pareto_figure(
        &figures_dir.join("m12_pareto_oncology.svg"),
        input.option2_ranked,
        &ga_ranked_for_figures,
        input.option2_pool_all,
        input.ga_pool_all,
    )?;
    write_pediatric_ecv_figure(
        &figures_dir.join("m12_pediatric_ecv_margin.svg"),
        input.option1_ranked.first(),
        first_ranked(input.option2_ranked),
        first_ranked(&ga_ranked_for_figures),
    )?;
    if !input.ga_best_per_gen.is_empty() {
        write_ga_convergence_figure(
            &figures_dir.join("m12_ga_convergence.svg"),
            input.ga_best_per_gen,
            input.fast_mode,
        )?;
    }

    // Dean-venturi placement figure: recompute per-placement metrics from the
    // GA top design to show how bend position influences cavitation.
    let option2_serpentine_focus = pick_serpentine_venturi_focus(input.option2_ranked)
        .unwrap_or(&input.option2_ranked[0]);
    let ga_serpentine_focus = pick_serpentine_venturi_focus(&ga_ranked_for_figures)
        .unwrap_or(&ga_ranked_for_figures[0]);
    let dean_venturi_points = extract_dean_venturi_points(Some(ga_serpentine_focus));
    write_dean_venturi_placement_figure(
        &figures_dir.join("m12_dean_venturi_placement.svg"),
        "Dean Number vs Cavitation at Serpentine Bend Apices",
        &dean_venturi_points,
    )?;

    let option1 = input.option1_ranked.first();
    let option2 = &input.option2_ranked[0];
    let ga_best = &ga_ranked_for_figures[0];

    write_treatment_lane_zoom_figure(
        &figures_dir.join("m12_treatment_lane_zoom.svg"),
        option2_serpentine_focus,
        ga_serpentine_focus,
    )?;

    let make_spec = |number, title: &str, file_name, prefix, caption: &str, alt, is_ph| {
        if is_ph {
            spec_placeholder(number, title, file_name, prefix, caption, alt)
        } else {
            spec(number, title, file_name, prefix, caption, alt)
        }
    };

    let mut specs = vec![
        option1.map_or_else(
            || {
                spec_placeholder(
                    4,
                    "Option 1 Unavailable Under Current Physics",
                    "selected_option1_schematic.svg",
                    figure_path_prefix,
                    "No selective acoustic design satisfied strict eligibility under the current physics regime, so Figure 4 is an explicit placeholder documenting the empty Option 1 shortlist.",
                    "Option 1 unavailable placeholder",
                )
            },
            |option1| {
                make_spec(
                    4,
                    &format!(
                        "Selected Option 1 Design - {} Selective Acoustic",
                        stage_sequence_label(option1)
                    ),
                    "selected_option1_schematic.svg",
                    figure_path_prefix,
                    &selected_schematic_caption(
                        option1,
                        "Option 1",
                        "selective acoustic center-treatment network",
                    ),
                    "Option 1 schematic",
                    option1_is_placeholder,
                )
            },
        ),
        make_spec(
            5,
            &format!(
                "Selected Option 2 Design - {} Combined Selective Venturi",
                stage_sequence_label(option2)
            ),
            "selected_option2_combined_schematic.svg",
            figure_path_prefix,
            &selected_schematic_caption(
                option2,
                "Option 2",
                "selective venturi treatment ranked by the combined score",
            ),
            "Option 2 schematic",
            option2_is_placeholder,
        ),
        make_spec(
            6,
            &format!(
                "Best HydroSDT GA-Optimized Design - {}",
                stage_sequence_label(ga_best)
            ),
            "top_hydrosdt_schematic.svg",
            figure_path_prefix,
            &format!(
                "Best HydroSDT GA-optimized design after architecture-preserving in-place refinement. Topology: {}. Visible split layers: {} ({}). Active venturi throats: {} with {} serial stage(s) per path. The treatment-path geometry shows where Dean-generating curvature and venturi placement were co-localised for the final GA-ranked lineage.",
                ga_best.topology_display_name(),
                visible_split_layers(ga_best),
                stage_sequence_label(ga_best),
                ga_best.metrics.active_venturi_throat_count,
                ga_best.metrics.serial_venturi_stages_per_path,
            ),
            "GA best schematic",
            ga_is_placeholder,
        ),
        spec(
            7,
            "Cross-Mode Scoring Comparison",
            "m12_cross_mode_scoring.svg",
            figure_path_prefix,
            "Cross-mode scoring comparison for the selected Option 1, Option 2, and GA tracks. Bars should be interpreted within-track because each optimization goal uses a distinct composite objective; the figure is intended to show relative ranking structure, not absolute score equivalence across tracks.",
            "Cross mode score bars",
        ),
        spec(
            8,
            "Head-to-Head Design Comparison",
            "m12_head_to_head_scores.svg",
            figure_path_prefix,
            "Head-to-head comparison of the top-ranked Option 2 candidates under the deterministic tie-break sequence. Differences reflect changes in asymmetric split-width partitioning, venturi stage count, and throat geometry rather than unrelated rendering variation.",
            "Option 2 top 5 score bars",
        ),
        spec(
            9,
            "Selected-Design Cavitation Number (σ)",
            "m12_cavitation_distribution.svg",
            figure_path_prefix,
            "Selected-design cavitation-number values for the top-ranked Option 2 and GA venturi designs, with reference lines at σ = 0 and σ = 1. The plot preserves linear detail in the cavitating regime (σ < 1) and log-compresses the non-cavitating branch (σ > 1) so extreme positive GA values do not visually collapse the clinically relevant range near inception.",
            "Selected design sigma scatter",
        ),
        spec(
            10,
            "Selected-Design Oncology Trade-Off Frontier",
            "m12_pareto_oncology.svg",
            figure_path_prefix,
            "Trade-off frontier across the full Option 2 eligible pool and the full HydroSDT-filtered GA pool, showing tumor-targeted cavitation intensity versus healthy-cell protection. Faint background points show the explored background population used for analysis, while the highlighted selected designs and frontier line identify where the shortlist sits inside that broader design landscape.",
            "Selected design oncology frontier",
        ),
        spec(
            11,
            "Pediatric Circuit Volume Margin",
            "m12_pediatric_ecv_margin.svg",
            figure_path_prefix,
            "Selected-design extracorporeal circuit volume as a percentage of the 3 kg neonatal 10% circuit-volume limit (25.5 mL). Lower is better; values below 100% satisfy the pediatric reference margin while preserving the selected treatment topology.",
            "Pediatric ECV margin bars",
        ),
    ];
    if !input.ga_best_per_gen.is_empty() {
        specs.push(spec(
            12,
            "GA Fitness Convergence",
            "m12_ga_convergence.svg",
            figure_path_prefix,
            "HydroSDT GA best-fitness trajectory across generations. The figure annotates the best-score change over the trailing generations so readers can distinguish a true plateau from a run that is still improving late through mutation and crossover.",
            "GA convergence line",
        ));
    }
    specs.push(spec(
        13,
        "Dean Number vs Cavitation at Serpentine Venturi Sites",
        "m12_dean_venturi_placement.svg",
        figure_path_prefix,
        "Per-venturi-placement comparison of Dean number (secondary flow intensity at the bend apex) and cavitation strength (1 - sigma) for the top GA-optimized serpentine design. Higher Dean numbers at bend apices pre-focus CTCs toward the outer wall via centrifugal secondary flow before they enter the venturi constriction, co-localizing inertial focusing with hydrodynamic cavitation. Throat velocity is annotated below each placement.",
        "Dean venturi placement dual-axis bars",
    ));
    specs.push(spec(
        14,
        "Treatment-Lane Geometry Zoom - Option 2 vs GA",
        "m12_treatment_lane_zoom.svg",
        figure_path_prefix,
        "Side-by-side treatment-lane zoom comparing the best serpentine-capable Option 2 venturi path with the top serpentine-capable GA-refined path. Purple polylines denote the treatment channel, orange segments denote venturi throats, and the panel framing is normalized per design so local serpentine curvature, throat clustering, and path compactness can be compared without the full-device split-tree dominating the visual. Panel headers report residence time, cavitation intensity, and local bend-radius callouts, while the bottom delta box summarizes how the GA changed ultrasound exposure residence and cavitation delivery relative to the deterministic serpentine Option 2 baseline.",
        "Treatment lane zoom comparison",
    ));
    Ok(specs)
}

fn spec(
    number: usize,
    title: &str,
    file_name: &str,
    figure_path_prefix: &str,
    caption: &str,
    alt: &str,
) -> NarrativeFigureSpec {
    NarrativeFigureSpec {
        number,
        title: title.to_string(),
        path: format!("{figure_path_prefix}/{file_name}"),
        caption: caption.to_string(),
        alt: alt.to_string(),
        is_placeholder: false,
    }
}

fn spec_placeholder(
    number: usize,
    title: &str,
    file_name: &str,
    figure_path_prefix: &str,
    caption: &str,
    alt: &str,
) -> NarrativeFigureSpec {
    NarrativeFigureSpec {
        number,
        title: title.to_string(),
        path: format!("{figure_path_prefix}/{file_name}"),
        caption: format!("{caption} *(Placeholder - source data not available in this run)*"),
        alt: alt.to_string(),
        is_placeholder: true,
    }
}

/// Returns `true` if a placeholder was generated (file did not already exist).
fn ensure_existing_or_placeholder(
    path: &Path,
    title: &str,
) -> Result<bool, Box<dyn std::error::Error>> {
    if path.exists() {
        Ok(false)
    } else {
        write_placeholder(
            path,
            title,
            "Source schematic not present in current run; placeholder generated.",
        )?;
        Ok(true)
    }
}

fn selected_schematic_caption(
    ranked: &Milestone12ReportDesign,
    option_label: &str,
    treatment_summary: &str,
) -> String {
    let candidate = &ranked.candidate;
    let ecv_ml = ranked.metrics.total_ecv_ml;
    let pediatric_limit_ml = pediatric_reference_ecv_limit_ml();
    let ecv_pct = 100.0 * ecv_ml / pediatric_limit_ml.max(1e-12);
    format!(
        "{} selected design: {}. Topology: {}. Visible split layers: {} ({}). ECV = {:.3} mL ({:.1}% of 3 kg neonatal 10% circuit-volume limit = {:.1} mL){}. Line thickness is proportional to channel width, and geometry-authored serpentines depict mirrored Dean-generating curvature rather than a single apex. Candidate: `{}`.",
        option_label,
        treatment_summary,
        ranked.topology_display_name(),
        visible_split_layers(ranked),
        stage_sequence_label(ranked),
        ecv_ml,
        ecv_pct,
        pediatric_limit_ml,
        venturi_caption_suffix(ranked),
        candidate.id,
    )
}

fn first_ranked(designs: &[Milestone12ReportDesign]) -> &Milestone12ReportDesign {
    &designs[0]
}

fn pediatric_reference_ecv_limit_ml() -> f64 {
    PEDIATRIC_REFERENCE_WEIGHT_KG * PEDIATRIC_BLOOD_VOLUME_ML_PER_KG * 0.10
}

fn pick_serpentine_venturi_focus(
    designs: &[Milestone12ReportDesign],
) -> Option<&Milestone12ReportDesign> {
    designs.iter().find(|design| {
        design.metrics.active_venturi_throat_count > 0
            && design
                .candidate
                .blueprint()
                .channels
                .iter()
                .any(|channel| {
                    matches!(
                        channel.visual_role,
                        Some(ChannelVisualRole::CenterTreatment | ChannelVisualRole::VenturiThroat)
                    ) && matches!(channel.channel_shape, ChannelShape::Serpentine { .. })
                })
    })
}

fn venturi_caption_suffix(ranked: &Milestone12ReportDesign) -> String {
    format!(
        ". Treatment mode: {}. Active venturi throats: {}",
        ranked.metrics.treatment_zone_mode, ranked.metrics.active_venturi_throat_count
    )
}

fn stage_sequence_label(ranked: &Milestone12ReportDesign) -> String {
    ranked.stage_sequence_label()
}

fn visible_split_layers(ranked: &Milestone12ReportDesign) -> usize {
    ranked.visible_split_layers()
}

/// Extract per-venturi-placement Dean number and cavitation data from a design.
///
/// The 1D lumped-element solver models the serpentine treatment channel as a
/// single resistance element, so all venturi placements resolve to the same
/// aggregate flow and pressure.  To produce per-position variation this
/// function applies two physical corrections:
///
/// 1. **Serial pressure decay**: each successive venturi throat consumes a
///    fraction of the available pressure budget.  The total treatment-channel
///    pressure drop (from_pressure - to_pressure) is distributed linearly
///    across N placements, so throat k sees upstream_pressure reduced by
///    k/N of the total drop.  This raises sigma (weakens cavitation) at
///    downstream positions.
///
/// 2. **Mirrored curvature variation**: real serpentines alternate between
///    inner (tighter) and outer (wider) bends.  The curvature radius at
///    bend k is modulated as R_k = R_base * (1 +/- 0.25 * cos(pi * k))
///    for mirrored U-turn geometry, producing alternating higher/lower
///    Dean numbers.
fn extract_dean_venturi_points(design: Option<&Milestone12ReportDesign>) -> Vec<DeanVenturiPoint> {
    use cfd_1d::{evaluate_venturi_screening, VenturiScreeningInput};

    let design = match design {
        Some(d) => d,
        None => return Vec::new(),
    };
    let candidate = &design.candidate;
    let topology = match candidate.topology_spec() {
        Ok(spec) => spec,
        Err(_) => return Vec::new(),
    };
    let n_placements = topology.venturi_placements.len();
    if n_placements == 0 {
        return Vec::new();
    }

    let solve = match crate::metrics::solve_blueprint_candidate(candidate) {
        Ok(s) => s,
        Err(_) => return Vec::new(),
    };

    // Find any venturi-carrying channel sample for base flow/pressure data.
    let base_sample = solve
        .channel_samples
        .iter()
        .find(|s| s.is_venturi_channel)
        .or_else(|| {
            let first_target = &topology.venturi_placements[0].target_channel_id;
            solve
                .channel_samples
                .iter()
                .find(|s| s.id.starts_with(first_target.as_str()))
        })
        .or_else(|| solve.channel_samples.iter().find(|s| s.flow_m3_s.abs() > 1e-18));
    let base_sample = match base_sample {
        Some(s) => s,
        None => return Vec::new(),
    };

    let flow_m3_s = base_sample.flow_m3_s.abs();
    let base_upstream_pa = candidate.operating_point.absolute_inlet_pressure_pa()
        + base_sample.from_pressure_pa.max(0.0);
    // Estimate per-placement pressure consumption.  The total device pressure
    // drop is split proportionally across the treatment channel venturi stages.
    // Each serial throat + its approach/diffuser consumes roughly equal pressure.
    let total_dp_pa = design.metrics.total_pressure_drop_pa.max(0.0);
    // Treatment-path share of total drop (treatment channels carry the venturi
    // constrictions, which dominate the resistance).  Use 70% as a conservative
    // estimate for venturi-dominated designs.
    let treatment_dp_pa = total_dp_pa * 0.70;

    // Base serpentine bend radius from the topology spec.
    let base_bend_radius_m = topology
        .venturi_placements
        .first()
        .and_then(|p| topology.channel_route(&p.target_channel_id))
        .and_then(|route| route.serpentine.as_ref())
        .map(|s| s.bend_radius_m)
        .or_else(|| {
            candidate
                .topology_spec()
                .ok()
                .into_iter()
                .flat_map(|spec| spec.treatment_channel_ids())
                .filter_map(|channel_id| topology.channel_route(&channel_id))
                .find_map(|route| route.serpentine.as_ref().map(|s| s.bend_radius_m))
        })
        .unwrap_or(1.5e-3);

    let kinematic_viscosity = BLOOD_VISCOSITY_PA_S / BLOOD_DENSITY_KG_M3;

    let mut points = Vec::with_capacity(n_placements);
    // Track cumulative state across serial throats:
    // - pressure consumed by each preceding throat
    // - nuclei generated by cavitation and partially dissolved in transit
    let mut cumulative_pressure_consumed = 0.0_f64;
    let mut upstream_nuclei = 0.0_f64;

    for (idx, placement) in topology.venturi_placements.iter().enumerate() {
        // Each preceding throat consumed its own pressure drop (Bernoulli +
        // friction + incomplete diffuser recovery).  This is much larger than
        // the simple linear fraction because the throat velocity head is high.
        let local_upstream_pa = base_upstream_pa - cumulative_pressure_consumed;

        // Mirrored curvature: alternating tighter/wider bends in a
        // serpentine U-turn pattern.  The +/-25% modulation is typical
        // for mirrored millifluidic serpentines where inner bends have
        // smaller radius than outer bends.
        let curvature_modulation = 1.0 + 0.25 * (std::f64::consts::PI * idx as f64).cos();
        let local_bend_radius_m = base_bend_radius_m * curvature_modulation;

        // Venturi screening with position-corrected pressure.
        let area_inlet = (placement.throat_geometry.inlet_width_m
            * placement.throat_geometry.throat_height_m)
            .max(1e-18);
        let area_throat = (placement.throat_geometry.throat_width_m
            * placement.throat_geometry.throat_height_m)
            .max(1e-18);
        let upstream_vel = flow_m3_s / area_inlet;
        let throat_vel = flow_m3_s / area_throat;
        let screening = evaluate_venturi_screening(VenturiScreeningInput {
            upstream_pressure_pa: local_upstream_pa,
            upstream_velocity_m_s: upstream_vel,
            throat_velocity_m_s: throat_vel,
            throat_hydraulic_diameter_m: 2.0
                * placement.throat_geometry.throat_width_m
                * placement.throat_geometry.throat_height_m
                / (placement.throat_geometry.throat_width_m
                    + placement.throat_geometry.throat_height_m)
                    .max(1e-18),
            throat_length_m: placement.throat_geometry.throat_length_m,
            density_kg_m3: BLOOD_DENSITY_KG_M3,
            viscosity_pa_s: BLOOD_VISCOSITY_PA_S,
            vapor_pressure_pa: BLOOD_VAPOR_PRESSURE_PA,
            vena_contracta_coeff: VENTURI_CC,
            diffuser_recovery_coeff: DIFFUSER_DISCHARGE_COEFF,
            upstream_nuclei_fraction: upstream_nuclei,
        });

        // Accumulate pressure consumed by this throat.
        // Each venturi consumes: bernoulli_drop + friction_drop - diffuser_recovery.
        // The net consumption is what reduces the available pressure for downstream throats.
        let net_throat_dp = screening.bernoulli_drop_pa
            + screening.throat_friction_drop_pa
            - screening.diffuser_recovery_pa;
        cumulative_pressure_consumed += net_throat_dp.max(0.0);

        // Propagate nuclei: the outlet nuclei fraction from this throat
        // becomes the upstream for the next, after partial dissolution
        // during the inter-throat transit.
        let inter_throat_transit_s = if throat_vel > 1e-6 {
            placement.throat_geometry.throat_length_m * 3.0 / throat_vel // ~3x throat length gap
        } else {
            0.01
        };
        let transport = cfd_core::physics::cavitation::nuclei_transport::NucleiTransport::new(
            cfd_core::physics::cavitation::nuclei_transport::NucleiTransportConfig::default(),
        );
        upstream_nuclei = transport.advect_1d_dissolution(
            screening.outlet_nuclei_fraction,
            inter_throat_transit_s,
        );

        // Dean number with position-dependent curvature.
        let d_h = base_sample.cross_section.hydraulic_diameter();
        let velocity = flow_m3_s / base_sample.cross_section.area().max(1e-18);
        let reynolds = velocity * d_h / kinematic_viscosity;
        let dean_number = if local_bend_radius_m > 0.0 {
            reynolds * (d_h / (2.0 * local_bend_radius_m)).sqrt()
        } else {
            0.0
        };

        let bend_type = if idx % 2 == 0 { "outer" } else { "inner" };
        points.push(DeanVenturiPoint {
            label: format!("Bend {} ({})", idx + 1, bend_type),
            dean_number,
            cavitation_number: screening.cavitation_number,
            throat_velocity_m_s: screening.effective_throat_velocity_m_s,
            upstream_pressure_kpa: local_upstream_pa * 1e-3,
            bend_radius_mm: local_bend_radius_m * 1e3,
        });
    }
    points
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::domain::fixtures::{operating_point, stage0_venturi_candidate};
    use crate::reporting::{compute_blueprint_report_metrics, Milestone12ReportDesign};
    use cfd_schematics::topology::VenturiPlacementMode;

    #[test]
    fn dean_venturi_points_from_serpentine_candidate() {
        // Higher flow + pressure for more realistic cavitation regime
        let op = operating_point(5.0e-6, 400_000.0, 0.45);
        let candidate =
            stage0_venturi_candidate("dean-test", op, VenturiPlacementMode::CurvaturePeakDeanNumber);
        let metrics = compute_blueprint_report_metrics(&candidate).expect("metrics should compute");
        let design = Milestone12ReportDesign::new(1, candidate, metrics, 0.5);
        let points = extract_dean_venturi_points(Some(&design));

        // The fixture has 2 venturi placements on a 6-segment serpentine.
        assert!(
            !points.is_empty(),
            "should produce at least one Dean-venturi data point"
        );
        for (i, pt) in points.iter().enumerate() {
            assert!(
                pt.dean_number >= 0.0,
                "Bend {} Dean number should be non-negative, got {}",
                i + 1,
                pt.dean_number
            );
            assert!(
                pt.cavitation_number.is_finite(),
                "Bend {} cavitation number should be finite",
                i + 1
            );
            assert!(
                pt.throat_velocity_m_s > 0.0,
                "Bend {} throat velocity should be positive, got {}",
                i + 1,
                pt.throat_velocity_m_s
            );
            eprintln!(
                "  Bend {}: De={:.1}, sigma={:.4}, v_throat={:.2} m/s",
                i + 1,
                pt.dean_number,
                pt.cavitation_number,
                pt.throat_velocity_m_s
            );
        }
    }
}
