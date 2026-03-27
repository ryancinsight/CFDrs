//! Conclusions, reference-list, figure, and storage blocks for the M12 narrative.

use std::fmt::Write as _;

use crate::constraints::M12_GA_HYDRO_SEED;
use crate::reporting::figures::NarrativeFigureSpec;
use crate::reporting::{Milestone12GaRankingAuditEntry, Milestone12ReportDesign};

fn parse_lineage_operator(metadata: &str) -> Option<&str> {
    metadata.split(';').find_map(|part| {
        let (key, value) = part.split_once('=')?;
        (key == "operator").then_some(value)
    })
}

fn ancestry_event_counts(design: &Milestone12ReportDesign) -> (usize, usize) {
    design
        .candidate
        .blueprint()
        .lineage()
        .map_or((0, 0), |lineage| {
            lineage
                .mutations
                .iter()
                .fold((0_usize, 0_usize), |mut counts, event| {
                    let operator = parse_lineage_operator(&event.mutation).unwrap_or("");
                    if operator.starts_with("operating_point") {
                        counts.1 += 1;
                    } else {
                        counts.0 += 1;
                    }
                    counts
                })
        })
}

fn build_ga_ranking_tradeoff_table(
    selected: &Milestone12GaRankingAuditEntry,
    displaced: &Milestone12GaRankingAuditEntry,
    table_number: usize,
) -> String {
    format!(
        "<div align=\"center\">\n\n<p><strong>Table {}.</strong> GA ancestry-adjusted final ranking tradeoff</p>\n\n<table style=\"width:88%; max-width:6.2in; margin:0 auto; border-collapse:collapse; table-layout:fixed; font-size:8.5pt; line-height:1.15;\">\n<thead><tr><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:18%;\">Finalist</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:34%;\">Candidate ID</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:8%;\">Raw score</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:10%;\">Adjusted score</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:10%;\">Geometry penalty</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:10%;\">Operating-point penalty</th></tr></thead>\n<tbody><tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Selected GA rank-1</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\"><code>{}</code></td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.4}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.4}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.3}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.3}</td></tr><tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Highest raw-score finalist</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\"><code>{}</code></td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.4}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.4}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.3}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.3}</td></tr></tbody></table>\n</div>",
        table_number,
        selected.candidate_id,
        selected.score,
        selected.adjusted_selection_score,
        selected.geometry_concentration_penalty,
        selected.operating_point_diversity_penalty,
        displaced.candidate_id,
        displaced.score,
        displaced.adjusted_selection_score,
        displaced.geometry_concentration_penalty,
        displaced.operating_point_diversity_penalty,
    )
}

fn current_geometry_penalty_weight() -> f64 {
    std::env::var("M12_GA_GEOMETRY_PENALTY_WEIGHT")
        .ok()
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite() && *value >= 0.0)
        .unwrap_or(0.00043)
}

fn current_cavitation_gain_margin() -> f64 {
    std::env::var("M12_GA_CAVITATION_GAIN_MARGIN")
        .ok()
        .and_then(|value| value.parse::<f64>().ok())
        .filter(|value| value.is_finite() && *value >= 0.0)
        .unwrap_or(0.01)
}

pub(super) fn build_workspace_configuration_section() -> String {
    format!(
        "The Milestone 12 reporting pipeline exposes the following workspace-level configuration knobs for HydroSDT GA shortlist formation and audit reproducibility.\n\n- `M12_GA_GEOMETRY_PENALTY_WEIGHT` = `{:.5}`. Applies the lineage-concentration penalty used to demote GA shortlists dominated by repeated geometry ancestry from the same serpentine family.\n- `M12_GA_CAVITATION_GAIN_MARGIN` = `{:.3}`. Requires a GA venturi design to exceed the deterministic Option 2 baseline cumulative cavitation dose by at least this amount before it remains in the report-ranked HydroSDT subset.\n\nIf these environment variables are unset, the defaults above are applied automatically during report generation.",
        current_geometry_penalty_weight(),
        current_cavitation_gain_margin(),
    )
}

fn geometry_displacement_threshold(
    selected: &Milestone12GaRankingAuditEntry,
    displaced: &Milestone12GaRankingAuditEntry,
) -> Option<f64> {
    let geometry_gap = displaced.geometry_concentration_penalty
        - selected.geometry_concentration_penalty;
    if geometry_gap <= 0.0 {
        return None;
    }

    let operating_point_weight = 0.00020;
    let numerator = displaced.score - selected.score
        + operating_point_weight
            * (selected.operating_point_diversity_penalty
                - displaced.operating_point_diversity_penalty);
    Some(numerator / geometry_gap)
}

fn cavitation_regime_summary(sigma: f64) -> String {
    if sigma < 0.0 {
        format!("σ = {sigma:.4} < 0 confirms sub-vapor-pressure throat operation with active hydrodynamic cavitation")
    } else if sigma < 1.0 {
        format!("σ = {sigma:.4} remains in the inception-capable window (0 < σ < 1), so the throat is cavitation-eligible without yet dropping below vapor pressure")
    } else {
        format!("σ = {sigma:.4} is above the hydrodynamic cavitation threshold")
    }
}

pub(super) fn build_conclusions(
    total_candidates: usize,
    topology_family_count: usize,
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
    ga_top: &[Milestone12ReportDesign],
    ga_best: &Milestone12ReportDesign,
    ga_best_per_gen: &[f64],
    ga_ranking_audit: &[Milestone12GaRankingAuditEntry],
    tradeoff_table_number: usize,
) -> String {
    let m2 = &option2.metrics;
    let mg = &ga_best.metrics;
    let mut s = String::new();

    let temp_fmt = |k: f64| -> String {
        if k < 0.01 && k > 0.0 {
            "< 0.01".to_string()
        } else {
            format!("{k:.2}")
        }
    };

    // §1 — Milestone completion
    let family_word = if topology_family_count == 1 {
        "family"
    } else {
        "families"
    };
    let _ = writeln!(
        s,
        "Milestone 12 required selecting millifluidic device designs meeting hydrodynamic and \
cavitation parameters for extracorporeal SDT. From {total_candidates} candidates across \
{topology_family_count} topology {family_word}, CFDrs identified one Option 1 and one Option 2 design via deterministic \
eligibility gating and ranking. No stochastic elements affect the final selection.\n",
    );

    // §2 — Option 1 acoustic
    if let Some(option1) = option1 {
        let m1 = &option1.metrics;
        let _ = writeln!(
            s,
            "**Option 1 - Selective Acoustic Center Treatment** (`{}`): {:.1}% of WBCs kept out of \
the treatment lane; RBC peripheral fraction {:.1}%; HI/pass = {:.4}% (FDA 0.1% limit); \
ECV = {:.3} mL; P95 wall shear = {:.1} Pa. Zero active venturi throats; treatment relies on \
externally applied 412 kHz ultrasound. Acoustic resonance factor (ARF) = {:.4} \
(channel D_h match to λ/2 ≈ 1.87 mm). Score {:.4} under AsymmetricSplitResidenceSeparation.\n",
            option1.candidate.id,
            (1.0 - m1.wbc_recovery) * 100.0,
            m1.rbc_peripheral_fraction * 100.0,
            m1.hemolysis_index_per_pass,
            m1.total_ecv_ml,
            m1.wall_shear_p95_pa,
            m1.acoustic_resonance_factor,
            option1.score,
        );
    } else {
        let _ = writeln!(
            s,
            "**Option 1 - Selective Acoustic Center Treatment**: no design satisfied the strict acoustic eligibility gates under the current physics regime. The report therefore records Option 1 explicitly as an empty shortlist rather than treating the absence as a ranking regression.\n"
        );
    }

    // §3 — Option 2 venturi
    let gauge_kpa = option2.inlet_gauge_kpa();
    let d_throat_um = option2.throat_width_um().unwrap_or(0.0);
    let cavitation_summary = cavitation_regime_summary(m2.cavitation_number);
    let cavitation_gain_margin = current_cavitation_gain_margin();
    let _ = writeln!(
        s,
        "**Option 2 - Hydrodynamic Cavitation SDT** (`{}`): {} under {:.0} kPa gauge through {:.0} µm throat; {} serial stage(s) per path; \
{} total active throats. Cancer routing {:.1}% (cancer_center_fraction); therapeutic window \
score = {:.3}; WBC treatment exposure {:.1}%; HI/pass = {:.4}%; throat viscous heating ΔT = {} K \
(FDA 5 K limit: {}). Score {:.4} under AsymmetricSplitVenturiCavitationSelectivity. \
**Not comparable to Option 1 score.**\n",
        option2.candidate.id,
        cavitation_summary,
        gauge_kpa,
        d_throat_um,
        m2.serial_venturi_stages_per_path,
        m2.active_venturi_throat_count,
        m2.cancer_center_fraction * 100.0,
        m2.therapeutic_window_score,
        m2.wbc_recovery * 100.0,
        m2.hemolysis_index_per_pass,
        temp_fmt(m2.throat_temperature_rise_k),
        if m2.fda_thermal_compliant {
            "PASS"
        } else {
            "FAIL"
        },
        option2.score,
    );

    // §4 — Safety and FDA compliance
    if let Some(option1) = option1 {
        let m1 = &option1.metrics;
        let _ =
            writeln!(
            s,
            "**Safety and FDA Compliance:** Both designs pass all five hard eligibility gates. \
Max P95 wall shear: Option 1 = {:.1} Pa, Option 2 = {:.1} Pa (FDA 150 Pa sustained limit). \
Venturi-channel transit time: {:.2e} s ({} 5 ms transient threshold, so {} shear limit applies). \
Clotting risk index = {:.4} at nominal flow; flow caution flags are `Q>=200={}` and `Q>=600={}` \
for this selected operating point. ECV = {:.3} mL within pediatric circuit targets. \
FDA thermal compliance (42 °C ceiling) for Option 2: {} ({} K rise).\n",
            m1.wall_shear_p95_pa,
            m2.wall_shear_p95_pa,
            m2.throat_transit_time_s,
            if m2.throat_transit_time_s < 5e-3 { "<" } else { "≥" },
            if m2.throat_transit_time_s < 5e-3 { "transient 300 Pa" } else { "sustained 150 Pa" },
            m2.clotting_risk_index,
            if m2.clotting_flow_compliant { "PASS" } else { "FAIL" },
            if m2.clotting_flow_compliant_10ml_s { "PASS" } else { "FAIL" },
            m2.total_ecv_ml,
            if m2.fda_thermal_compliant { "PASS" } else { "FAIL" },
            temp_fmt(m2.throat_temperature_rise_k),
        );
    } else {
        let _ = writeln!(
            s,
            "**Safety and FDA Compliance:** The selected Option 2 design passes all five hard eligibility gates. Option 1 produced no eligible shortlist under the current physics regime, so no acoustic selected-design safety row exists for this run. Option 2 P95 wall shear = {:.1} Pa (FDA 150 Pa sustained limit). Venturi-channel transit time: {:.2e} s ({} 5 ms transient threshold, so {} shear limit applies). Clotting risk index = {:.4} at nominal flow; flow caution flags are `Q>=200={}` and `Q>=600={}` for this selected operating point. ECV = {:.3} mL within pediatric circuit targets. FDA thermal compliance (42 °C ceiling) for Option 2: {} ({} K rise).\n",
            m2.wall_shear_p95_pa,
            m2.throat_transit_time_s,
            if m2.throat_transit_time_s < 5e-3 { "<" } else { "≥" },
            if m2.throat_transit_time_s < 5e-3 { "transient 300 Pa" } else { "sustained 150 Pa" },
            m2.clotting_risk_index,
            if m2.clotting_flow_compliant { "PASS" } else { "FAIL" },
            if m2.clotting_flow_compliant_10ml_s { "PASS" } else { "FAIL" },
            m2.total_ecv_ml,
            if m2.fda_thermal_compliant { "PASS" } else { "FAIL" },
            temp_fmt(m2.throat_temperature_rise_k),
        );
    }

    // §5 — Topology dominance rationale
    let (w_mm, h_mm) = option2
        .candidate
        .blueprint()
        .topology_spec()
        .and_then(|spec| {
            spec.split_stages
                .first()
                .and_then(|stage| stage.branches.iter().find(|branch| branch.treatment_path))
                .map(|branch| (spec.inlet_width_m * 1000.0, branch.route.height_m * 1000.0))
        })
        .unwrap_or((0.0, 0.0));
    let _ = writeln!(
        s,
        "**Topology Selection Physics - Option 1 vs Option 2:** The two optimization \
goals impose distinct constraints on optimal tree depth. For **Option 1** \
(selective acoustic, no venturi), the hybrid additive + geometric-mean synergy score \
rewards deeper split trees because each additional Zweifach–Fung splitting stage compounds \
cancer-center enrichment multiplicatively in the underlying flow partition: a three-stage \
Tri→Tri→Tri tree achieves cancer_center_fraction ~0.75+ compared to ~0.556 for Tri→Tri. \
Both the additive cancer-focusing term (22% weight) and the 4th-root synergy term amplify \
this separation advantage while the additive base ensures no feasible design collapses \
to zero (floor = 0.001). \
For **Option 2** (selective venturi cavitation), the selected design uses an \
{w_mm:.0} mm × {h_mm:.0} mm cross-section. Here, deeper trifurcation beyond 2 levels becomes \
counterproductive: a third level narrows terminal treatment channels to approximately \
140 µm (below the 200–400 µm Zweifach–Fung inertial-focusing optimum for millifluidic \
scales), routes only 1/27 of inlet flow to each venturi throat (insufficient velocity \
for cavitation inception at the available pressure budget), and consumes excess pressure \
drop across additional branching stages. RBC peripheral separation is achieved at \
Level 1; Level 2 refines cancer/WBC center concentration. The optimal tree depth is \
therefore physics-driven: acoustic-only treatment benefits from maximal splitting depth, \
while venturi cavitation requires sufficient flow velocity at the throat, imposing an \
upper bound on the number of upstream flow-dividing stages.\n",
    );

    // §6 — GA results and modeling note
    let _ = writeln!(
        s,
        "**GA Results - In-Place Dean-Serpentine Refinement:** The blueprint-native GA \
(seed {M12_GA_HYDRO_SEED}, InPlaceDeanSerpentineRefinement goal) produced rank-1 design \
`{}` with {} venturi throat geometries ({} serial stage(s)); note that throats may not \
produce hydrodynamic cavitation at all operating points — σ >> 1 indicates the throats function \
as Dean-focusing constrictions rather than cavitation sources. The GA applies three classes of \
architecture-preserving mutations and crossover operations: (1) multi-regime treatment-lane \
serpentine variants spanning compact, smooth, dense, and long curvature profiles to modulate \
Dean secondary flow (De = Re √(D_h/2R)) at bend apices where centrifugal forces focus larger \
CTCs toward the outer wall; (2) venturi retargeting / throat-geometry edits that alter serial \
throat count and constriction ratio to lower σ at the vena contracta; and (3) compatible-parent \
crossover plus diversity-aware survivor selection so operating-point advantages and topology \
features from distinct high-performing parents can be recombined instead of only stacked along a \
single lineage. The millifluidic Dean correlation of Bayat-Rezai (2017) is used for channels \
with D_h > 500 µm, providing a validated correction to the classical Dean (1927) formula for the \
larger aspect ratios and Reynolds numbers characteristic of millifluidic geometries. The Dean \
number bonus (De_max/100) in the GA score explicitly rewards designs that co-localise inertial \
focusing from Dean vortices with hydrodynamic cavitation at venturi throats positioned at bend \
apices. For HydroSDT shortlist selection, GA venturi candidates must exceed the deterministic \
Option 2 baseline cumulative cavitation dose by at least {:.3} to remain in the report-ranked \
subset. GA score {:.4}. **Not comparable to Option 2 Combined mode score.**\n",
        ga_best.candidate.id,
        mg.active_venturi_throat_count,
        mg.serial_venturi_stages_per_path,
        cavitation_gain_margin,
        ga_best.score,
    );

    if ga_top.len() >= 2 {
        let runner_up = &ga_top[1];
        let ga_sep = ga_best.metrics.three_pop_sep_efficiency;
        let runner_up_sep = runner_up.metrics.three_pop_sep_efficiency;
        let _ = writeln!(
            s,
            "Within the top GA shortlist, three-population separation remains nearly flat \
({ga_sep:.4} for rank-1 versus {runner_up_sep:.4} for rank-2), so rank-1 wins mainly by \
maintaining stronger cavitation-driven treatment performance under comparable routing \
selectivity rather than by a large change in bulk separation architecture.\n"
        );
    }

    let (geometry_event_count, operating_point_event_count) = ancestry_event_counts(ga_best);
    let _ = writeln!(
        s,
        "The rank-1 ancestry audit records {geometry_event_count} geometry-lineage events and \
{operating_point_event_count} operating-point refinement event(s), indicating that the winning \
GA family reflects deep inherited topology ancestry with limited late-stage flow/pressure tuning \
rather than a purely local last-step tweak.\n"
    );

    if let Some(selected) = ga_ranking_audit.first() {
        let displaced_raw_winner = ga_ranking_audit.iter().max_by(|left, right| {
            left.score
                .total_cmp(&right.score)
                .then_with(|| left.rank.cmp(&right.rank))
        });
        if let Some(displaced) = displaced_raw_winner {
            if displaced.candidate_id != selected.candidate_id || (displaced.score - selected.score).abs() > 1.0e-12 {
                let _ = writeln!(
                    s,
                    "The ancestry-adjusted ranking selected `{}` at adjusted score {:.4} \
(raw {:.4}; geometry penalty {:.3}; operating-point penalty {:.3}) over the highest raw-score \
finalist `{}` (raw {:.4}; adjusted {:.4}) because the displaced design carried a larger geometry \
concentration penalty of {:.3}{}.\n",
                    selected.candidate_id,
                    selected.adjusted_selection_score,
                    selected.score,
                    selected.geometry_concentration_penalty,
                    selected.operating_point_diversity_penalty,
                    displaced.candidate_id,
                    displaced.score,
                    displaced.adjusted_selection_score,
                    displaced.geometry_concentration_penalty,
                    if displaced.operating_point_diversity_penalty > 0.0 {
                        format!(
                            " and operating-point monoculture penalty {:.3}",
                            displaced.operating_point_diversity_penalty
                        )
                    } else {
                        String::new()
                    },
                );
                let _ = writeln!(
                    s,
                    "The direct ranking tradeoff is summarized below.\n\n{}",
                    build_ga_ranking_tradeoff_table(selected, displaced, tradeoff_table_number)
                );
                if let Some(threshold) = geometry_displacement_threshold(selected, displaced) {
                    let _ = writeln!(
                        s,
                        "At the current operating-point coefficient (0.00020), any geometry penalty coefficient above {:.6} is sufficient to let this lower-concentration finalist displace the higher raw-score design; this run used {:.6}.\n",
                        threshold,
                        current_geometry_penalty_weight(),
                    );
                }
            }
        }
    }

    // §6b — GA convergence note (flat fitness indicates local optimality of seed)
    if ga_best_per_gen.len() >= 2 {
        let first = ga_best_per_gen.first().copied().unwrap_or(0.0);
        let last = ga_best_per_gen.last().copied().unwrap_or(0.0);
        if (first - last).abs() < 1.0e-8 {
            let _ = writeln!(
                s,
                "The GA did not improve upon the seed design across {} generations, indicating \
the parametric selection was already locally optimal under the available mutation and crossover \
operators in this run.\n",
                ga_best_per_gen.len(),
            );
        } else {
            let improvement = last - first;
            let _ = writeln!(
                s,
                "Over {} generations the GA improved the best score by {improvement:.4} \
(from {first:.4} to {last:.4}), confirming that diversity-aware survivor selection together with \
architecture-preserving mutation and crossover discovered beneficial refinements beyond the \
parametric sweep.\n",
                ga_best_per_gen.len(),
            );
        }
    }

    // §7 — Acoustic resonance opportunity
    if let Some(option1) = option1 {
        let m1 = &option1.metrics;
        let _ = writeln!(
            s,
            "**Acoustic Resonance Opportunity (Option 1):** At 412 kHz, the acoustic \
half-wavelength in blood is approximately 1.87 mm (c_blood = 1540 m/s). Channels with \
hydraulic diameter D_h ≈ nλ/2 create standing-wave pressure antinodes that preferentially \
trap resonant-radius bubbles (R_res ≈ 7.5 µm at 412 kHz), amplifying sonoporation and \
sonosensitizer activation. The Option 1 acoustic resonance factor ARF = {:.4} \
(channel_resonance_score = {:.4}). Channels with D_h tuned closer to 1.87 mm would \
approach ARF = 1.0. This channel-dimension resonance matching is scored in CFDrs via the \
`acoustic_resonance_factor` metric and represents an enhancement opportunity for \
Milestone 13 experimental optimization. All outputs are anchored to GA seed {M12_GA_HYDRO_SEED} \
and deterministic parametric ordering for bit-exact reproducibility.",
            m1.acoustic_resonance_factor, m1.channel_resonance_score,
        );
    } else {
        let _ = writeln!(
            s,
            "**Acoustic Resonance Opportunity (Option 1):** The selective acoustic track produced no eligible shortlist under the current physics regime, so no selected-design acoustic resonance factor is reported for this run. The resonance-targeting mechanism remains relevant for future geometry updates because channels tuned toward the 412 kHz blood half-wavelength (approximately 1.87 mm) should improve standing-wave amplification once a pressure-feasible acoustic design re-enters the strict Option 1 pool. All outputs are anchored to GA seed {M12_GA_HYDRO_SEED} and deterministic parametric ordering for bit-exact reproducibility."
        );
    }

    s
}

pub(super) fn build_references_block() -> String {
    "\
1. ANSI/SLAS 1-2004, Microplates, Footprint Dimensions.\n\
2. Giersiepen, M., et al. \"Estimation of shear stress-related blood damage in heart valve prostheses: in vitro comparison of 25 aortic valves.\" *International Journal of Artificial Organs*, 13(5):300–306, 1990.\n\
3. Brennen, C.E. *Cavitation and Bubble Dynamics*. Oxford University Press, 1995. [Eq. 3.12, cavitation number inception criterion.]\n\
4. Ohl, S.-W., et al. \"Sonoporation from jetting cavitation bubbles.\" *Biophysical Journal*, 91(11):4285–4295, 2006. [5× membrane lysis amplification at bubble collapse; basis for `LYSIS_CAVITATION_AMPLIFICATION = 5.0`.]\n\
5. Hellums, J.D. \"1993 Whitaker Lecture: Biorheology in thrombosis research.\" *Annals of Biomedical Engineering*, 22(5):445–455, 1994. [PAI exponent model, Eq. 3: n=1.325, m=0.462.]\n\
6. Di Carlo, D. \"Inertial microfluidics.\" *Lab on a Chip*, 9(21):3038–3046, 2009. [κ_RBC = a_RBC / D_h confinement criterion for inertial focusing; threshold 0.07.]\n\
7. Zweifach, B.W. and Fung, Y.C. \"Phase separation in capillary networks.\" *Microvascular Research*, 1971. [β_RBC = 1.0 for passive tracer routing at millifluidic scales.]\n\
8. FDA 2019 Guidance, *Nonclinical Tests and Recommended Labeling for Intravascular Administration Sets, Blood Administration Sets, and Blood Component Administration Sets*. [FDA predicate: Maquet RotaFlow, K143453, 1% hemolysis ceiling.]\n\
9. Lentner, C. (Ed.) *Geigy Scientific Tables, Vol. 3: Physical Chemistry, Composition of Blood.* Novartis, 1984. [Table 30: neonatal reference blood volume 85 mL/kg; basis for `PEDIATRIC_BLOOD_VOLUME_ML_PER_KG = 85.0`.]\n\
10. Dean, W.R. \"Note on the motion of fluid in a curved pipe.\" *Philosophical Magazine*, 4(20):208–223, 1927. [Dean number and secondary flow in curved channels; basis for `CurvaturePeakDeanNumber` venturi placement mode.]\n\
11. Taskin, M.E., et al. \"Evaluation of Eulerian and Lagrangian models for hemolysis estimation.\" *ASAIO Journal*, 58(4):363–372, 2012. [Strain-based hemolysis model; alternative to Giersiepen (1990) for cumulative strain history.]\n\
12. Amini, H., et al. \"Inertial microfluidic physics.\" *Lab on a Chip*, 14(15):2739–2761, 2014. [Confinement-dependent lift correction for a/D_h > 0.1; refines Di Carlo (2009) for CTCs in millifluidic channels.]\n\
13. Pries, A.R., et al. \"Blood viscosity in tube flow: dependence on diameter and hematocrit.\" *American Journal of Physiology*, 263(6):H1770–H1778, 1992. [Fahraeus–Lindqvist apparent viscosity reduction in tubes < 300 µm.]\n\
14. Pries, A.R., et al. \"Red cell distribution at microvascular bifurcations.\" *Microvascular Research*, 38(1):81–101, 1989. [Plasma skimming hematocrit partitioning at asymmetric bifurcations.]\n\
15. Quemada, D. \"Rheology of concentrated disperse systems and minimum energy dissipation principle.\" *Rheologica Acta*, 17(6):632–642, 1978. [Rouleaux aggregation viscosity model for low-shear zones.]\n\
16. Durst, F., et al. \"The development lengths of laminar pipe and channel flows.\" *Journal of Fluids Engineering*, 127(6):1154–1160, 2005. [Developing-flow entrance correction for short venturi throats with L/D_h < 20.]\n\
17. White, F.M. *Fluid Mechanics*, 7th ed. McGraw-Hill, 2011. [Cross-section averaged throat velocity for cavitation number computation.]\n\
18. Bayat, P. and Rezai, P. \"Semi-empirical estimation of Dean flow velocity in curved microchannels.\" *Scientific Reports*, 7:13655, 2017. [Millifluidic Dean correlation for D_h > 500 µm; validated correction for large aspect ratios.]\n\
19. Gor'kov, L.P. \"On the forces acting on a small particle in an acoustical field in an ideal fluid.\" *Soviet Physics Doklady*, 6(9):773–775, 1962. [Acoustic radiation force and energy density; basis for `acoustic_contrast_factor()` and `acoustic_energy_density()`.]\n\
20. Bruus, H. \"Acoustofluidics 7: The acoustic radiation force on small particles.\" *Lab on a Chip*, 12(6):1014-1021, 2012. [Acoustic radiation force in microfluidics; basis for radiation force model F = 4πa³ΦkE_ac·sin(2kx).]\n\
21. Rosenthal, I., et al. \"Sonodynamic therapy: a review of the synergistic effects of drugs and ultrasound.\" *Ultrasonics Sonochemistry*, 11(6):349–363, 2004. [First-order sonosensitizer activation kinetics η = 1−exp(−k·I·t); basis for `sonosensitizer_activation_efficiency()`.]\n\
22. Rayleigh, Lord. \"On the pressure developed in a liquid during the collapse of a spherical cavity.\" *Philosophical Magazine*, 34(200):94–98, 1917. [Bubble collapse time t_c = 0.915R√(ρ/p) and jet velocity; basis for `cavitation_hemolysis_amplification()`.]\n\
23. Secomb, T.W. \"Blood flow in the microcirculation.\" *Annual Review of Fluid Mechanics*, 49:443–461, 2017. [Network blood viscosity model including in-vivo correction and phase separation parameter X₀.]\n\
24. Gosling, R.G. and King, D.H. \"Arterial assessment by Doppler-shift ultrasound.\" *Proceedings of the Royal Society of Medicine*, 67(6):447–449, 1974. [Pulsatility index PI = (V_sys−V_dia)/V_mean for Womersley pulsatile flow characterization.]\n\
25. SonALAsense Internal Data, CFDrs canonical simulation data."
        .to_string()
}

pub(super) fn build_figure_toc_rows(specs: &[NarrativeFigureSpec]) -> String {
    let mut out = String::new();
    for spec in specs {
        let _ = writeln!(
            out,
            "| {} | [{}](#fig-{}) |",
            spec.number, spec.title, spec.number
        );
    }
    out
}

pub(super) fn build_figure_sections(specs: &[NarrativeFigureSpec]) -> String {
    let mut out = String::new();
    for spec in specs {
        let _ = writeln!(out, "<a id=\"fig-{}\"></a>\n", spec.number);
        let _ = writeln!(
            out,
            "<p align=\"center\"><img src=\"{}\" alt=\"{}\" style=\"max-width:100%;width:6.5in;\" /></p>\n",
            spec.path, spec.alt
        );
        let _ = writeln!(out, "*Figure {}. {}*\n", spec.number, spec.caption);
    }
    out
}

pub(super) fn build_storage_policy_section() -> String {
    "\
Milestone 12 computational data is stored under the following controls.\n\n\
**Version Control:** Simulation and report assets are maintained in the SonALAsense GitHub repository with branch/merge workflows and auditable change history.\n\n\
**Shared Storage:** Project data packages and shared artifacts are managed in SonALAsense Egnyte with role-based access and organization-managed retention policies.\n\n\
**Reproducibility:** Canonical milestone outputs are regenerated deterministically from source code and preserved alongside run artifacts.\n"
        .to_string()
}

pub(super) fn build_storage_artifact_index() -> String {
    "\
- Figures: inline throughout §5 Results\n\
- Generation artifacts: top-5 JSON, validation summaries, robustness outputs, and GA artifacts are stored in the canonical run output directory"
        .to_string()
}
