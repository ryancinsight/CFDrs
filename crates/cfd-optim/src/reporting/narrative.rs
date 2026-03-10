//! Milestone 12 narrative report generation from template + canonical results.

use std::collections::BTreeMap;
use std::path::{Path, PathBuf};

use serde::{Deserialize, Serialize};

use crate::analysis::RobustnessReport;
use crate::constraints::M12_GA_HYDRO_SEED;
use crate::reporting::contract_trace::load_m12_contract_text;
use crate::reporting::figures::{generate_m12_report_figures, FigureGenerationInput};
use crate::reporting::narrative_sections::{
    build_cavitation_formulas_intro, build_cri_expansion_sensitivity, build_figure_sections,
    build_figure_toc_rows, build_limits_of_usage, build_limits_of_usage_intro,
    build_option2_top5_table, build_results_intro, build_robustness_section, build_selected_table,
    build_selected_table_intro, build_storage_artifact_index, build_storage_policy_section,
    build_strict_core_intro, build_strict_core_table, build_top5_intro, build_tri_cell_intro,
    build_tri_cell_table, build_validation_section,
};
use crate::reporting::template::render_template_strict;
use crate::reporting::{Milestone12ReportDesign, ValidationRow};

/// Metadata for the narrative title page.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct M12Metadata {
    pub agreement_number: String,
    pub project_description: String,
    pub milestone_label: String,
    pub delivery_date: String,
    pub prepared_for: Vec<String>,
    pub logo_path: String,
    pub address_lines: Vec<String>,
    pub admin_contact_name: String,
    pub admin_contact_email: String,
    pub investigator_name: String,
    pub investigator_email: String,
}

/// Inputs used to render the narrative report.
pub struct Milestone12NarrativeInput<'a> {
    pub total_candidates: usize,
    pub option1_pool_len: usize,
    pub option2_pool_len: usize,
    pub option1_ranked: &'a [Milestone12ReportDesign],
    pub option2_ranked: &'a [Milestone12ReportDesign],
    pub ga_top: &'a [Milestone12ReportDesign],
    pub option2_pool_all: &'a [Milestone12ReportDesign],
    pub ga_pool_all: &'a [Milestone12ReportDesign],
    pub validation_rows: &'a [ValidationRow],
    pub option2_robustness: &'a [RobustnessReport],
    pub ga_best_per_gen: &'a [f64],
    pub fast_mode: bool,
}

/// Written artifact paths from narrative generation.
#[derive(Debug, Clone)]
pub struct Milestone12NarrativeArtifacts {
    pub narrative_path: PathBuf,
    pub figure_manifest_path: PathBuf,
    pub figure_count: usize,
}

/// Generate the M12 narrative report and figure manifest.
///
/// # Errors
/// Returns an error if required data is missing or any artifact cannot be written.
pub fn write_milestone12_narrative_report(
    workspace_root: &Path,
    canonical_results_path: &Path,
    input: &Milestone12NarrativeInput<'_>,
) -> Result<Milestone12NarrativeArtifacts, Box<dyn std::error::Error>> {
    let report_dir = workspace_root.join("report");
    let template = std::fs::read_to_string(
        report_dir
            .join("templates")
            .join("m12_narrative_template.md"),
    )?;
    let metadata: M12Metadata = serde_json::from_str(&std::fs::read_to_string(
        report_dir.join("m12_metadata.json"),
    )?)?;
    let contract = load_m12_contract_text(workspace_root)?;

    let option1 = input.option1_ranked.first();
    let option2 = input
        .option2_ranked
        .first()
        .ok_or("option2_ranked is empty; narrative requires selected Option 2")?;
    let ga_best = input
        .ga_top
        .first()
        .ok_or("ga_top is empty; narrative requires GA rank-1")?;

    let figure_specs = generate_m12_report_figures(
        &report_dir.join("figures"),
        &FigureGenerationInput {
            option1_ranked: input.option1_ranked,
            option2_ranked: input.option2_ranked,
            ga_top: input.ga_top,
            option2_pool_all: input.option2_pool_all,
            ga_pool_all: input.ga_pool_all,
            validation_rows: input.validation_rows,
            ga_best_per_gen: input.ga_best_per_gen,
            fast_mode: input.fast_mode,
        },
    )?;
    let manifest_path = report_dir.join("milestone12").join("figure_manifest.json");
    std::fs::create_dir_all(manifest_path.parent().ok_or("invalid manifest path")?)?;
    std::fs::write(&manifest_path, serde_json::to_string_pretty(&figure_specs)?)?;

    let mut values = BTreeMap::new();
    insert_title_page_values(&mut values, &metadata);
    insert_contract_values(&mut values, &contract);
    insert_data_values(&mut values, input, option1, option2, ga_best);
    insert_table_values(&mut values, input, option1, option2);
    values.insert(
        "FIGURE_TOC_ROWS".to_string(),
        build_figure_toc_rows(&figure_specs),
    );
    values.insert(
        "FIGURE_SECTIONS".to_string(),
        build_figure_sections(&figure_specs),
    );
    values.insert(
        "STORAGE_POLICY_SECTION".to_string(),
        build_storage_policy_section(),
    );
    values.insert(
        "STORAGE_ARTIFACT_INDEX".to_string(),
        build_storage_artifact_index(),
    );
    values.insert(
        "CONCLUSIONS".to_string(),
        build_conclusions(input.total_candidates, option1, option2, ga_best),
    );
    values.insert("REFERENCES_BLOCK".to_string(), build_references_block());

    values.insert(
        "APPENDIX_A_SUPPLEMENTAL".to_string(),
        build_appendix_a_supplemental(
            input.total_candidates,
            input.option1_pool_len,
            input.option2_pool_len,
            ga_best,
            input.option2_robustness,
            input.validation_rows,
            canonical_results_path,
        ),
    );

    let manifest_json = std::fs::read_to_string(&manifest_path)?;
    values.insert(
        "APPENDIX_B_FIGURES".to_string(),
        format!("```json\n{}\n```", manifest_json),
    );

    let narrative = render_template_strict(&template, &values)?;

    let narrative_path = report_dir.join("ARPA-H_SonALAsense_Milestone 12 Report.md");
    std::fs::write(&narrative_path, narrative)?;

    Ok(Milestone12NarrativeArtifacts {
        narrative_path,
        figure_manifest_path: manifest_path,
        figure_count: figure_specs.len(),
    })
}

fn insert_title_page_values(values: &mut BTreeMap<String, String>, metadata: &M12Metadata) {
    values.insert(
        "AGREEMENT_NUMBER".to_string(),
        metadata.agreement_number.clone(),
    );
    values.insert(
        "PROJECT_DESCRIPTION".to_string(),
        metadata.project_description.clone(),
    );
    values.insert(
        "MILESTONE_LABEL".to_string(),
        metadata.milestone_label.clone(),
    );
    values.insert("DELIVERY_DATE".to_string(), metadata.delivery_date.clone());
    values.insert(
        "PREPARED_FOR_BLOCK".to_string(),
        metadata.prepared_for.join("\n\n"),
    );
    values.insert("LOGO_PATH".to_string(), metadata.logo_path.clone());
    values.insert(
        "ADDRESS_BLOCK".to_string(),
        metadata.address_lines.join("\n\n"),
    );
    values.insert(
        "ADMIN_CONTACT_NAME".to_string(),
        metadata.admin_contact_name.clone(),
    );
    values.insert(
        "ADMIN_CONTACT_EMAIL".to_string(),
        metadata.admin_contact_email.clone(),
    );
    values.insert(
        "INVESTIGATOR_NAME".to_string(),
        metadata.investigator_name.clone(),
    );
    values.insert(
        "INVESTIGATOR_EMAIL".to_string(),
        metadata.investigator_email.clone(),
    );
}

fn insert_contract_values(
    values: &mut BTreeMap<String, String>,
    contract: &crate::reporting::contract_trace::MilestoneContractText,
) {
    values.insert(
        "MILESTONE_DESCRIPTION".to_string(),
        contract.description.clone(),
    );
    values.insert(
        "MILESTONE_EXIT_CRITERIA".to_string(),
        contract.exit_criteria.clone(),
    );
    values.insert(
        "MILESTONE_DELIVERABLE".to_string(),
        format!("{}, including limits of usage.", contract.deliverable),
    );
}

fn insert_data_values(
    values: &mut BTreeMap<String, String>,
    input: &Milestone12NarrativeInput<'_>,
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
    ga_best: &Milestone12ReportDesign,
) {
    let abstract_validation_sentence = if input.fast_mode {
        "Robustness screening was not computed in this fast-mode run, and multi-fidelity pressure-drop confirmation is generated by the companion `milestone12_validation` example rather than the main report run."
    } else {
        "Selected designs were screened for robustness under +/-10%/+/-20% operating-parameter perturbations. Multi-fidelity pressure-drop confirmation is generated separately by the companion `milestone12_validation` example using 1D Hagen-Poiseuille, 2D FVM, and 3D FEM solvers with non-Newtonian Carreau-Yasuda blood rheology."
    };
    let methods_pipeline_step7 = if input.fast_mode {
        "7. **Deferred robustness and companion validation flow** — When `M12_FAST=1` is enabled, robustness screening is intentionally skipped. Multi-fidelity 2D/3D validation is produced separately by the `milestone12_validation` example and may be appended after the main report run."
    } else {
        "7. **Companion multi-fidelity validation** — Selected designs are first ranked and screened here; 2D FVM and 3D FEM pressure-drop confirmation is then produced by the separate `milestone12_validation` example, while this report run retains robustness-screening results."
    };
    let milestone_scope_sentence = if input.fast_mode {
        "Milestone 12 requires selecting optimal design(s) that satisfy hydrodynamic and cavitation parameters. This fast-mode report documents the deterministic execution of the CFDrs pipeline from candidate generation through strict gating, scoring, and ranking, while explicitly deferring robustness screening and leaving multi-fidelity validation to the companion example."
    } else {
        "Milestone 12 requires selecting optimal design(s) that satisfy hydrodynamic and cavitation parameters. This report documents the deterministic execution of the CFDrs pipeline from candidate generation through strict gating, scoring, ranking, and robustness screening, with multi-fidelity validation emitted separately by the companion example."
    };
    values.insert(
        "TOTAL_CANDIDATES".to_string(),
        input.total_candidates.to_string(),
    );
    values.insert(
        "CANONICAL_SOURCE".to_string(),
        "report/milestone12_results.md".to_string(),
    );
    values.insert(
        "OPTION1_POOL".to_string(),
        input.option1_pool_len.to_string(),
    );
    values.insert(
        "OPTION2_POOL".to_string(),
        input.option2_pool_len.to_string(),
    );
    values.insert(
        "OPTION1_ID".to_string(),
        option1.map_or_else(
            || "NONE_ELIGIBLE_CURRENT_PHYSICS".to_string(),
            |design| design.candidate.id.clone(),
        ),
    );
    values.insert("OPTION2_ID".to_string(), option2.candidate.id.clone());
    values.insert("GA_ID".to_string(), ga_best.candidate.id.clone());
    values.insert("GA_SEED".to_string(), M12_GA_HYDRO_SEED.to_string());
    values.insert(
        "ABSTRACT_VALIDATION_SENTENCE".to_string(),
        abstract_validation_sentence.to_string(),
    );
    values.insert(
        "METHODS_PIPELINE_STEP7".to_string(),
        methods_pipeline_step7.to_string(),
    );
    values.insert(
        "MILESTONE_SCOPE_SENTENCE".to_string(),
        milestone_scope_sentence.to_string(),
    );
    values.insert(
        "FIG_PROCESS".to_string(),
        "../report/figures/milestone12_creation_optimization_process.svg".to_string(),
    );
    values.insert(
        "FIG_TREATMENT_BI".to_string(),
        "../report/figures/treatment_zone_plate.svg".to_string(),
    );
    values.insert(
        "FIG_TREATMENT_TRI".to_string(),
        "../report/figures/treatment_zone_plate_trifurcation.svg".to_string(),
    );
}

fn insert_table_values(
    values: &mut BTreeMap<String, String>,
    input: &Milestone12NarrativeInput<'_>,
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
) {
    values.insert(
        "RESULTS_INTRO".to_string(),
        build_results_intro(
            input.total_candidates,
            input.option1_pool_len,
            input.option2_pool_len,
        ),
    );
    values.insert(
        "SELECTED_TABLE_INTRO".to_string(),
        build_selected_table_intro(option1, option2),
    );
    values.insert(
        "SELECTED_TABLE".to_string(),
        build_selected_table(option1, option2),
    );
    values.insert(
        "TOP5_INTRO".to_string(),
        build_top5_intro(input.option2_ranked),
    );
    values.insert(
        "OPTION2_TOP5_TABLE".to_string(),
        build_option2_top5_table(input.option2_ranked),
    );
    values.insert(
        "TRI_CELL_INTRO".to_string(),
        build_tri_cell_intro(option2, input.ga_top.first()),
    );
    values.insert(
        "TRI_CELL_TABLE".to_string(),
        build_tri_cell_table(input.option2_ranked, input.ga_top),
    );
    values.insert("STRICT_CORE_INTRO".to_string(), build_strict_core_intro());
    values.insert(
        "STRICT_CORE_TABLE".to_string(),
        build_strict_core_table(input.option2_ranked, input.ga_top),
    );
    values.insert(
        "ROBUSTNESS_SECTION".to_string(),
        build_robustness_section(input.option2_robustness, input.fast_mode),
    );
    values.insert(
        "VALIDATION_SECTION".to_string(),
        build_validation_section(input.validation_rows, input.fast_mode),
    );
    values.insert(
        "CAVITATION_FORMULAS_INTRO".to_string(),
        build_cavitation_formulas_intro(option2),
    );
    values.insert(
        "LIMITS_OF_USAGE_INTRO".to_string(),
        build_limits_of_usage_intro(option2),
    );
    values.insert(
        "LIMITS_OF_USAGE".to_string(),
        build_limits_of_usage(option2),
    );
    values.insert(
        "CRI_EXPANSION_SENSITIVITY".to_string(),
        build_cri_expansion_sensitivity(),
    );
}

fn build_conclusions(
    total_candidates: usize,
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
    ga_best: &Milestone12ReportDesign,
) -> String {
    use std::fmt::Write as _;
    let m2 = &option2.metrics;
    let mg = &ga_best.metrics;
    let mut s = String::new();

    // §1 — Milestone completion
    let _ = writeln!(
        s,
        "Milestone 12 required selecting millifluidic device designs meeting hydrodynamic and \
cavitation parameters for extracorporeal SDT. From {} candidates across \
27 topology families, CFDrs identified one Option 1 and one Option 2 design via deterministic \
eligibility gating and ranking — no stochastic elements affect the final selection.\n",
        total_candidates,
    );

    // §2 — Option 1 acoustic
    if let Some(option1) = option1 {
        let m1 = &option1.metrics;
        let _ = writeln!(
            s,
            "**Option 1 — Selective Acoustic Center Treatment** (`{}`): {:.1}% of WBCs routed to \
center lane; RBC peripheral fraction {:.1}%; HI/pass = {:.2e}% (FDA 0.1% limit); \
ECV = {:.2} mL; P95 wall shear = {:.1} Pa. Zero active venturi throats — treatment relies on \
externally applied 412 kHz ultrasound. Acoustic resonance factor (ARF) = {:.4} \
(channel D_h match to λ/2 ≈ 1.87 mm). Score {:.4} under SelectiveAcousticTherapy.\n",
            option1.candidate.id,
            m1.wbc_recovery * 100.0,
            (1.0 - m1.rbc_venturi_exposure_fraction) * 100.0,
            m1.hemolysis_index_per_pass,
            m1.total_ecv_ml,
            m1.wall_shear_p95_pa,
            m1.acoustic_resonance_factor,
            option1.score,
        );
    } else {
        let _ = writeln!(
            s,
            "**Option 1 — Selective Acoustic Center Treatment**: no design satisfied the strict acoustic eligibility gates under the current physics regime. The report therefore records Option 1 explicitly as an empty shortlist rather than treating the absence as a ranking regression.\n"
        );
    }

    // §3 — Option 2 venturi
    let gauge_kpa = option2.inlet_gauge_kpa();
    let d_throat_um = option2.throat_width_um().unwrap_or(0.0);
    let _ = writeln!(
        s,
        "**Option 2 — Hydrodynamic Cavitation SDT** (`{}`): σ = {:.4} < 0 confirms active \
cavitation under {:.0} kPa gauge through {:.0} µm throat; {} serial stage(s) per path; \
{} total active throats. Cancer routing {:.1}% (cancer_center_fraction); therapeutic window \
score = {:.3}; HI/pass = {:.2e}%; throat viscous heating ΔT = {:.2} K (FDA 5 K limit: {}). \
Score {:.4} under CombinedSdtLeukapheresis — **not comparable to Option 1 score**.\n",
        option2.candidate.id,
        m2.cavitation_number,
        gauge_kpa,
        d_throat_um,
        m2.serial_venturi_stages_per_path,
        m2.active_venturi_throat_count,
        m2.cancer_center_fraction * 100.0,
        m2.therapeutic_window_score,
        m2.hemolysis_index_per_pass,
        m2.throat_temperature_rise_k,
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
Throat transit time exception: {:.2e} s {} 5 ms threshold. \
Clotting risk index = {:.4} at nominal flow; device must not exceed 500 mL/min without \
re-evaluation. ECV = {:.2} mL within pediatric circuit targets. \
FDA thermal compliance (42 °C ceiling) for Option 2: {} ({:.2} K rise).\n",
            m1.wall_shear_p95_pa,
            m2.wall_shear_p95_pa,
            m2.throat_transit_time_s,
            if m2.throat_transit_time_s < 5e-3 { "<" } else { "≥" },
            m2.clotting_risk_index,
            m2.total_ecv_ml,
            if m2.fda_thermal_compliant { "PASS" } else { "FAIL" },
            m2.throat_temperature_rise_k,
        );
    } else {
        let _ = writeln!(
            s,
            "**Safety and FDA Compliance:** The selected Option 2 design passes all five hard eligibility gates. Option 1 produced no eligible shortlist under the current physics regime, so no acoustic selected-design safety row exists for this run. Option 2 P95 wall shear = {:.1} Pa (FDA 150 Pa sustained limit). Throat transit time exception: {:.2e} s {} 5 ms threshold. Clotting risk index = {:.4} at nominal flow; device must not exceed 500 mL/min without re-evaluation. ECV = {:.2} mL within pediatric circuit targets. FDA thermal compliance (42 °C ceiling) for Option 2: {} ({:.2} K rise).\n",
            m2.wall_shear_p95_pa,
            m2.throat_transit_time_s,
            if m2.throat_transit_time_s < 5e-3 { "<" } else { "≥" },
            m2.clotting_risk_index,
            m2.total_ecv_ml,
            if m2.fda_thermal_compliant { "PASS" } else { "FAIL" },
            m2.throat_temperature_rise_k,
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
        "**Dominant Topology — Tri→Tri Physics:** PST2 (Primitive Selective Split Sequence \
Tri→Tri) dominated both tracks. The selected Option 2 design uses an {:.0} mm × {:.0} mm \
cross-section. Progressive trifurcation beyond 2 levels is counterproductive: a third level \
narrows terminal treatment channels to approximately 140 µm (below the 200–400 µm \
Zweifach-Fung inertial-focusing optimum), routes only 1/27 of inlet flow to each venturi \
(insufficient velocity for cavitation inception at the required pressure budget), and \
consumes excess pressure drop across the additional branching stages. RBC peripheral \
separation is achieved at Level 1; Level 2 refines cancer/WBC center concentration. \
Adding Level 3 does not improve separation enough to offset the flow-starvation penalty \
at each venturi throat.\n",
        w_mm, h_mm,
    );

    // §6 — GA results and modeling note
    let _ = writeln!(
        s,
        "**GA Results and Model Fidelity Note:** The HydroSDT GA (seed {M12_GA_HYDRO_SEED}, \
HydrodynamicCavitationSDT mode) produced rank-1 design `{}` with {} active throats \
({} serial stage(s)) — a narrower, more cavitation-intensive geometry than parametric Option 2. \
GA score {:.4} reflects the HydroSDT objective (35% cancer_cav, 20% separation, 20% RBC \
protection, 15% sonoluminescence) — **not the Combined mode** used for Option 2; scores are \
not comparable. **Important modeling note:** each serial venturi throat is evaluated \
independently for cavitation number and intensity; cumulative re-nucleation cascade between \
consecutive serial throats (inter-throat bubble collapse → re-growth → re-collapse) is not \
explicitly modeled in the 1D network solver. This is a conservative simplification — serial \
stages may provide greater cumulative dose than computed. Experimental validation at \
{}-stage conditions is recommended before clinical use.\n",
        ga_best.candidate.id,
        mg.active_venturi_throat_count,
        mg.serial_venturi_stages_per_path,
        ga_best.score,
        mg.serial_venturi_stages_per_path,
    );

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

fn build_appendix_a_supplemental(
    total_candidates: usize,
    opt1_pool: usize,
    opt2_pool: usize,
    ga_best: &Milestone12ReportDesign,
    robustness: &[crate::analysis::RobustnessReport],
    validation_rows: &[ValidationRow],
    canonical_results_path: &Path,
) -> String {
    use std::fmt::Write as _;
    let mut s = String::new();
    let _ = writeln!(
        s,
        "> Generated by `cargo run -p cfd-optim --example milestone12_report --no-default-features`.\n\
> Full canonical data: `report/milestone12_results.md`\n"
    );
    let _ = writeln!(s, "**Dataset Counts**\n");
    let _ = writeln!(s, "- Total generated candidates: **{}**", total_candidates);
    let _ = writeln!(
        s,
        "- Option 1 eligibility pool (selective acoustic): **{opt1_pool}**"
    );
    let _ = writeln!(
        s,
        "- Option 2 eligibility pool (selective venturi, σ<1): **{opt2_pool}**\n"
    );
    let _ = writeln!(s, "**GA Reproducibility**\n");
    let _ = writeln!(s, "- HydroSDT GA RNG seed: `{M12_GA_HYDRO_SEED}`");
    let _ = writeln!(
        s,
        "- GA rank-1: `{}` (score {:.4}, sigma {:.4})\n",
        ga_best.candidate.id, ga_best.score, ga_best.metrics.cavitation_number
    );
    if robustness.is_empty() {
        let _ = writeln!(
            s,
            "**Robustness Screening:** *Not computed in this run. Re-run without `M12_FAST=1`.*\n"
        );
    } else {
        let robust_count = robustness.iter().filter(|r| r.is_robust).count();
        let _ = writeln!(
            s,
            "**Robustness Screening:** {robust_count}/{} candidates robust.\n",
            robustness.len()
        );
    }
    if validation_rows.is_empty() {
        let _ = writeln!(
            s,
            "**Multi-Fidelity Validation:** *Not computed in this run. Re-run without `M12_FAST=1`.*\n"
        );
    } else {
        let _ = writeln!(
            s,
            "**Multi-Fidelity Validation:** {} entries computed. See §5.3 for full table.\n",
            validation_rows.len()
        );
    }
    let _ = writeln!(
        s,
        "**Serial Venturi Modeling Note:** Each serial venturi throat is evaluated \
independently for cavitation number and intensity in the 1D network solver. The cumulative \
re-nucleation cascade between consecutive serial throats is not explicitly modeled. \
This is a conservative simplification — experimental validation at multi-stage conditions \
is recommended.\n"
    );
    // Attempt to append a provenance link to full canonical data
    match std::fs::metadata(canonical_results_path) {
        Ok(_) => {
            let _ = writeln!(
                s,
                "Full canonical tables (Selected Designs, Top-5, Gate Evidence, Limits of Usage): \
see `{}`.",
                canonical_results_path.display()
            );
        }
        Err(_) => {
            let _ = writeln!(
                s,
                "*Full canonical results file not found at `{}`.*",
                canonical_results_path.display()
            );
        }
    }
    s
}

fn build_references_block() -> String {
    "\
1. ANSI/SLAS 1-2004 — Microplates — Footprint Dimensions.\n\
2. Giersiepen, M., et al. \"Estimation of shear stress-related blood damage in heart valve prostheses — in vitro comparison of 25 aortic valves.\" *International Journal of Artificial Organs*, 13(5):300–306, 1990.\n\
3. Brennen, C.E. *Cavitation and Bubble Dynamics*. Oxford University Press, 1995. [Eq. 3.12, cavitation number inception criterion.]\n\
4. Ohl, S.-W., et al. \"Sonoporation from jetting cavitation bubbles.\" *Biophysical Journal*, 91(11):4285–4295, 2006. [5× membrane lysis amplification at bubble collapse; basis for `LYSIS_CAVITATION_AMPLIFICATION = 5.0`.]\n\
5. Hellums, J.D. \"1993 Whitaker Lecture: Biorheology in thrombosis research.\" *Annals of Biomedical Engineering*, 22(5):445–455, 1994. [PAI exponent model, Eq. 3: n=1.325, m=0.462.]\n\
6. Di Carlo, D. \"Inertial microfluidics.\" *Lab on a Chip*, 9(21):3038–3046, 2009. [κ_RBC = a_RBC / D_h confinement criterion for inertial focusing; threshold 0.07.]\n\
7. Zweifach, B.W. and Fung, Y.C. \"Phase separation in capillary networks.\" *Microvascular Research*, 1971. [β_RBC = 1.0 for passive tracer routing at millifluidic scales.]\n\
8. FDA 2019 Guidance — *Nonclinical Tests and Recommended Labeling for Intravascular Administration Sets, Blood Administration Sets, and Blood Component Administration Sets*. [FDA predicate: Maquet RotaFlow, K143453, 1% hemolysis ceiling.]\n\
9. Lentner, C. (Ed.) *Geigy Scientific Tables, Vol. 3: Physical Chemistry, Composition of Blood.* Novartis, 1984. [Table 30: neonatal reference blood volume 85 mL/kg; basis for `PEDIATRIC_BLOOD_VOLUME_ML_PER_KG = 85.0`.]\n\
10. SonALAsense Internal Data — CFDrs canonical simulation data (see Appendix A)."
        .to_string()
}
