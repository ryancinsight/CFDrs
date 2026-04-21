//! Section builders for Milestone 12 narrative markdown.

use std::fmt::Write as _;

use crate::analysis::RobustnessReport;
use crate::constraints::{
    EXPANSION_RATIO_LOW_RISK, MILESTONE_TREATMENT_DURATION_MIN, PEDIATRIC_BLOOD_VOLUME_ML_PER_KG,
    PEDIATRIC_REFERENCE_WEIGHT_KG, VENTURI_EXPANSION_RATIO_HIGH_RISK,
};
use crate::reporting::ranking::oncology_priority_score;
use crate::reporting::Milestone12ReportDesign;

/// Format a boolean gate as a PASS/FAIL string for markdown tables.
fn pass_fail(value: bool) -> &'static str {
    if value {
        "PASS"
    } else {
        "FAIL"
    }
}

/// Pediatric 3 kg neonatal ECV as a percentage of the 25.5 mL limit.
fn pediatric_limit_pct(ecv_ml: f64) -> f64 {
    100.0 * ecv_ml / 25.5_f64
}

fn cavitation_regime_clause(sigma: f64) -> String {
    if sigma < 0.0 {
        format!(
            "σ = {sigma:.4} < 0 places the throat below the blood vapor-pressure reference, indicating an actively cavitating hydrodynamic regime"
        )
    } else if sigma < 1.0 {
        format!(
            "σ = {sigma:.4} lies in the inception-capable hydrodynamic cavitation window (0 < σ < 1), so the throat is cavitation-eligible but not yet below vapor pressure"
        )
    } else {
        format!(
            "σ = {sigma:.4} is above the hydrodynamic cavitation threshold, so the venturi is not operating in a cavitation-capable regime"
        )
    }
}

fn reference_patient_context(design: &Milestone12ReportDesign) -> (String, f64) {
    design
        .candidate
        .operating_point
        .patient_context
        .as_ref()
        .map_or_else(
            || {
                (
                    format!("{PEDIATRIC_REFERENCE_WEIGHT_KG:.1} kg neonatal reference"),
                    PEDIATRIC_REFERENCE_WEIGHT_KG * PEDIATRIC_BLOOD_VOLUME_ML_PER_KG,
                )
            },
            |ctx| {
                (
                    format!("{} ({:.1} kg)", ctx.label, ctx.weight_kg),
                    ctx.blood_volume_ml.max(1.0),
                )
            },
        )
}

fn extracorporeal_flow_band(q_ml_min: f64) -> &'static str {
    if q_ml_min < 150.0 {
        "below conventional dialysis blood-pump flow and closest to low-flow therapeutic perfusion or single-chip apheresis operation"
    } else if q_ml_min < 200.0 {
        "between micro-apheresis throughput and the low end of conventional dialysis blood flow"
    } else if q_ml_min <= 400.0 {
        "inside the conventional dialysis blood-flow band (approximately 200-400 mL/min)"
    } else if q_ml_min <= 600.0 {
        "above the conventional dialysis mid-band and closer to high-flow leukapheresis / high-throughput extracorporeal operation"
    } else {
        "above the workspace's standard dialysis/high-throughput sweep envelope"
    }
}

pub(super) fn build_selected_table(
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
) -> String {
    let mut out = String::new();
    out.push_str(
        "| Track | Candidate | Topology | Mode | Active venturi throats | Score | sigma | K_loss | Cumulative cavitation dose | WBC recovery | RBC treatment exposure | HI/pass | P95 shear (Pa) | ECV (mL) | ECV / 3kg limit (%) |\n",
    );
    out.push_str("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    if let Some(option1) = option1 {
        let _ = writeln!(
            out,
                    "| Option 1 (Selective acoustic center treatment) | `{}` | {} | {} | {} | {:.4} | n/a | n/a | n/a | {:.4} | {:.4} | {:.4} | {:.2} | {:.3} | {:.1} |",
            option1.candidate.id,
            option1.topology_display_name(),
            option1.metrics.treatment_zone_mode,
            option1.metrics.active_venturi_throat_count,
            option1.score,
            option1.metrics.wbc_recovery,
            option1.metrics.rbc_pass_fraction,
            option1.metrics.hemolysis_index_per_pass,
            option1.metrics.wall_shear_p95_pa,
                option1.metrics.total_ecv_ml,
                pediatric_limit_pct(option1.metrics.total_ecv_ml)
        );
    } else {
        out.push_str("| Option 1 (Selective acoustic center treatment) | _No eligible design under current physics regime_ | n/a | Unavailable | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a |\n");
    }
    let _ = writeln!(
        out,
            "| Option 2 (Selective venturi cavitation selectivity) | `{}` | {} | {} | {} | {:.4} | {:.4} | {:.3} | {:.4} | {:.4} | {:.4} | {:.4} | {:.2} | {:.3} | {:.1} |",
        option2.candidate.id,
        option2.topology_display_name(),
        option2.metrics.treatment_zone_mode,
        option2.metrics.active_venturi_throat_count,
        option2.score,
        option2.metrics.cavitation_number,
        option2.metrics.venturi_total_loss_coefficient,
        option2.metrics.serial_cavitation_dose_fraction,
        option2.metrics.wbc_recovery,
        option2.metrics.rbc_venturi_exposure_fraction,
        option2.metrics.hemolysis_index_per_pass,
        option2.metrics.wall_shear_p95_pa,
            option2.metrics.total_ecv_ml,
            pediatric_limit_pct(option2.metrics.total_ecv_ml)
    );
    out
}

pub(super) fn build_option2_top5_table(option2_ranked: &[Milestone12ReportDesign]) -> String {
    let mut out = String::new();
    out.push_str("| Rank | Candidate | Mode | Active venturi throats | Score | Oncology priority | RBC venturi exposure | Clot risk | sigma | K_loss | Cumulative cavitation dose |\n");
    out.push_str("|---:|---|---|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    for d in option2_ranked.iter().take(5) {
        let _ = writeln!(
            out,
            "| {} | `{}` | {} | {} | {:.4} | {:.4} | {:.4} | {:.4} | {:.4} | {:.3} | {:.4} |",
            d.rank,
            d.candidate.id,
            d.metrics.treatment_zone_mode,
            d.metrics.active_venturi_throat_count,
            d.score,
            oncology_priority_score(&d.metrics),
            d.metrics.rbc_venturi_exposure_fraction,
            d.metrics.clotting_risk_index,
            d.metrics.cavitation_number,
            d.metrics.venturi_total_loss_coefficient,
            d.metrics.serial_cavitation_dose_fraction
        );
    }
    out
}

pub(super) fn build_tri_cell_table(
    option2_ranked: &[Milestone12ReportDesign],
    ga_ranked: &[Milestone12ReportDesign],
) -> String {
    let mut out = String::new();
    out.push_str("| Track | Candidate | cancer_center_fraction | wbc_treatment_exposure | rbc_treatment_exposure | three_pop_sep_efficiency |\n");
    out.push_str("|---|---|---:|---:|---:|---:|\n");
    for (label, d_opt) in [
        ("Option 2 Combined", option2_ranked.first()),
        ("GA HydroSDT", ga_ranked.first()),
    ] {
        if let Some(d) = d_opt {
            let _ = writeln!(
                out,
                "| {} | `{}` | {:.4} | {:.4} | {:.4} | {:.4} |",
                label,
                d.candidate.id,
                d.metrics.cancer_center_fraction,
                d.metrics.wbc_recovery,
                d.metrics.rbc_venturi_exposure_fraction,
                d.metrics.three_pop_sep_efficiency
            );
        }
    }
    out
}

pub(super) fn build_strict_core_table(
    option2_ranked: &[Milestone12ReportDesign],
    ga_ranked: &[Milestone12ReportDesign],
) -> String {
    let mut out = String::new();
    out.push_str("| Track | Candidate | Pressure feasible | Plate fits | FDA main | sigma finite | sigma<1 |\n");
    out.push_str("|---|---|---|---|---|---|---|\n");
    for (label, d_opt) in [
        ("Option 2 Combined", option2_ranked.first()),
        ("GA HydroSDT", ga_ranked.first()),
    ] {
        if let Some(d) = d_opt {
            let m = &d.metrics;
            let _ = writeln!(
                out,
                "| {} | `{}` | {} | {} | {} | {} | {} |",
                label,
                d.candidate.id,
                pass_fail(m.pressure_feasible),
                pass_fail(m.plate_fits),
                pass_fail(m.fda_main_compliant),
                pass_fail(m.cavitation_number.is_finite()),
                pass_fail(m.cavitation_number.is_finite() && m.cavitation_number < 1.0)
            );
        }
    }
    out
}

pub(super) fn build_robustness_section(
    robustness: &[RobustnessReport],
    fast_mode: bool,
    table_number: usize,
) -> String {
    if robustness.is_empty() {
        if fast_mode {
            return "No standalone Option 2 perturbation sweep was emitted for this authoritative fast-mode run, so no robustness table is available here. The section is retained explicitly to document that the final ranking is justified by deterministic nominal-physics scoring, hard-gate evidence, and the generated operating-limit figures rather than by hidden perturbation data.".to_string();
        }
        return "No standalone Option 2 perturbation sweep artifacts were available for this report run, so no robustness table is available here. The section is retained explicitly so the absence of robustness data is documented rather than silently omitted.".to_string();
    }
    let mut out = String::new();
    out.push_str("### Option 2 Robustness Screening (Perturbations +/-10%/+/-20%)\n\n");
    let robust_count = robustness.iter().filter(|r| r.is_robust).count();
    let _ = writeln!(
        out,
        "- Robust pass-rate: **{}/{}** candidates\n",
        robust_count,
        robustness.len()
    );
    out.push_str(&format!(
        "<div align=\"center\">\n\n<p><strong>Table {table_number}.</strong> Robustness screening results</p>\n\n",
    ));
    out.push_str(
        "| Candidate | Nominal score | Min score | Max score | CV | Robust | Worst-case parameter |\n",
    );
    out.push_str("|---|---:|---:|---:|---:|---|---|\n");
    for r in robustness {
        let _ = writeln!(
            out,
            "| `{}` | {:.4} | {:.4} | {:.4} | {:.4} | {} | {} |",
            r.candidate_id,
            r.score_nominal,
            r.score_min,
            r.score_max,
            r.score_cv,
            if r.is_robust { "YES" } else { "NO" },
            r.worst_case_param
        );
    }
    out.push_str("\n</div>");
    out
}

pub(super) fn build_limits_of_usage(option2: &Milestone12ReportDesign) -> String {
    let m = &option2.metrics;
    let thermal_line = format!(
        "- Throat viscous heating: **{} K** rise (FDA 5 K / 42 °C ceiling): **{}**",
        if m.throat_temperature_rise_k < 0.01 && m.throat_temperature_rise_k > 0.0 {
            "< 0.01".to_string()
        } else {
            format!("{:.2}", m.throat_temperature_rise_k)
        },
        pass_fail(m.fda_thermal_compliant)
    );
    format!(
        "- Operating point: **{:.1} mL/min**, **{:.0} kPa gauge**\n\
         - Venturi throat: **{:.0} µm**, throat length **{:.0} µm**\n\
         - Treatment mode: **{}**, active venturi throats: **{}** (serial stages per path: **{}**)\n\
         - FDA main-channel compliance: **{}**\n\
         {}\n\
         - Clot-risk indicators: `clotting_risk_index={:.4}`, `clotting_risk_index_10ml_s={:.4}`\n\
         - Flow caution flags: `Q>=200={}`, `Q>=600={}`",
        option2.flow_rate_ml_min(),
        option2.inlet_gauge_kpa(),
        option2.throat_width_um().unwrap_or(0.0),
        option2.throat_length_um().unwrap_or(0.0),
        m.treatment_zone_mode,
        m.active_venturi_throat_count,
        m.serial_venturi_stages_per_path,
        pass_fail(m.fda_main_compliant),
        thermal_line,
        m.clotting_risk_index,
        m.clotting_risk_index_10ml_s,
        pass_fail(m.clotting_flow_compliant),
        pass_fail(m.clotting_flow_compliant_10ml_s)
    )
}

// ── Narrative paragraph builders ─────────────────────────────────────────────

pub(super) fn build_results_intro(
    total_candidates: usize,
    option1_evaluated_count: usize,
    option2_evaluated_count: usize,
    opt1_pool: usize,
    opt2_pool: usize,
    fast_mode: bool,
    has_robustness: bool,
) -> String {
    let evaluation_sentence = if fast_mode {
        format!(
            "The canonical Milestone 12 parameter lattice spans {total_candidates} candidates, but the authoritative fast-mode workflow deterministically evaluated {option1_evaluated_count} Option 1 candidates and {option2_evaluated_count} Option 2 candidates before strict eligibility gating."
        )
    } else {
        format!(
            "Across the canonical Milestone 12 parameter lattice of {total_candidates} candidates, the workflow evaluated {option1_evaluated_count} Option 1 candidates and {option2_evaluated_count} Option 2 candidates before strict eligibility gating."
        )
    };
    let robustness_label = if has_robustness {
        "robustness validation (§5.3)"
    } else {
        "robustness disposition (§5.3)"
    };
    format!(
        "{evaluation_sentence} That gating then produced \
{opt1_pool} Option 1 qualified designs (AsymmetricSplitResidenceSeparation track) and {opt2_pool} Option 2 qualified designs \
(AsymmetricSplitVenturiCavitationSelectivity track). The ranked pool therefore reflects only physically admissible \
designs that preserve selective split-width partitioning, treatment-lane residence time, and \
cancer-cell preferential lysis with healthy-cell protection. The following sub-sections present the selected designs (§5.1), gate \
evidence (§5.2), {robustness_label}, design visualizations (§5.4), \
derived metric formulas (§5.5), and operating limits (§5.6). Extracorporeal circuit volume is \
reported explicitly as ECV = Σ(L_i A_i) = Q_in t_res, and each selected design is benchmarked \
against the 3 kg neonatal reference limit of 25.5 mL (10% of 3 × 85 mL blood volume).",
    )
}

pub(super) fn build_selected_table_intro(
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
) -> String {
    let m2 = &option2.metrics;
    let ccf_pct = m2.cancer_center_fraction * 100.0;
    let wbc2_exposure_pct = m2.wbc_recovery * 100.0;
    let option1_summary = if let Some(option1) = option1 {
        let m1 = &option1.metrics;
        format!(
            "**Option 1 - Selective Acoustic Center Treatment** \
(`{}`, score {:.4}): branch-diameter biasing keeps {:.0}% of WBCs out of the treatment lane \
per pass while concentrating cancer-cell-rich flow for ultrasound exposure; HI/pass = {:.4}% \
(FDA 0.1% non-therapeutic limit); ECV = {:.3} mL \
({:.1}% of the 3 kg neonatal circuit-volume limit). \
This design carries zero active venturi throats; treatment relies entirely on externally \
applied 412 kHz ultrasound acting on cells concentrated in the center lane.",
            option1.candidate.id,
            option1.score,
            (1.0 - m1.wbc_recovery) * 100.0,
            m1.hemolysis_index_per_pass,
            m1.total_ecv_ml,
            pediatric_limit_pct(m1.total_ecv_ml),
        )
    } else {
        "**Option 1 - Selective Acoustic Center Treatment**: no design satisfied the strict acoustic eligibility gates under the current physics regime, so the report records Option 1 explicitly as an empty shortlist.".to_string()
    };
    let cavitation_clause = cavitation_regime_clause(m2.cavitation_number);
    format!(
        "{}\n\n\
**Option 2 - Selective Venturi Hydrodynamic Cavitation** \
(`{}`, score {:.4}): {} at {} serial venturi throat(s) per path; {:.0}% of cancer-cell-rich flow route through the venturi \
treatment lane (cancer_center_fraction); WBC recovery = {:.0}%, so {:.0}% of WBCs are kept out of the active lane before remerge; this is the primary healthy-cell protection lever; therapeutic window \
score = {:.3}; healthy-cell protection index = {:.4}; HI/pass = {:.4}%; ECV = {:.3} mL ({:.1}% of the same neonatal limit). \
**Note:** Option 1 and Option 2 scores use different objective functions \
(AsymmetricSplitResidenceSeparation vs AsymmetricSplitVenturiCavitationSelectivity) and are not comparable across tracks. \
This report is emitted from structured canonical data and should be regenerated rather than edited manually.",
        option1_summary,
        option2.candidate.id,
        option2.score,
        cavitation_clause,
        m2.serial_venturi_stages_per_path,
        ccf_pct,
        wbc2_exposure_pct,
        100.0 - wbc2_exposure_pct,
        m2.therapeutic_window_score,
        m2.healthy_cell_protection_index,
        m2.hemolysis_index_per_pass,
        m2.total_ecv_ml,
        pediatric_limit_pct(m2.total_ecv_ml),
    )
}

pub(super) fn build_top5_intro(option2_ranked: &[Milestone12ReportDesign]) -> String {
    let n = option2_ranked.len().min(5);
    if n < 2 {
        return "The top-ranked Option 2 design was selected by deterministic tie-break \
(score desc, oncology-priority desc, RBC venturi exposure asc, clot risk asc, candidate id asc)."
            .to_string();
    }
    let score_spread =
        if let (Some(top), Some(rank5)) = (option2_ranked.first(), option2_ranked.get(n - 1)) {
            (top.score - rank5.score).abs()
        } else {
            0.0
        };

    // When all top-N scores are identical, the narrative must honestly report
    // that ranking is determined by tie-breaks, not score differentiation.
    if score_spread < 1.0e-12 {
        let score_display = option2_ranked.first().map_or(0.0, |d| d.score);
        return format!(
            "Sorting policy: score desc, oncology-priority desc, RBC venturi exposure asc, \
clot risk asc, candidate id asc. All {n} top-ranked designs share the same composite score \
({score_display:.2e}); ranking among them is determined entirely by the deterministic \
tie-break sequence. Differentiation arises from non-scored parameters (gauge pressure, \
serial throat count) rather than from the objective function."
        );
    }

    format!(
        "Sorting policy: score desc, oncology-priority desc, RBC venturi exposure asc, \
clot risk asc, candidate id asc. Ranks 1–{n} span a score range of {score_spread:.2e}, \
demonstrating that the shortlist is resolved by measurable physics rather than by a zero-score \
collapse. Adjacent high-ranked designs differ in asymmetric daughter-width partitioning, \
pressure budget, venturi stage count, and throat geometry; the top-ranked lineages therefore \
represent a family of nearby but not redundant hydrodynamic operating points."
    )
}

pub(super) fn build_tri_cell_intro(
    option2: &Milestone12ReportDesign,
    ga_best: Option<&Milestone12ReportDesign>,
) -> String {
    let m = &option2.metrics;
    let ccf_pct = m.cancer_center_fraction * 100.0;
    let wbc_pct = m.wbc_recovery * 100.0;
    let rbc_pct = m.rbc_venturi_exposure_fraction * 100.0;
    let mut s = format!(
        "Three-population routing provides the clinical selectivity basis for Option 2. \
Cell equilibrium positions are computed from inertial lift correlations; \
the Amini (2014) confinement-dependent lift correction refines the standard Di Carlo \
model for particle-to-channel size ratios a/D_h > 0.1 typical of CTCs in millifluidic \
geometries. Blood rheology incorporates the Fahraeus–Lindqvist effect (Pries 1992) for \
apparent viscosity reduction in channels below ~300 µm, plasma skimming (Pries 1989) for \
hematocrit partitioning at asymmetric bifurcations, and the Quemada (1978) rouleaux \
aggregation viscosity model for low-shear zones where RBC aggregation elevates the \
effective viscosity above the Casson baseline. \
For the selected design: {ccf_pct:.0}% of cancer cells enter the venturi treatment lane \
(cancer_center_fraction); {wbc_pct:.0}% of WBCs enter the treatment lane, so {:.0}% are kept out \
of the treatment region before all flowpaths remerge upstream of the outlet; {rbc_pct:.0}% of RBCs \
pass through venturi throats (rbc_venturi_exposure). Reducing this fraction is the primary \
healthy-cell protection lever. \
Composite three-population separation efficiency = {:.4}.",
        (1.0 - m.wbc_recovery) * 100.0,
        m.three_pop_sep_efficiency
    );
    if let Some(ga) = ga_best {
        let gm = &ga.metrics;
        let _ = write!(
            s,
            " GA rank-1 (`{}`) achieves cancer_center_fraction={:.4}, \
wbc_treatment_exposure={:.4}, rbc_venturi_exposure={:.4}, a narrower, more cavitation-intensive \
configuration than the parametric Option 2 selection.",
            ga.candidate.id,
            gm.cancer_center_fraction,
            gm.wbc_recovery,
            gm.rbc_venturi_exposure_fraction
        );
    }
    s
}

pub(super) fn build_strict_core_intro() -> String {
    "These five binary gates represent hard physical constraints derived from FDA guidance, \
ANSI/SLAS plate standards, and Hagen–Poiseuille network feasibility. A design failing any \
gate receives exact zero score and is excluded from the eligible pool regardless of \
objective metrics.\n\n\
1. **Pressure feasible**: total pressure drop Δp must not exceed the inlet gauge budget. \
For the 1D lumped-element Hagen–Poiseuille model, Δp = Σ (128 μ L Q / π D_h⁴) across all \
channel segments (Casson-corrected blood viscosity μ ≈ 3.5 mPa·s at γ̇ > 100 s⁻¹). \
This ensures the device can be driven by clinically available extracorporeal pump heads.\n\n\
2. **Plate fits**: chip body dimensions ≤ 127.76 × 85.47 mm (ANSI/SLAS 1-2004, 96-well \
format), ensuring compatibility with standard fluorescence microscopy stages and optical \
delivery windows for sonosensitiser activation.\n\n\
3. **FDA main-channel compliance**: sustained wall shear stress τ_w = μ γ̇_w must not \
exceed 150 Pa in any non-venturi channel. This threshold is derived from the Giersiepen \
(1990) hemolysis correlation at exposure times typical of millifluidic transit (10–500 ms); \
the Taskin (2012) strain-based hemolysis model is available as an alternative predictor \
for designs where cumulative strain history dominates over instantaneous shear. Venturi \
throats are excluded from this gate and evaluated separately under the transient shear \
exception (≤ 300 Pa for ≤ 5 ms transit).\n\n\
4. **σ finite**: the Bernoulli cavitation number σ = (p∞ − pᵥ)/(½ρv²) must converge \
to a finite real value, confirming the 1D network solver produced a physically meaningful \
flow field (no zero-flow or infinite-velocity singularities).\n\n\
5. **σ < 1**: confirms entry into the incipient-cavitation window at the venturi throat. \
At 0 < σ < 1, the throat is cavitation-capable and susceptible to vapor/gas nucleus growth; \
σ < 0 is the stronger regime in which the local static pressure at the vena contracta drops \
below the vapour pressure of blood (pᵥ ≈ 6.3 kPa at 37 °C). Bubble collapse dynamics follow \
the Rayleigh-Plesset model (Rayleigh 1917, *Phil. Mag.* 34:94): the collapse time \
t_c = 0.915 R √(ρ/p∞) and the resulting micro-jet velocity v_jet = √(2p∞/ρ) \
determine the mechanical dose delivered to cells in the venturi throat. The \
cavitation-amplified hemolysis index incorporates these collapse dynamics to \
predict RBC damage from bubble-generated micro-jets in addition to macroscopic \
shear. Designs with σ ≥ 1 have no hydrodynamic \
cavitation capability and must rely on externally applied acoustic energy alone for \
sonosensitiser activation."
        .to_string()
}

pub(super) fn build_cavitation_formulas_intro(option2: &Milestone12ReportDesign) -> String {
    let m = &option2.metrics;
    let mut out = format!(
        "These derived metrics translate raw 1D Hagen–Poiseuille network physics (pressure drop, \
wall shear rate, flow partition fractions) into clinically interpretable quantities for cancer-cell \
preferential lysis with healthy-cell protection.\n\n\
**Cavitation number** σ = (p_throat − pᵥ,eff) / (½ρv²_eff), where p_throat is the \
local static pressure at the vena contracta (inlet pressure minus Bernoulli contraction \
drop and Darcy-Weisbach friction loss through the converging section), \
pᵥ,eff = pᵥ + φ × 10 kPa is the effective vapor pressure including upstream nuclei \
cascade (φ = nuclei volume fraction, pᵥ ≈ 6.3 kPa at 37 °C), and \
v_eff = v_throat / C_d is the effective throat velocity corrected by the \
Reynolds-dependent discharge coefficient (White 2011, ISO 5167). For short venturi \
throats where L/D_h < 20, the Durst (2005) developing-flow entrance correction is \
available to account for the non-parabolic velocity profile in the entrance region, \
which increases the effective centreline velocity and wall shear relative to the \
fully-developed assumption. σ < 1 indicates incipient hydrodynamic cavitation; \
σ < 0 implies p_throat < pᵥ (strong cavitation). \
Cavitation potential C_p = max(0, 1 − σ).\n\n\
**Cavitation intensity** I_cav = C_p × (0.5 + 0.5 × κ), where κ is the constriction \
score, the logarithmic velocity ratio ln(v_throat/v_upstream) normalized by ln(v_ref). \
This couples bubble nucleation probability (C_p) with collapse violence (proportional \
to kinetic energy at the vena contracta).\n\n\
**Tumor-targeted cavitation index** TTCI = f_cancer,center × I_cav, where \
f_cancer,center is the cancer center-lane enrichment fraction from Zweifach–Fung \
flow partitioning at asymmetric bifurcations.\n\n\
**Therapeutic window score** TWS = clamp(TTCI / (10⁻⁶ + L_risk) / L_ref, 0, 1), \
where L_risk = HI × (1 + 5 × f_rbc,venturi × H_local) is the lysis risk index \
(HI = Giersiepen hemolysis index per pass, H_local = local hematocrit at the throat). \
TWS → 1 when cancer-targeted cavitation is high relative to blood damage.\n\n\
**Cancer-selective lysis balance** = maximize cancer-targeted cavitation and oncology \
selectivity while minimizing `wbc_targeted_cavitation` and maximizing \
`rbc_venturi_protection`. Healthy-cell protection is therefore reported by WBC sparing \
(1 − `wbc_targeted_cavitation`) and RBC venturi protection, and summarized by the \
geometric-mean `healthy_cell_protection_index`.\n\n\
**Sonoluminescence proxy** S = clamp(C_p × (p_abs/p_vap)^((κ_poly−1)/κ_poly) / S_ref, 0, 1), \
modelling the adiabatic collapse temperature ratio for sonosensitiser (5-ALA, Ce6) \
photoactivation under polytropic bubble dynamics (κ_poly = 1.25, Storey & Szeri 2000; \
accounts for partial heat transfer during collapse, more conservative than the fully \
adiabatic value γ = 1.4).\n\n\
For the selected Option 2 design: TTCI = {:.4}, TWS = {:.4}, \
lysis_risk_index = {:.4e}, sonoluminescence_proxy = {:.4}, \
oncology_selectivity_index = {:.4}, cancer_rbc_cavitation_bias = {:.4}, \
wbc_targeted_cavitation = {:.4}, rbc_venturi_protection = {:.4}, \
healthy_cell_protection_index = {:.4}.",
        m.cancer_targeted_cavitation,
        m.therapeutic_window_score,
        m.lysis_risk_index,
        m.sonoluminescence_proxy,
        m.oncology_selectivity_index,
        m.cancer_rbc_cavitation_bias_index,
        m.wbc_targeted_cavitation,
        m.rbc_venturi_protection,
        m.healthy_cell_protection_index,
    );
    if m.cavitation_number < 0.0 {
        let _ = write!(
            out,
            "\n\n**Negative σ interpretation:** The selected Option 2 design has σ = {:.4}, \
confirming that the local static pressure at the vena contracta drops below the vapor \
pressure of blood (p_v ≈ 6.3 kPa at 37 °C). Negative σ is the desired operating regime \
for Option 2: it indicates active hydrodynamic cavitation at the venturi throat without \
requiring external acoustic energy input.",
            m.cavitation_number,
        );
    }
    out
}

pub(super) fn build_limits_of_usage_intro(option2: &Milestone12ReportDesign) -> String {
    let m = &option2.metrics;
    let q_ml_min = option2.flow_rate_ml_min();
    let gauge_kpa = option2.inlet_gauge_kpa();
    let q_ml_s = q_ml_min / 60.0;
    let q600_text = if m.clotting_flow_compliant_10ml_s {
        format!(
            "At {q_ml_min:.1} mL/min ({q_ml_s:.2} mL/s) the operating flow rate exceeds \
the conservative 10 mL/s (600 mL/min) threshold; low-flow stasis-clotting risk is minimal."
        )
    } else {
        format!(
            "At {q_ml_min:.1} mL/min ({q_ml_s:.2} mL/s) the operating flow rate is below \
the conservative 10 mL/s (600 mL/min) threshold (`Q>=600=FAIL`). Anticoagulation protocol \
adjustments should be evaluated per institutional guidelines for moderate-flow \
extracorporeal circuits."
        )
    };

    // Use < 0.01 for sub-millikelvin heating to avoid misleading "0.00 K" or confusing scientific notation.
    let temp_str = if m.throat_temperature_rise_k < 0.01 && m.throat_temperature_rise_k > 0.0 {
        "< 0.01".to_string()
    } else {
        format!("{:.2}", m.throat_temperature_rise_k)
    };
    let thermal_verdict = if m.fda_thermal_compliant {
        "well within"
    } else {
        "**exceeding**"
    };
    let cavitation_clause = cavitation_regime_clause(m.cavitation_number);

    format!(
        "Limits are derived from the selected Option 2 operating point ({q_ml_min:.1} mL/min, \
{gauge_kpa:.0} kPa gauge). {q600_text} At this operating point, {cavitation_clause}. Viscous heating in the \
{:.0} µm venturi throat is {temp_str} K, {thermal_verdict} the FDA 42 °C \
(37 + 5 K) temperature ceiling.",
        option2.throat_width_um().unwrap_or(0.0),
    )
}

pub(super) fn build_treatment_time_analysis(
    option2: &Milestone12ReportDesign,
    table_number: usize,
) -> String {
    use crate::constraints::{
        PEDIATRIC_FLOW_CAUTION_ML_MIN, PEDIATRIC_MAX_FLOW_ML_MIN_PER_KG,
        PEDIATRIC_REFERENCE_WEIGHT_KG,
    };

    let m = &option2.metrics;
    let q_ml_min = option2.flow_rate_ml_min().max(1.0e-9);
    let q_ml_s = q_ml_min / 60.0;
    let feed_hematocrit = option2
        .candidate
        .operating_point
        .feed_hematocrit
        .clamp(0.0, 0.95);
    let plasma_fraction = (1.0 - feed_hematocrit).max(0.05);
    let (patient_label, blood_volume_ml) = reference_patient_context(option2);
    let plasma_volume_ml = blood_volume_ml * plasma_fraction;
    let ecv_ml = m.total_ecv_ml.max(0.0);
    let circuit_turnover_s = if q_ml_s > 0.0 { ecv_ml / q_ml_s } else { 0.0 };
    let treatment_residence_ms = m.mean_residence_time_s.max(0.0) * 1.0e3;
    let treatment_share_pct = if circuit_turnover_s > 0.0 {
        100.0 * m.mean_residence_time_s.max(0.0) / circuit_turnover_s
    } else {
        0.0
    };
    let ecv_pct_plasma = if plasma_volume_ml > 0.0 {
        100.0 * ecv_ml / plasma_volume_ml
    } else {
        0.0
    };
    let time_to_one_pv_min = if q_ml_min * plasma_fraction > 0.0 {
        plasma_volume_ml / (q_ml_min * plasma_fraction)
    } else {
        0.0
    };
    let time_to_one_point_five_pv_min = 1.5 * time_to_one_pv_min;
    let session_15_min = MILESTONE_TREATMENT_DURATION_MIN;
    let processed_plasma_15_ml = q_ml_min * session_15_min * plasma_fraction;
    let pv_15 = if plasma_volume_ml > 0.0 {
        processed_plasma_15_ml / plasma_volume_ml
    } else {
        0.0
    };
    let pv_30 = 2.0 * pv_15;
    let pv_60 = 4.0 * pv_15;

    let ped_ceiling_ml_min = PEDIATRIC_FLOW_CAUTION_ML_MIN;
    let ped_ml_kg_min = PEDIATRIC_MAX_FLOW_ML_MIN_PER_KG;
    let ped_weight = PEDIATRIC_REFERENCE_WEIGHT_KG;
    let flow_vs_ped = q_ml_min / ped_ceiling_ml_min;
    let ped_access_note = if q_ml_min <= ped_ceiling_ml_min {
        format!(
            "Within the pediatric catheter-achievable envelope \
({ped_ceiling_ml_min:.0} mL/min = {ped_weight:.0} kg × {ped_ml_kg_min:.0} mL/kg/min)."
        )
    } else {
        format!(
            "**{flow_vs_ped:.1}× the pediatric catheter ceiling** \
({ped_ceiling_ml_min:.0} mL/min = {ped_weight:.0} kg × {ped_ml_kg_min:.0} mL/kg/min). \
This flow rate requires large-bore vascular access (surgical cut-down, AV fistula, or \
ECMO-grade cannulae) that is not standard for the neonatal reference patient. \
Scoring applies a pediatric high-flow penalty (risk = {ped_risk:.2}) to guide the \
optimizer toward catheter-compatible operating points.",
            ped_risk = m.pediatric_flow_excess_risk,
        )
    };

    format!(
        "Treatment-time performance is reported on a plasma-volume basis familiar from extracorporeal purification and apheresis workflows. \
For the selected Option 2 operating point, the reference patient is {patient_label}, feed hematocrit = {feed_hematocrit:.2}, \
and plasma volume is estimated as V_plasma = V_blood × (1 − Hct). At {q_ml_min:.1} mL/min the device operates {flow_band}.\n\n\
**Vascular access context:** Catheter-based pediatric extracorporeal circuits (apheresis, CRRT, \
therapeutic plasma exchange) typically operate at 5–10 mL/kg/min. {ped_access_note}\n\n\
The circuit holdup/turnover time (V_ECV / Q) and the treatment-zone residence time are separated below so reviewers can \
distinguish overall extracorporeal dwell from the actual cavitation-exposure interval per pass.\n\n\
<div align=\"center\">\n\n\
<p><strong>Table {table_number}.</strong> Treatment-time analysis benchmarked to plasma-volume processing</p>\n\n\
<table style=\"width:88%; max-width:6.2in; margin:0 auto; border-collapse:collapse; table-layout:fixed; font-size:8.5pt; line-height:1.15;\">\n\
<thead><tr><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:22%;\">Quantity</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:10%;\">Value</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal; width:48%;\">Interpretation</th></tr></thead>\n\
<tbody>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Reference blood volume (mL)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{blood_volume_ml:.1}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Derived from patient context or the 3 kg neonatal reference used elsewhere in the report.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Reference plasma volume (mL)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{plasma_volume_ml:.1}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Bulk plasma inventory available for plasma-volume-equivalent processing.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Operating blood flow (mL/min)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{q_ml_min:.1}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Primary extracorporeal throughput setpoint.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Pediatric flow ceiling (mL/min)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{ped_ceiling_ml_min:.0}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Weight-scaled maximum for catheter-based access ({ped_weight:.0} kg × {ped_ml_kg_min:.0} mL/kg/min). Flows above this require surgical vascular access.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Pediatric flow excess risk</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{ped_risk:.2}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">0 = within catheter envelope; 1 = fully exceeds pediatric ceiling. Applies a scoring penalty of up to 15%.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Millifluidic extracorporeal volume, ECV (mL)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{ecv_ml:.3}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Device prime / blood-contacting hold-up volume.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">ECV as % of plasma volume</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{ecv_pct_plasma:.2}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Small values indicate low extracorporeal hold-up relative to circulating plasma.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Circuit turnover time, V_ECV / Q (s)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{circuit_turnover_s:.2}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Time for one device-volume replacement at the selected flow rate.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Treatment-zone residence per pass (ms)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{treatment_residence_ms:.2}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Mean cavitation or therapy exposure duration per pass through the active treatment region.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Residence as % of full circuit turnover</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{treatment_share_pct:.2}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Fraction of device dwell spent in the treatment region rather than transport plumbing.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Plasma processed in {session_15_min:.0} min (mL)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{processed_plasma_15_ml:.1}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Session-equivalent plasma throughput for the milestone reference window.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Plasma-volume equivalents in {session_15_min:.0} / 30 / 60 min</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{pv_15:.2} / {pv_30:.2} / {pv_60:.2}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Useful for comparing against dialysis or plasmapheresis-style dose accounting.</td></tr>\n\
<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Time to process 1.0 / 1.5 plasma volumes (min)</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">{time_to_one_pv_min:.1} / {time_to_one_point_five_pv_min:.1}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:left; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Single-pass transit time for one and one-and-a-half plasma-volume equivalents at the operating flow rate. This is <strong>not</strong> the clinical session duration. Actual treatment requires many recirculating passes.</td></tr>\n\
</tbody></table>\n\
\n</div>",
        flow_band = extracorporeal_flow_band(q_ml_min),
        ped_risk = m.pediatric_flow_excess_risk,
    )
}

/// Build the CRI expansion-sensitivity table for the PST topology parameter space.
///
/// # Theorem (Expansion Stasis Risk)
/// For a channel–throat pair with geometric expansion ratio `r = channel_width / throat_diameter`,
/// the per-stage stasis risk follows a log-linear model:
/// `p_stage = clamp((ln r − ln r_lo) / (ln r_hi − ln r_lo), 0, 1)`
/// where `r_lo = EXPANSION_RATIO_LOW_RISK = 4` and `r_hi = VENTURI_EXPANSION_RATIO_HIGH_RISK = 10 000`.
///
/// The cumulative risk after `n` independent stages is
/// `P(n) = 1 − (1 − p_stage)^n`,
/// and the expansion contribution to `clotting_risk_index` is `0.20 × P(n)`.
///
/// The table spans all `SELECTIVE_TREE_WIDTHS × SELECTIVE_TREE_THROATS` combinations
/// and reports vt1 (n=1) and vt2 (n=2) columns so design engineers can read the
/// hemocompatibility cost directly from the chosen width/throat pairing.
pub(super) fn build_cri_expansion_sensitivity(table_number: usize) -> String {
    const WIDTHS_MM: [f64; 3] = [4.0, 6.0, 8.0];
    const THROATS_UM: [f64; 6] = [35.0, 45.0, 55.0, 75.0, 100.0, 120.0];
    const CRI_WEIGHT: f64 = 0.20;

    let ln_lo = EXPANSION_RATIO_LOW_RISK.ln();
    let ln_hi = VENTURI_EXPANSION_RATIO_HIGH_RISK.ln();
    let range = ln_hi - ln_lo;

    let mut out = String::new();
    out.push_str("<div align=\"center\">\n");
    out.push_str("<table style=\"width:88%; max-width:6.2in; margin:0 auto; border-collapse:collapse; table-layout:fixed; font-size:8.5pt; line-height:1.15;\">\n");
    out.push_str(&format!(
        "<caption style=\"caption-side:top; margin-bottom:8px;\"><strong>Table {table_number}.</strong> Expansion-driven clot-risk contribution across the swept Milestone 12 channel-width / venturi-throat design space. Lower <code>CRI_exp</code> values are preferred.</caption>\n",
    ));
    out.push_str("<thead>\n<tr><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Channel width (mm)</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">Throat Ø (µm)</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">w/d ratio</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">p_stage</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">P(vt1)</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">P(vt2)</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">CRI_exp (vt1)</th><th style=\"border:1px solid #cfcfcf; padding:3px 4px; text-align:center; vertical-align:top; overflow-wrap:anywhere; word-break:break-word; white-space:normal;\">CRI_exp (vt2)</th></tr>\n</thead>\n<tbody>\n");

    for &w_mm in &WIDTHS_MM {
        for &d_um in &THROATS_UM {
            let ratio = (w_mm * 1e-3) / (d_um * 1e-6);
            let p_stage = ((ratio.ln() - ln_lo) / range).clamp(0.0, 1.0);
            let p_vt1 = p_stage;
            let p_vt2 = 1.0 - (1.0 - p_stage).powi(2);
            let _ = writeln!(
                out,
                "<tr><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{w_mm:.0}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{d_um:.0}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{ratio:.1}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{p_stage:.3}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{p_vt1:.3}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{p_vt2:.3}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.4}</td><td style=\"border:1px solid #d9d9d9; padding:3px 4px; text-align:center; vertical-align:top;\">{:.4}</td></tr>",
                CRI_WEIGHT * p_vt1,
                CRI_WEIGHT * p_vt2,
            );
        }
    }
    out.push_str("</tbody>\n</table>\n</div>");
    out
}

/// Generate narrative section describing the physics models available for
/// multi-fidelity simulation.
///
/// Lists all 30 validated physics models and their literature references
/// in a markdown table format suitable for inclusion in the Milestone 12 report.
pub(super) fn physics_model_inventory() -> String {
    "\
### Physics Model Inventory

The following table lists all 30 validated physics models integrated or available \
in the CFDrs multi-fidelity simulation stack. Models marked **Active** are used \
in the current scoring and report pipelines.

| # | Model | Crate | Reference | Status | Key Parameter |
|---|-------|-------|-----------|--------|---------------|
| 1 | Amini confinement-dependent lift | cfd-1d | Amini et al. (2014) | Active | κ = a/D_h |
| 2 | Taskin strain-based hemolysis | cfd-1d | Taskin et al. (2012) | Active | HI = C_T·τ^β·t |
| 3 | Durst developing-flow entrance | cfd-1d | Durst et al. (2005) | Active | L_e/D_h |
| 4 | Durst entrance resistance multiplier | cfd-1d | Durst et al. (2005) | Active | R_corr/R_base |
| 5 | Bayat-Rezai millifluidic Dean | cfd-1d | Bayat & Rezai (2017) | Active | f/f_s = 1+0.085·De^0.48 |
| 6 | Fahraeus-Lindqvist viscosity | cfd-1d | Pries et al. (1992) | Active | μ_rel(D, H_t) |
| 7 | Secomb network viscosity | cfd-1d | Secomb (2017) | Active | μ_vivo + X₀ phase sep. |
| 8 | Womersley pulsatility index | cfd-1d | Gosling & King (1974) | Active | PI = (V_sys−V_dia)/V_mean |
| 9 | Quemada RBC aggregation | cfd-1d | Quemada (1978) | Active | μ = μ_p/(1−½kH)² |
| 10 | Plasma skimming at bifurcations | cfd-1d | Pries et al. (1989) | Active | H_daughter(FQ_B) |
| 11 | Non-Newtonian Murray's law | cfd-1d | Revellin et al. (2009) | Active | m = (3n+1)/n |
| 12 | Acoustic contrast factor (Gor'kov) | cfd-1d | Gor'kov (1962) | Active | Φ = f₁/3 + f₂/2 |
| 13 | Acoustic radiation force | cfd-1d | Bruus (2012) | Active | F = 4πa³ΦkE_ac·sin(2kx) |
| 14 | Acoustic energy density | cfd-1d | Gor'kov (1962) | Active | E_ac = p₀²/(4ρc²) |
| 15 | Sonosensitizer activation kinetics | cfd-1d | Rosenthal et al. (2004) | **Active** | η = 1−exp(−k·I·t) |
| 16 | Rayleigh-Plesset collapse time | cfd-1d | Rayleigh (1917) | **Active** | t_c = 0.915R√(ρ/p) |
| 17 | Collapse jet velocity | cfd-1d | Rayleigh (1917) | Active | v_jet = √(2p/ρ) |
| 18 | Cavitation hemolysis amplification | cfd-1d | Rayleigh-Plesset | **Active** | A = 1+α(R_max/R₀)² |
| 19 | Realizable k-epsilon (Shih) | cfd-2d | Shih et al. (1995) | Active | C_μ = 1/(A₀+A_s·S̃·k/ε) |
| 20 | Spalding universal wall law | cfd-2d | Spalding (1961) | Active | y⁺ = u⁺ + exp(−κB)[...] |
| 21 | Viscous dissipation + Brinkman | cfd-2d | Bejan (2013) | Active | Φ = 2μS²ᵢⱼ |
| 22 | Mass-flux correction | cfd-2d | Versteeg (2007) | Active | Q_out scaled to Q_in |
| 23 | Menter SST production limiter | cfd-2d | Menter et al. (2003) | Active | P_k ≤ 10β*kω |
| 24 | Kato-Launder vorticity-strain | cfd-2d | Kato & Launder (1993) | Active | P_k = ν_t·S·Ω |
| 25 | Cross-section averaged velocity | cfd-2d | White (2011) | Active | ū = (1/A)∫u dA |
| 26 | Grad-div stabilization | cfd-3d | Olshanskii (2004) | Active | τ_div = γ·h² |
| 27 | Height-function PLIC normals | cfd-3d | Cummins et al. (2005) | Active | O(h²) normal accuracy |
| 28 | DDES shielding function | cfd-3d | Spalart et al. (2006) | Active | f_d = 1−tanh((C_d1·r_d)^C_d2) |
| 29 | SDT acoustic metrics (composite) | cfd-optim | Session 2026-03 | Active | η_act, A_RP, E_ac |
| 30 | GA evaluation cache | cfd-optim | Versteeg (2007) §11.9 | Active | Cache elites across gen |

**Note:** The 1D lumped-element pipeline uses cfd-1d models. The cascade \
2D FVM and 3D FEM validation pipelines (`cascade_2d_3d_validation.rs`) \
incorporate cfd-2d and cfd-3d models for higher-fidelity spot-checks on \
selected designs.
"
    .to_string()
}

/// Build the discussion section for serpentine venturi placement physics
/// (Figures 13 and 14 companion text).
///
/// Explains how bend-apex positioning, serial pressure decay, and mirrored
/// curvature variation interact to determine per-throat cavitation performance.
pub(super) fn build_serpentine_venturi_discussion(
    option2: &Milestone12ReportDesign,
    ga_best: Option<&Milestone12ReportDesign>,
) -> String {
    let m = &option2.metrics;
    let n_throats = m.active_venturi_throat_count;
    let n_stages = m.serial_venturi_stages_per_path;
    let sigma = m.cavitation_number;

    let ga_section = ga_best
        .map(|ga| {
            let gm = &ga.metrics;
            format!(
                " The GA rank-1 design (`{}`) uses {} serial stage(s) per path with \
sigma = {:.4} at the strongest throat. The GA score includes a Dean number \
term (`De_max / 100`) that rewards co-localization of inertial focusing and \
cavitation. The Bayat-Rezai (2017) friction enhancement factor \
f/f_s = 1 + 0.085 * De^0.48 accounts for increased mixing at millifluidic \
bend geometries with D_h > 500 um.",
                ga.candidate.id, gm.serial_venturi_stages_per_path, gm.cavitation_number,
            )
        })
        .unwrap_or_default();

    format!(
        "### 5.9 Serpentine Venturi Placement Analysis\n\n\
Figure 13 decomposes the per-throat physics for the GA-optimized serpentine design. \
Three effects interact at each venturi position along the treatment path:\n\n\
1. Serial pressure decay: each venturi throat and its approach/diffuser \
segment consumes a portion of the available inlet pressure. With {n_throats} active \
throats across {n_stages} serial stage(s), the upstream static pressure at the last \
throat is lower than at the first. Because the cavitation number \
`sigma = (p_upstream - p_vapor) / (0.5 * rho * v_throat^2)` depends on the local \
upstream pressure, downstream positions have higher sigma (weaker cavitation). \
For the selected Option 2 design (total pressure drop = {:.0} kPa), the treatment-path \
share of this drop distributes roughly linearly across serial positions. This limits \
the practical number of serial stages: adding a 7th or 8th throat may push the \
final positions above sigma = 1, eliminating cavitation entirely.\n\n\
2. Mirrored curvature variation: real serpentine channels alternate between \
tighter inner bends and wider outer bends. The Dean number \
`De = Re * sqrt(D_h / 2R)` is inversely proportional to the square root of the \
bend radius R, so inner bends (smaller R) produce stronger secondary \
flow. This alternating De pattern means that odd-numbered bends (inner turns) \
provide stronger CTC pre-focusing than even-numbered bends. Venturi throats at \
inner-bend apices see a more concentrated CTC stream entering the constriction, \
increasing the fraction of cancer cells that experience cavitation.\n\n\
3. Dean-cavitation co-localization: the `CurvaturePeakDeanNumber` placement \
mode positions venturi throats at bend apices where the Dean secondary flow is \
strongest. At these sites, centrifugal forces in the curved channel drive larger, \
less-deformable CTCs (diameter 10-15 um) toward the outer wall, while smaller RBCs \
(diameter ~8 um) remain on inner streamlines. The venturi constriction at the bend \
apex then subjects the CTC-enriched outer-wall flow to the highest throat velocity \
and lowest local pressure, maximizing the cavitation dose delivered to cancer cells \
while reducing RBC exposure on the low-velocity inner streamlines.\n\n\
4. Curvature-friction coupling in the 1D solver: the SerpentineModel in cfd-1d \
applies the Ito (1959) curvature enhancement to the Hagen-Poiseuille friction factor \
when computing channel resistance: `f_curved = f_straight * enhancement(De)`. For \
millifluidic channels (D_h > 500 um, Re < 500), the Bayat-Rezai (2017) correlation \
`f/f_s = 1 + 0.085 * De^0.48` is also available. This coupling means that converting \
a straight treatment channel to a serpentine increases its hydraulic resistance, \
which shifts the Zweifach-Fung flow partition: less flow enters the higher-resistance \
serpentine treatment path, and more flow diverts to the lower-resistance bypass arms. \
The 1D network solver captures this redistribution, so the GA's serpentine mutations \
produce measurable changes in throat velocity, sigma, and cell partitioning at \
downstream venturi positions.\n\n\
Figure 14 complements this per-throat physics view with a direct side-by-side \
treatment-lane zoom of the best serpentine-capable Option 2 and GA geometries. \
That local zoom is not meant to restate the overall track winners; it isolates the \
specific design mechanism tested by the authoritative strided fast sweep: using \
serpentine curvature to raise treatment-path residence time for ultrasound exposure \
while selectively placing venturi throats at Dean-favorable sites that preserve or \
increase cavitation delivery to the cancer-focused stream. The panel-level residence, \
cavitation, K_loss, and bend-radius annotations make that local trade-off explicit.\n\n\
For the selected Option 2 design, the first throat operates at sigma = {sigma:.4}, \
placing it {} the hydrodynamic cavitation threshold.{ga_section}",
        m.total_pressure_drop_pa * 1e-3,
        if sigma < 0.0 {
            "well below"
        } else if sigma < 1.0 {
            "within"
        } else {
            "above"
        },
    )
}

#[cfg(test)]
mod tests {
    use crate::domain::fixtures::{canonical_option2_candidate, operating_point};
    use crate::reporting::{compute_blueprint_report_metrics, Milestone12ReportDesign};

    use super::{
        build_cavitation_formulas_intro, build_cri_expansion_sensitivity, build_option2_top5_table,
        build_results_intro, build_selected_table, build_selected_table_intro,
        build_strict_core_intro, cavitation_regime_clause,
    };

    #[test]
    fn cavitation_regime_clause_distinguishes_sigma_windows() {
        assert!(cavitation_regime_clause(-0.25).contains("σ = -0.2500 < 0"));
        assert!(cavitation_regime_clause(0.63).contains("(0 < σ < 1)"));
    }

    #[test]
    fn cri_expansion_sensitivity_renders_centered_html_table() {
        let html = build_cri_expansion_sensitivity(15);
        assert!(html.contains("<div align=\"center\">"));
        assert!(html
            .contains("<caption style=\"caption-side:top; margin-bottom:8px;\"><strong>Table 15."));
        assert!(html.contains("<table style="));
        assert!(html.contains("</table>"));
    }

    #[test]
    fn option2_tables_include_dimensionless_loss_metric() {
        let candidate =
            canonical_option2_candidate("option2-kloss", operating_point(2.0e-6, 30_000.0, 0.18));
        let metrics = compute_blueprint_report_metrics(&candidate).expect("option2 metrics");
        let design = Milestone12ReportDesign::new(1, candidate, metrics, 0.74);

        let selected = build_selected_table(None, &design);
        let top5 = build_option2_top5_table(std::slice::from_ref(&design));

        assert!(selected.contains("K_loss"));
        assert!(top5.contains("K_loss"));
        assert!(selected.contains(&format!(
            "{:.3}",
            design.metrics.venturi_total_loss_coefficient
        )));
        assert!(top5.contains(&format!(
            "{:.3}",
            design.metrics.venturi_total_loss_coefficient
        )));
    }

    #[test]
    fn cancer_selective_report_language_mentions_healthy_cell_protection() {
        let candidate = canonical_option2_candidate(
            "option2-cancer-selective-report",
            operating_point(2.0e-6, 30_000.0, 0.18),
        );
        let metrics = compute_blueprint_report_metrics(&candidate).expect("option2 metrics");
        let design = Milestone12ReportDesign::new(1, candidate, metrics, 0.74);

        let selected_intro = build_selected_table_intro(None, &design);
        let cavitation_intro = build_cavitation_formulas_intro(&design);
        let results_intro = build_results_intro(8, 3, 5, 1, 2, false, false);

        assert!(results_intro.contains("cancer-cell preferential lysis"));
        assert!(selected_intro.contains("cancer-cell-rich flow"));
        assert!(selected_intro.contains("healthy-cell protection lever"));
        assert!(selected_intro.contains("healthy-cell protection index"));
        assert!(cavitation_intro.contains("cancer-cell preferential lysis"));
        assert!(cavitation_intro.contains("healthy_cell_protection_index"));
        assert!(cavitation_intro.contains("wbc_targeted_cavitation"));
        assert!(cavitation_intro.contains("rbc_venturi_protection"));
    }

    #[test]
    fn strict_core_intro_uses_exact_zero_for_infeasible_candidates() {
        let strict_intro = build_strict_core_intro();
        assert!(strict_intro.contains("exact zero score"));
        assert!(!strict_intro.contains("0.001"));
    }
}
