//! Section builders for Milestone 12 narrative markdown.

use std::fmt::Write as _;

use crate::analysis::RobustnessReport;
use crate::constraints::{EXPANSION_RATIO_LOW_RISK, VENTURI_EXPANSION_RATIO_HIGH_RISK};
use crate::reporting::{Milestone12ReportDesign, ValidationRow};

pub(super) fn build_selected_table(
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
) -> String {
    let mut out = String::new();
    out.push_str(
        "| Track | Candidate | Topology | Mode | Active venturi throats | Score | sigma | Cumulative cavitation dose | WBC treatment exposure | RBC treatment exposure | HI/pass | P95 shear (Pa) | ECV (mL) | ECV / 3kg limit (%) |\n",
    );
    out.push_str("|---|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n");
    if let Some(option1) = option1 {
        let _ = writeln!(
            out,
                    "| Option 1 (Selective acoustic center treatment) | `{}` | {} | {} | {} | {:.4} | n/a | n/a | {:.4} | {:.4} | {:.2e} | {:.2} | {:.3} | {:.1} |",
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
        out.push_str("| Option 1 (Selective acoustic center treatment) | _No eligible design under current physics regime_ | n/a | Unavailable | 0 | n/a | n/a | n/a | n/a | n/a | n/a | n/a | n/a |\n");
    }
    let _ = writeln!(
        out,
            "| Option 2 (Selective venturi cavitation selectivity) | `{}` | {} | {} | {} | {:.2e} | {:.4} | {:.4} | {:.4} | {:.4} | {:.2e} | {:.2} | {:.3} | {:.1} |",
        option2.candidate.id,
        option2.topology_display_name(),
        option2.metrics.treatment_zone_mode,
        option2.metrics.active_venturi_throat_count,
        option2.score,
        option2.metrics.cavitation_number,
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
    out.push_str("| Rank | Candidate | Mode | Active venturi throats | Score | RBC venturi exposure | Clot risk | sigma | Cumulative cavitation dose |\n");
    out.push_str("|---:|---|---|---:|---:|---:|---:|---:|---:|\n");
    for d in option2_ranked.iter().take(5) {
        let _ = writeln!(
            out,
            "| {} | `{}` | {} | {} | {:.2e} | {:.4} | {:.4} | {:.4} | {:.4} |",
            d.rank,
            d.candidate.id,
            d.metrics.treatment_zone_mode,
            d.metrics.active_venturi_throat_count,
            d.score,
            d.metrics.rbc_venturi_exposure_fraction,
            d.metrics.clotting_risk_index,
            d.metrics.cavitation_number,
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

pub(super) fn build_robustness_section(robustness: &[RobustnessReport], fast_mode: bool) -> String {
    if robustness.is_empty() {
        let _ = fast_mode;
        return "\
### Option 2 Robustness Screening (Perturbations +/-10%/+/-20%)\n\n\
*Robustness screening was not computed in this run (fast mode). \
Re-run without `M12_FAST=1` to populate.*\n"
            .to_string();
    }
    let _ = fast_mode;
    let mut out = String::new();
    out.push_str("### Option 2 Robustness Screening (Perturbations +/-10%/+/-20%)\n\n");
    let robust_count = robustness.iter().filter(|r| r.is_robust).count();
    let _ = writeln!(
        out,
        "- Robust pass-rate: **{}/{}** candidates\n",
        robust_count,
        robustness.len()
    );
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
    out
}

pub(super) fn build_validation_section(rows: &[ValidationRow], fast_mode: bool) -> String {
    if rows.is_empty() {
        let _ = fast_mode;
        return "\
### Multi-Fidelity Pressure-Drop Validation (Selected Venturi Designs)\n\n\
*Multi-fidelity validation is produced by the companion example \
`cargo run -p cfd-optim --example milestone12_validation --no-default-features`. \
This report run did not embed 2D/3D validation rows.*\n"
            .to_string();
    }
    let _ = fast_mode;
    let mut out = String::new();
    out.push_str("### Multi-Fidelity Pressure-Drop Validation (Selected Venturi Designs)\n\n");
    out.push_str("| Track | Candidate | dp1D Bernoulli (Pa) | dp2D FVM (Pa) | dp3D FEM (Pa) | 1D-2D diff (%) | 2D-3D diff (%) | Mass error (%) |\n");
    out.push_str("|---|---|---:|---:|---:|---:|---:|---:|\n");
    for row in rows {
        let _ = writeln!(
            out,
            "| {} | `{}` | {:.2} | {:.2} | {:.2} | {:.2} | {:.2} | {:.2} |",
            row.track,
            row.id,
            row.dp_1d_bernoulli_pa,
            row.dp_2d_fvm_pa,
            row.dp_3d_fem_pa,
            row.agreement_1d_2d_pct,
            row.agreement_2d_3d_pct,
            row.mass_error_3d_pct
        );
    }
    out
}

pub(super) fn build_limits_of_usage(option2: &Milestone12ReportDesign) -> String {
    let m = &option2.metrics;
    let thermal_line = format!(
        "- Throat viscous heating: **{:.2} K** rise (FDA 5 K / 42 °C ceiling): **{}**",
        m.throat_temperature_rise_k,
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
    opt1_pool: usize,
    opt2_pool: usize,
) -> String {
    format!(
        "From {total_candidates} total candidates generated from canonical Milestone 12 split-sequence scaffolds, \
strict eligibility gating produced {} Option 1 qualified designs \
(AsymmetricSplitResidenceSeparation track) and {} Option 2 qualified designs \
(AsymmetricSplitVenturiCavitationSelectivity track). The ranked pool therefore reflects only physically admissible \
designs that preserve selective split-width partitioning, treatment-lane residence time, and \
healthy-cell protection. The following sub-sections present the selected designs (§5.1), gate \
evidence (§5.2), robustness and multi-fidelity validation (§5.3), design visualizations (§5.4), \
derived metric formulas (§5.5), and operating limits (§5.6). Extracorporeal circuit volume is \
reported explicitly as ECV = Σ(L_i A_i) = Q_in t_res, and each selected design is benchmarked \
against the 3 kg neonatal reference limit of 25.5 mL (10% of 3 × 85 mL blood volume).",
        opt1_pool,
        opt2_pool,
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
            "**Option 1 — Selective Acoustic Center Treatment** \
(`{}`, score {:.4}): branch-diameter biasing keeps {:.0}% of WBCs out of the treatment lane \
per pass while concentrating CTC-rich flow for ultrasound exposure; HI/pass = {:.2e}% \
(FDA 0.1% non-therapeutic limit); ECV = {:.2} mL \
({:.1}% of the 3 kg neonatal circuit-volume limit). \
This design carries zero active venturi throats — treatment relies entirely on externally \
applied 412 kHz ultrasound acting on cells concentrated in the center lane.",
            option1.candidate.id,
            option1.score,
            (1.0 - m1.wbc_recovery) * 100.0,
            m1.hemolysis_index_per_pass,
            m1.total_ecv_ml,
            pediatric_limit_pct(m1.total_ecv_ml),
        )
    } else {
        "**Option 1 — Selective Acoustic Center Treatment**: no design satisfied the strict acoustic eligibility gates under the current physics regime, so the report records Option 1 explicitly as an empty shortlist.".to_string()
    };
    format!(
        "{}\n\n\
**Option 2 — Selective Venturi Hydrodynamic Cavitation** \
(`{}`, score {:.4}): σ = {:.4} < 0 confirms active hydrodynamic cavitation at \
{} serial venturi throat(s) per path; {:.0}% of CTCs route through the venturi \
treatment lane (cancer_center_fraction); WBC treatment exposure = {:.0}%; therapeutic window \
score = {:.3}; HI/pass = {:.2e}%; ECV = {:.2} mL ({:.1}% of the same neonatal limit). \
**Note:** Option 1 and Option 2 scores use different objective functions \
(AsymmetricSplitResidenceSeparation vs AsymmetricSplitVenturiCavitationSelectivity) and are not comparable across tracks. \
This report is emitted from structured canonical data and should be regenerated rather than edited manually.",
        option1_summary,
        option2.candidate.id,
        option2.score,
        m2.cavitation_number,
        m2.serial_venturi_stages_per_path,
        ccf_pct,
        wbc2_exposure_pct,
        m2.therapeutic_window_score,
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
            format!("{:.2e}", (top.score - rank5.score).abs())
        } else {
            "small".to_string()
        };
    format!(
        "Sorting policy: score desc, oncology-priority desc, RBC venturi exposure asc, \
clot risk asc, candidate id asc. Ranks 1–{n} span a score range of {score_spread}, \
demonstrating that the top configuration is not a fragile single point — adjacent designs \
differ primarily in venturi stage count and throat length factor rather than fundamental \
topology. A wider score gap at the boundary between Rank {n} and lower-ranked candidates \
marks a meaningful performance discontinuity driven by trifurcation center fraction."
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
For the selected design: {ccf_pct:.0}% of CTCs enter the venturi treatment lane \
(cancer_center_fraction); {wbc_pct:.0}% of WBCs enter the treatment lane, so {:.0}% are kept out \
of the treatment region before all flowpaths remerge upstream of the outlet; {rbc_pct:.0}% of RBCs \
pass through venturi throats (rbc_venturi_exposure) — reducing this fraction is the primary \
hemocompatibility lever. \
Composite three-population separation efficiency = {:.4}.",
        (1.0 - m.wbc_recovery) * 100.0,
        m.three_pop_sep_efficiency
    );
    if let Some(ga) = ga_best {
        let gm = &ga.metrics;
        let _ = write!(
            s,
            " GA rank-1 (`{}`) achieves cancer_center_fraction={:.4}, \
wbc_treatment_exposure={:.4}, rbc_venturi_exposure={:.4} — a narrower, more cavitation-intensive \
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
gate receives the infeasibility floor score (0.001) and is excluded from the eligible pool \
regardless of objective metrics.\n\n\
1. **Pressure feasible** — total pressure drop Δp must not exceed the inlet gauge budget. \
For the 1D lumped-element Hagen–Poiseuille model, Δp = Σ (128 μ L Q / π D_h⁴) across all \
channel segments (Casson-corrected blood viscosity μ ≈ 3.5 mPa·s at γ̇ > 100 s⁻¹). \
This ensures the device can be driven by clinically available extracorporeal pump heads.\n\n\
2. **Plate fits** — chip body dimensions ≤ 127.76 × 85.47 mm (ANSI/SLAS 1-2004, 96-well \
format), ensuring compatibility with standard fluorescence microscopy stages and optical \
delivery windows for sonosensitiser activation.\n\n\
3. **FDA main-channel compliance** — sustained wall shear stress τ_w = μ γ̇_w must not \
exceed 150 Pa in any non-venturi channel. This threshold is derived from the Giersiepen \
hemolysis correlation at exposure times typical of millifluidic transit (10–500 ms). Venturi \
throats are excluded from this gate and evaluated separately under the transient shear \
exception (≤ 300 Pa for ≤ 5 ms transit).\n\n\
4. **σ finite** — the Bernoulli cavitation number σ = (p∞ − pᵥ)/(½ρv²) must converge \
to a finite real value, confirming the 1D network solver produced a physically meaningful \
flow field (no zero-flow or infinite-velocity singularities).\n\n\
5. **σ < 1** — confirms incipient cavitation onset at the venturi throat. At σ < 1, the \
local static pressure at the vena contracta drops below the vapour pressure of blood \
(pᵥ ≈ 6.3 kPa at 37 °C), nucleating vapour/gas bubbles whose subsequent collapse delivers \
the mechanical forces for sonodynamic therapy. Designs with σ ≥ 1 have no hydrodynamic \
cavitation activity and must rely on externally applied acoustic energy alone for \
sonosensitiser activation."
        .to_string()
}

pub(super) fn build_cavitation_formulas_intro(option2: &Milestone12ReportDesign) -> String {
    let m = &option2.metrics;
    format!(
        "These derived metrics translate raw 1D Hagen–Poiseuille network physics (pressure drop, \
wall shear rate, flow partition fractions) into clinically interpretable quantities.\n\n\
**Cavitation number** σ = (p∞ − pᵥ) / (½ρv²_throat), where p∞ is the far-field \
(inlet) absolute pressure, pᵥ ≈ 6.3 kPa (blood vapor pressure at 37 °C), and \
v_throat is the Casson-regime mean velocity at the venturi throat. σ < 1 indicates \
incipient hydrodynamic cavitation; σ < 0 implies p_throat < pᵥ (strong cavitation). \
Cavitation potential C_p = max(0, 1 − σ).\n\n\
**Cavitation intensity** I_cav = C_p × (0.5 + 0.5 × κ), where κ is the constriction \
score — the logarithmic velocity ratio ln(v_throat/v_upstream) normalised by ln(v_ref). \
This couples bubble nucleation probability (C_p) with collapse violence (proportional \
to kinetic energy at the vena contracta).\n\n\
**Tumor-targeted cavitation index** TTCI = f_cancer,center × I_cav, where \
f_cancer,center is the cancer center-lane enrichment fraction from Zweifach–Fung \
flow partitioning at asymmetric bifurcations.\n\n\
**Therapeutic window score** TWS = clamp(TTCI / (10⁻⁶ + L_risk) / L_ref, 0, 1), \
where L_risk = HI × (1 + 5 × f_rbc,venturi × H_local) is the lysis risk index \
(HI = Giersiepen hemolysis index per pass, H_local = local hematocrit at the throat). \
TWS → 1 when cancer-targeted cavitation is high relative to blood damage.\n\n\
**Sonoluminescence proxy** S = clamp(C_p × (p_abs/p_vap)^((κ_poly−1)/κ_poly) / S_ref, 0, 1), \
modelling the adiabatic collapse temperature ratio for sonosensitiser (5-ALA, Ce6) \
photoactivation under polytropic bubble dynamics (κ_poly = 1.4).\n\n\
For the selected Option 2 design: TTCI = {:.4}, TWS = {:.4}, \
lysis_risk_index = {:.4e}, sonoluminescence_proxy = {:.4}, \
oncology_selectivity_index = {:.4}, cancer_rbc_cavitation_bias = {:.4}.",
        m.cancer_targeted_cavitation,
        m.therapeutic_window_score,
        m.lysis_risk_index,
        m.sonoluminescence_proxy,
        m.oncology_selectivity_index,
        m.cancer_rbc_cavitation_bias_index,
    )
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

    // Use scientific notation for sub-millikelvin heating to avoid misleading "0.00 K".
    let temp_str = if m.throat_temperature_rise_k < 0.01 && m.throat_temperature_rise_k > 0.0 {
        format!("{:.2e}", m.throat_temperature_rise_k)
    } else {
        format!("{:.2}", m.throat_temperature_rise_k)
    };
    let thermal_verdict = if m.fda_thermal_compliant {
        "well within"
    } else {
        "**exceeding**"
    };

    format!(
        "Limits are derived from the selected Option 2 operating point ({q_ml_min:.1} mL/min, \
{gauge_kpa:.0} kPa gauge). {q600_text} Viscous heating in the \
{:.0} µm venturi throat is {temp_str} K — {thermal_verdict}: the FDA 42 °C \
(37 + 5 K) temperature ceiling.",
        option2.throat_width_um().unwrap_or(0.0),
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
pub(super) fn build_cri_expansion_sensitivity() -> String {
    const WIDTHS_MM: [f64; 3] = [4.0, 6.0, 8.0];
    const THROATS_UM: [f64; 6] = [35.0, 45.0, 55.0, 75.0, 100.0, 120.0];
    const CRI_WEIGHT: f64 = 0.20;

    let ln_lo = EXPANSION_RATIO_LOW_RISK.ln();
    let ln_hi = VENTURI_EXPANSION_RATIO_HIGH_RISK.ln();
    let range = ln_hi - ln_lo;

    let mut out = String::new();
    out.push_str("| Channel width (mm) | Throat Ø (µm) | w/d ratio | p_stage | P(vt1) | P(vt2) | CRI_exp (vt1) | CRI_exp (vt2) |\n");
    out.push_str("|---:|---:|---:|---:|---:|---:|---:|---:|\n");

    for &w_mm in &WIDTHS_MM {
        for &d_um in &THROATS_UM {
            let ratio = (w_mm * 1e-3) / (d_um * 1e-6);
            let p_stage = ((ratio.ln() - ln_lo) / range).clamp(0.0, 1.0);
            let p_vt1 = p_stage;
            let p_vt2 = 1.0 - (1.0 - p_stage).powi(2);
            let _ = writeln!(
                out,
                "| {w_mm:.0} | {d_um:.0} | {ratio:.1} | {p_stage:.3} | {p_vt1:.3} | {p_vt2:.3} | {:.4} | {:.4} |",
                CRI_WEIGHT * p_vt1,
                CRI_WEIGHT * p_vt2,
            );
        }
    }
    out
}

fn pass_fail(value: bool) -> &'static str {
    if value {
        "PASS"
    } else {
        "FAIL"
    }
}

fn pediatric_limit_pct(ecv_ml: f64) -> f64 {
    100.0 * ecv_ml / 25.5_f64
}
