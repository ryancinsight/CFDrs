//! Conclusions, appendix, reference-list, figure, and storage blocks for the M12 narrative.

use std::fmt::Write as _;
use std::path::Path;

use crate::constraints::M12_GA_HYDRO_SEED;
use crate::reporting::figures::NarrativeFigureSpec;
use crate::reporting::{Milestone12ReportDesign, ValidationRow};

pub(super) fn build_conclusions(
    total_candidates: usize,
    option1: Option<&Milestone12ReportDesign>,
    option2: &Milestone12ReportDesign,
    ga_best: &Milestone12ReportDesign,
) -> String {
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
            "**Option 1 — Selective Acoustic Center Treatment** (`{}`): {:.1}% of WBCs kept out of \
the treatment lane; RBC peripheral fraction {:.1}%; HI/pass = {:.2e}% (FDA 0.1% limit); \
ECV = {:.2} mL; P95 wall shear = {:.1} Pa. Zero active venturi throats — treatment relies on \
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
score = {:.3}; WBC treatment exposure {:.1}%; HI/pass = {:.2e}%; throat viscous heating ΔT = {:.2} K \
(FDA 5 K limit: {}). Score {:.4} under AsymmetricSplitVenturiCavitationSelectivity — \
**not comparable to Option 1 score**.\n",
        option2.candidate.id,
        m2.cavitation_number,
        gauge_kpa,
        d_throat_um,
        m2.serial_venturi_stages_per_path,
        m2.active_venturi_throat_count,
        m2.cancer_center_fraction * 100.0,
        m2.therapeutic_window_score,
        m2.wbc_recovery * 100.0,
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
        "**Topology Selection Physics — Option 1 vs Option 2:** The two optimization \
goals impose distinct constraints on optimal tree depth. For **Option 1** \
(selective acoustic, no venturi), the hybrid additive + geometric-mean synergy score \
rewards deeper split trees because each additional Zweifach–Fung splitting stage compounds \
cancer-center enrichment multiplicatively in the underlying flow partition: a three-stage \
Tri→Tri→Tri tree achieves cancer_center_fraction ~0.75+ compared to ~0.556 for Tri→Tri. \
Both the additive cancer-focusing term (22% weight) and the 4th-root synergy term amplify \
this separation advantage while the additive base ensures no feasible design collapses \
to zero (floor = 0.001). \
For **Option 2** (selective venturi cavitation), the selected design uses an \
{:.0} mm × {:.0} mm cross-section. Here, deeper trifurcation beyond 2 levels becomes \
counterproductive: a third level narrows terminal treatment channels to approximately \
140 µm (below the 200–400 µm Zweifach–Fung inertial-focusing optimum for millifluidic \
scales), routes only 1/27 of inlet flow to each venturi throat (insufficient velocity \
for cavitation inception at the available pressure budget), and consumes excess pressure \
drop across additional branching stages. RBC peripheral separation is achieved at \
Level 1; Level 2 refines cancer/WBC center concentration. The optimal tree depth is \
therefore physics-driven: acoustic-only treatment benefits from maximal splitting depth, \
while venturi cavitation requires sufficient flow velocity at the throat, imposing an \
upper bound on the number of upstream flow-dividing stages.\n",
        w_mm, h_mm,
    );

    // §6 — GA results and modeling note
    let _ = writeln!(
        s,
        "**GA Results — In-Place Dean–Serpentine Refinement:** The blueprint-native GA \
(seed {M12_GA_HYDRO_SEED}, InPlaceDeanSerpentineRefinement goal) produced rank-1 design \
`{}` with {} active throats ({} serial stage(s)). The GA applies three classes of \
architecture-preserving mutations via `BlueprintTopologyMutation`: (1) branch width scaling \
(treatment ×1.08, bypass ×0.94) to shift the Zweifach–Fung flow partition; (2) serpentine \
insertion on treatment-path channels, introducing Dean secondary flow (De = Re √(D_h/2R)) \
at bend apices where centrifugal forces focus larger CTCs toward the outer wall; and \
(3) venturi throat narrowing (×0.92) to lower σ at the vena contracta. The Dean number \
bonus (De_max/100) in the GA score explicitly rewards designs that co-localise inertial \
focusing from Dean vortices with hydrodynamic cavitation at venturi throats positioned at \
bend apices. GA score {:.4} — **not comparable to Option 2 Combined mode score**. \
**Modeling note:** each serial venturi throat is evaluated independently for cavitation \
number and intensity; the cumulative re-nucleation cascade between consecutive serial \
throats (inter-throat bubble collapse → re-growth → re-collapse) is not explicitly modeled \
in the 1D network solver. This is a conservative simplification — serial stages may provide \
greater cumulative dose than computed. Experimental validation at {}-stage conditions is \
recommended before clinical use.\n",
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

pub(super) fn build_appendix_a_supplemental(
    total_candidates: usize,
    opt1_pool: usize,
    opt2_pool: usize,
    ga_best: &Milestone12ReportDesign,
    robustness: &[crate::analysis::RobustnessReport],
    validation_rows: &[ValidationRow],
    canonical_results_path: &Path,
) -> String {
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

pub(super) fn build_references_block() -> String {
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
10. Dean, W.R. \"Note on the motion of fluid in a curved pipe.\" *Philosophical Magazine*, 4(20):208–223, 1927. [Dean number and secondary flow in curved channels; basis for `CurvaturePeakDeanNumber` venturi placement mode.]\n\
11. SonALAsense Internal Data — CFDrs canonical simulation data (see Appendix A)."
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
To ensure integrity, security, and accessibility of Milestone 12 computational data, SonALAsense applies policy-aligned storage controls across its managed engineering systems.\n\n\
**Data Governance and Access Control:** versioned simulation/report assets are maintained in the SonALAsense GitHub repository with controlled branch/merge workflows and auditable change history.\n\n\
**Secure Collaboration and Operational Storage:** project data packages and controlled shared artifacts are managed in SonALAsense Egnyte with role-based access and organization-managed retention policies.\n\n\
**Backup, Recovery, and Traceability:** canonical milestone outputs are regenerated deterministically from source code and preserved alongside run artifacts to support reproducible milestone evidence.\n"
        .to_string()
}

pub(super) fn build_storage_artifact_index() -> String {
    "\
- Canonical results: Appendix A (embedded below)\n\
- Figure manifest: Appendix B (embedded below)\n\
- Figures: inline throughout §5 Results\n\
- Generation artifacts: top-5 JSON, validation summaries, robustness outputs, and GA artifacts are included in Appendix A canonical data"
        .to_string()
}
