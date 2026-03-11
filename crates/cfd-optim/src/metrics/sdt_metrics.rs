//! The [`SdtMetrics`] output struct — all physics-derived metrics for one
//! [`BlueprintCandidate`](crate::design::BlueprintCandidate).

use serde::{Deserialize, Deserializer, Serialize};

/// Per-channel hemolysis decomposition for a single channel segment.
///
/// Each entry records the local Giersiepen HI contribution, shear conditions,
/// and flow fraction for one channel in the solved 1D network.  This enables
/// identification of which arm (treatment vs bypass) dominates hemolytic damage.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChannelHemolysis {
    /// Channel identifier (e.g. `"center_lv0"`, `"L_lv0"`, `"throat"`).
    pub channel_id: String,

    /// Whether this channel is a venturi throat segment.
    pub is_venturi_throat: bool,

    /// Flow-weighted Giersiepen HI contribution from this segment.
    pub hi_contribution: f64,

    /// Wall shear stress in this segment [Pa].
    pub wall_shear_pa: f64,

    /// Transit time through this segment [s].
    pub transit_time_s: f64,

    /// Fraction of total inlet flow carried by this channel.
    pub flow_fraction: f64,
}

/// All physics-derived metrics for one [`BlueprintCandidate`](crate::design::BlueprintCandidate).
#[derive(Debug, Clone, Default, Serialize, Deserialize)]
pub struct SdtMetrics {
    // ── Cavitation ──
    /// Cavitation number σ at the primary venturi throat.
    /// `< 1` implies active cavitation inception; lower values = more intense.
    /// `f64::INFINITY` when the topology has no venturi.
    #[serde(deserialize_with = "deserialize_f64_null_as_infinity")]
    pub cavitation_number: f64,

    /// Dimensionless cavitation potential: `max(0, 1 − σ)`.
    /// `0` = no cavitation; approaching `1` = very strong cavitation.
    pub cavitation_potential: f64,

    /// Wall shear rate at the venturi throat [1/s].
    pub throat_shear_rate_inv_s: f64,

    /// Estimated wall shear stress at the venturi throat [Pa].
    /// May greatly exceed the FDA 150 Pa limit; transit time is < 1 μs–1 ms.
    pub throat_shear_pa: f64,

    /// Whether the throat shear stress exceeds FDA 150 Pa.
    /// Expected `true` for effective SDT; cumulative haemolysis is tracked
    /// separately via [`hemolysis_index_per_pass`].
    pub throat_exceeds_fda: bool,

    // ── Main-channel safety ──
    /// Maximum wall shear stress in sustained-contact (non-throat) channels [Pa].
    /// Must be < 150 Pa for FDA compliance.
    pub max_main_channel_shear_pa: f64,

    /// FDA compliance of main distribution / serpentine channels.
    /// `true` iff `max_main_channel_shear_pa < 150 Pa`.
    pub fda_main_compliant: bool,

    // ── Cumulative haemolysis ──
    /// Unadjusted (bulk-flow) Giersiepen haemolysis index per pass.
    ///
    /// This is the raw value prior to selective venturi-composition correction.
    #[serde(default)]
    pub bulk_hemolysis_index_per_pass: f64,

    /// Giersiepen (1990) haemolysis index per device pass (dimensionless fraction).
    /// Sum of throat + main-channel contributions after selective venturi
    /// composition correction.
    pub hemolysis_index_per_pass: f64,

    // ── Distribution / exposure ──
    /// Flow uniformity across parallel outlet branches: `1 − CV(Q_outlets)`.
    /// `1.0` = perfectly uniform.  Symmetric designs reach 1.0 by construction.
    pub flow_uniformity: f64,

    /// Fraction of the 36 treatment wells covered by channel passes.
    /// `1.0` = all 36 wells are in the channel path.
    pub well_coverage_fraction: f64,

    /// Mean fluid residence time in the treatment zone [s].
    pub mean_residence_time_s: f64,

    // ── System-level ──
    /// Total pressure drop across the device [Pa].
    pub total_pressure_drop_pa: f64,

    /// Total channel path length [mm].
    pub total_path_length_mm: f64,

    /// `true` if the total ΔP < available gauge pressure (syringe pump can
    /// drive the flow without stalling).
    pub pressure_feasible: bool,

    // ── Cell separation (2-population: cancer vs RBC) ──
    /// Inertial focusing separation efficiency for MCF-7 cancer cells vs RBCs.
    ///
    /// Defined as `|x̃_cancer − x̃_rbc|` ∈ [0, 1] where `x̃` is the normalised
    /// lateral equilibrium position (0 = center, 1 = wall).
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub cell_separation_efficiency: f64,

    /// Fraction of MCF-7 cancer cells collected in the center channel.
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub cancer_center_fraction: f64,

    /// Fraction of RBCs collected in the peripheral (wall) channels.
    ///
    /// `0.0` for topologies without a cell separation stage.
    pub rbc_peripheral_fraction: f64,

    // ── Three-population separation (WBC+cancer→center, RBC→periphery) ──
    /// Three-population separation efficiency.
    ///
    /// Defined as `x̃_rbc − max(x̃_cancer, x̃_wbc)`, quantifying how well
    /// RBCs are pushed to the wall relative to both WBCs and cancer cells.
    /// Higher values mean both WBC and cancer stay central while RBCs migrate
    /// to the wall.
    ///
    /// `0.0` for topologies without a three-population separation stage.
    pub three_pop_sep_efficiency: f64,

    /// Fraction of WBCs collected in the center channel.
    ///
    /// `0.0` for topologies without a three-population separation stage.
    pub wbc_center_fraction: f64,

    /// Fraction of RBCs collected in the peripheral channels in the
    /// three-population separation topology.
    ///
    /// `0.0` for topologies without a three-population separation stage.
    pub rbc_peripheral_fraction_three_pop: f64,

    /// Normalised lateral equilibrium position of WBCs (0=center, 1=wall).
    ///
    /// `0.5` (dispersed) for topologies without a separation stage.
    pub wbc_equilibrium_pos: f64,

    /// Normalised lateral equilibrium position of cancer cells (0=center, 1=wall).
    ///
    /// `0.5` (dispersed) for topologies without a separation stage.
    pub cancer_equilibrium_pos: f64,

    /// Normalised lateral equilibrium position of RBCs (0=center, 1=wall).
    ///
    /// `0.5` (dispersed) for topologies without a separation stage.
    pub rbc_equilibrium_pos: f64,

    // ── Leukapheresis-style parallel or recirculating blueprints ──
    /// Fraction of WBCs focused into the center outlet channel.
    ///
    /// Computed using [`enhanced_lateral_equilibrium`] with CFL correction at
    /// the candidate's `feed_hematocrit`.  `0.0` for non-leukapheresis topologies.
    pub wbc_recovery: f64,

    /// Fraction of RBCs that pass into the center outlet channel (want **low**).
    ///
    /// Low `rbc_pass_fraction` means better RBC removal.
    /// `0.0` for non-leukapheresis topologies.
    pub rbc_pass_fraction: f64,

    /// WBC purity in the center outlet: `WBC_center / (WBC_center + RBC_center)`.
    ///
    /// `0.0` for non-leukapheresis topologies.
    pub wbc_purity: f64,

    /// Total extracorporeal volume [mL] (sum of all channel volumes).
    ///
    /// Used to verify that ECV ≤ 10% of patient blood volume.
    /// `0.0` for non-leukapheresis topologies.
    pub total_ecv_ml: f64,

    // ── Candidate flow rate (convenience) ──
    /// Candidate volumetric flow rate in mL/min.
    ///
    /// Converted from `flow_rate_m3_s` (m³/s → mL/min) for use in scoring
    /// throughput against clinical targets (e.g. 10 mL/min leukapheresis).
    pub flow_rate_ml_min: f64,

    // ── Plate-fit and outlet-port count ──
    /// `true` if the estimated channel-layout footprint fits within the
    /// 96-well plate footprint (127.76 × 85.47 mm).
    ///
    /// All symmetric millifluidic topologies fit within the 45 × 45 mm treatment
    /// zone by construction. Parallel-path layouts can violate the plate
    /// boundary when `n_channels × pitch > plate_height`.
    pub plate_fits: bool,

    /// Number of external outlet ports for this topology.
    ///
    /// * `1` — SDT / exposure topologies: all branches merge back to a single
    ///   outlet via the outlet trunk segment.
    /// * `2` — Cell-separation and leukapheresis topologies: centre outlet
    ///   (WBC / cancer collection) and peripheral outlet (RBC removal).
    pub n_outlet_ports: usize,

    /// Platelet activation index per device pass (Hellums 1994 power-law model).
    ///
    /// `PAI = 1.8×10⁻⁸ × τ^1.325 × t^0.462`
    /// where τ = peak shear stress [Pa], t = exposure duration [s].
    /// Threshold: [`PAI_PASS_LIMIT`](crate::constraints::PAI_PASS_LIMIT) = 5×10⁻⁴.
    pub platelet_activation_index: f64,

    /// Estimated main-channel shear rate [1/s].
    ///
    /// Computed from the main-channel shear stress and reference blood viscosity:
    /// `γ̇_main ≈ τ_main / μ_ref`.
    #[serde(default)]
    pub main_channel_shear_rate_inv_s: f64,

    /// Low-flow stasis contribution to clotting risk ∈ [0, 1].
    ///
    /// Derived from blood flow rate relative to literature-informed extracorporeal
    /// operating bands. Higher values indicate greater stasis-driven clotting risk.
    #[serde(default)]
    pub low_flow_stasis_risk: f64,

    /// Minimum per-channel shear-rate stasis risk ∈ [0, 1].
    ///
    /// Uses the *minimum* shear rate across all flow-carrying channels (excluding
    /// zero-flow segments). Identifies the worst stasis zone in the device rather
    /// than the overall max-shear regime. Wider or slower channels produce lower
    /// minimum shear rates and higher stasis risk.
    #[serde(default)]
    pub min_shear_stasis_risk: f64,

    /// Maximum per-channel residence-time stasis risk ∈ [0, 1].
    ///
    /// Uses the longest individual channel transit time rather than the mean chip
    /// residence. The slowest flow path determines the stasis bottleneck and the
    /// greatest local fibrin deposition potential.
    #[serde(default)]
    pub max_residence_stasis_risk: f64,

    /// Post-expansion recirculation stasis risk ∈ [0, 1].
    ///
    /// Area expansion from venturi throat to downstream diffuser creates
    /// Borda-Carnot recirculation eddies. Risk scales with the log of the
    /// expansion area ratio (Idelchik 1994, Diagram 4-1; Borda-Carnot loss
    /// coefficient K_BC = (1 − A_small/A_large)²).
    #[serde(default)]
    pub expansion_stasis_risk: f64,

    /// Dead-volume stasis risk ∈ [0, 1].
    ///
    /// Fraction of total chip channel volume residing in channels where local
    /// wall shear rate falls below the platelet-adhesion threshold (200 /s;
    /// Folie & McIntire 1989). Higher dead-volume fraction indicates more
    /// stagnant blood volume susceptible to thrombus formation.
    #[serde(default)]
    pub dead_volume_stasis_risk: f64,

    /// Composite clotting risk index ∈ [0, 1].
    ///
    /// Geometry-sensitive weighted blend of five stasis sub-terms:
    ///
    /// `clotting_risk_index = 0.20·low_flow + 0.20·min_shear + 0.20·max_residence
    ///                      + 0.20·expansion + 0.20·dead_volume`
    ///
    /// Unlike the earlier 3-term model (which was degenerate for candidates at
    /// the same flow rate), every sub-term depends on per-channel solved-flow
    /// geometry, producing design-discriminating clot risk even within a single
    /// operating point.
    #[serde(default)]
    pub clotting_risk_index: f64,

    /// `true` if blood flow is at/above the clotting-caution threshold.
    ///
    /// Threshold uses [`CLOTTING_BFR_CAUTION_ML_MIN`](crate::constraints::CLOTTING_BFR_CAUTION_ML_MIN).
    #[serde(default)]
    pub clotting_flow_compliant: bool,

    /// Conservative low-flow clotting risk assuming a strict 10 mL/s flow claim.
    ///
    /// Uses the same risk blend as `clotting_risk_index`, but replaces the flow
    /// term with a strict-flow sensitivity anchor at 600 mL/min.
    #[serde(default)]
    pub clotting_risk_index_10ml_s: f64,

    /// `true` if blood flow is at/above 10 mL/s (600 mL/min).
    ///
    /// This field supports sensitivity analysis for strict low-flow clotting
    /// assumptions and is reported separately from default compliance gates.
    #[serde(default)]
    pub clotting_flow_compliant_10ml_s: bool,

    /// True when the treatment zone uses venturi throats (hydrodynamic SDT).
    #[serde(default)]
    pub venturi_treatment_enabled: bool,

    /// Human-readable treatment-zone actuation mode.
    ///
    /// Values: `"UltrasoundOnly"` or `"VenturiThroats"`.
    #[serde(default)]
    pub treatment_zone_mode: String,

    /// Total active venturi throat count across the device.
    #[serde(default)]
    pub active_venturi_throat_count: usize,

    /// Number of serial venturi throats per treatment path.
    #[serde(default)]
    pub serial_venturi_stages_per_path: usize,

    /// Fraction of total inlet flow that traverses venturi throat sections.
    ///
    /// Usually `1.0`; less than `1.0` for selective-separation topologies with
    /// bypass arms.
    #[serde(default = "default_one")]
    pub venturi_flow_fraction: f64,

    /// Model-predicted center-stream flow fraction for staged primitive selective split trees.
    ///
    /// `q_model = q_center(center_frac)^n_levels`
    #[serde(default)]
    pub cct_model_venturi_flow_fraction: f64,

    /// Solved-network venturi flow fraction for staged primitive selective split trees.
    ///
    /// Non-selective topologies keep this at `0.0`.
    #[serde(default)]
    pub cct_solved_venturi_flow_fraction: f64,

    /// Mean solved center-flow split across staged primitive trifurcation segments.
    ///
    /// Non-selective topologies keep this at `0.0`.
    #[serde(default)]
    pub cct_stage_center_qfrac_mean: f64,

    /// Model-predicted venturi flow fraction for staged CIF topologies.
    ///
    /// `q_model = q_pretri^n_pretri × q_terminal_tri × q_bi_treat`
    #[serde(default)]
    pub cif_model_venturi_flow_fraction: f64,

    /// Solved-network venturi flow fraction for staged CIF topologies.
    ///
    /// Non-CIF topologies keep this at `0.0`.
    #[serde(default)]
    pub cif_solved_venturi_flow_fraction: f64,

    /// Mean solved center-flow split across CIF pre-trifurcation stages.
    ///
    /// Non-CIF topologies keep this at `0.0`.
    #[serde(default)]
    pub cif_pretri_qfrac_mean: f64,

    /// Solved terminal-trifurcation center-flow fraction for CIF.
    ///
    /// Falls back to the model value when the solved split cannot be extracted.
    #[serde(default)]
    pub cif_terminal_tri_qfrac: f64,

    /// Solved terminal-bifurcation treatment-flow fraction for CIF.
    ///
    /// Falls back to the model value when the solved split cannot be extracted.
    #[serde(default)]
    pub cif_terminal_bi_qfrac: f64,

    /// CIF post-remerge outlet-tail length [mm].
    ///
    /// Shorter values indicate "remerge near outlet" layouts.
    #[serde(default)]
    pub cif_outlet_tail_length_mm: f64,

    /// Fraction of inlet RBC population that traverses the venturi path.
    ///
    /// Lower values indicate stronger RBC protection via peripheral skimming.
    #[serde(default = "default_one")]
    pub rbc_venturi_exposure_fraction: f64,

    /// Effective hematocrit at the venturi throat.
    ///
    /// For primitive selective split-sequence designs,
    /// this is lower than the bulk feed hematocrit because the Zweifach-Fung
    /// cascade routes RBCs preferentially to the peripheral bypass arms.
    /// HI and PAI are computed using this LOCAL concentration at the venturi.
    ///
    /// Equals `candidate.feed_hematocrit` for all other topologies.
    pub local_hematocrit_venturi: f64,

    /// Fraction of cancer cells that pass through the venturi throat and
    /// receive full cavitation / SDT dose.
    ///
    /// `cancer_dose_fraction = cancer_center_fraction × cavitation_potential`
    ///
    /// Non-zero for selective-separation topologies where the venturi is placed
    /// on the enriched treatment arm.
    pub cancer_dose_fraction: f64,

    // ── Cancer-targeted hydrodynamic cavitation SDT metrics ──────────────────
    /// Enhanced cavitation intensity combining cavitation depth and throat constriction.
    ///
    /// `cavitation_intensity = cavitation_potential × (0.5 + 0.5 × constriction_score)`
    ///
    /// where `constriction_score = clamp(ln(v_throat / v_inlet) / ln(80), 0, 1)`.
    /// A score of 1.0 requires both active cavitation (σ < `σ_crit`) and a narrow
    /// throat (50 µm in a 4 mm channel gives the reference constriction of 80×).
    /// Zero for topologies without a venturi.
    #[serde(default)]
    pub cavitation_intensity: f64,

    /// Cancer-targeted cavitation: fraction of cancer cells receiving cavitation dose,
    /// weighted by cavitation intensity.
    ///
    /// `cancer_targeted_cavitation = cancer_center_fraction × cavitation_intensity`
    ///
    /// Maximised when cancer cells are strongly enriched in the venturi center
    /// arm AND the throat constriction generates intense cavitation.
    #[serde(default)]
    pub cancer_targeted_cavitation: f64,

    /// WBC-targeted cavitation: fraction of WBCs receiving cavitation dose,
    /// weighted by cavitation intensity.
    ///
    /// `wbc_targeted_cavitation = wbc_center_fraction × cavitation_intensity`
    #[serde(default)]
    pub wbc_targeted_cavitation: f64,

    /// RBC venturi protection: RBCs shielded from the cavitating core.
    ///
    /// `rbc_venturi_protection = rbc_peripheral_fraction_three_pop
    ///                           × (1 − cavitation_intensity × rbc_venturi_exposure_fraction)`
    ///
    /// High values indicate most RBCs bypass the venturi throat entirely via
    /// peripheral skimming arms, avoiding haemolysis from cavitation collapse.
    #[serde(default)]
    pub rbc_venturi_protection: f64,

    /// Sonoluminescence proxy based on adiabatic Rayleigh-Plesset collapse temperature.
    ///
    /// `T_collapse / T_blood ≈ (P_abs / P_vapor)^((γ−1)/γ)`
    ///
    /// Normalised so that 3 bar gauge at full cavitation → 1.0.
    /// Higher values indicate more energetic bubble collapse and stronger
    /// sonosensitiser activation.  Zero when no cavitation occurs.
    #[serde(default)]
    pub sonoluminescence_proxy: f64,

    // ── RBC safety / FDA transit-time analysis ───────────────────────────────
    /// Transit time of blood through the venturi throat [s].
    ///
    /// Computed as `throat_length_m / v_throat`.
    /// Near-zero when the topology has no venturi.  Used to assess whether the
    /// FDA transit-time exception for high shear applies.
    #[serde(default)]
    pub throat_transit_time_s: f64,

    /// Overall FDA compliance considering the throat transit-time exception.
    ///
    /// `true` if:
    /// - `fda_main_compliant` is `true` (main channels ≤ 150 Pa), **AND**
    /// - Either the throat shear is also ≤ 150 Pa, **OR** the throat transit
    ///   time is < 5 ms **and** throat shear is < 300 Pa (brief transient
    ///   exemption per published VAD device-exemption practice).
    ///
    /// Use this as the primary FDA compliance indicator in clinical reports.
    #[serde(default)]
    pub fda_overall_compliant: bool,

    /// Composite lysis risk index (dimensionless).
    ///
    /// `lysis_risk_index = hemolysis_index_per_pass
    ///                   × (1 + 5 × rbc_venturi_exposure_fraction × local_hematocrit_venturi)`
    ///
    /// The amplification factor accounts for additional haemolysis caused by
    /// cavitation-bubble collapse in the venturi throat when RBCs are present.
    /// A factor of 5× is a conservative estimate based on acoustic micro-streaming
    /// and shockwave-induced lysis literature.
    ///
    /// Lower is safer; clinical target: < 0.001.
    #[serde(default)]
    pub lysis_risk_index: f64,

    /// Therapeutic window score ∈ [0, 1].
    ///
    /// `therapeutic_window_score
    ///   = cancer_targeted_cavitation / (1×10⁻⁶ + lysis_risk_index)
    ///     / THERAPEUTIC_WINDOW_REF`
    ///
    /// Quantifies cancer treatment intensity relative to RBC lysis risk.
    /// Higher values indicate designs that effectively treat cancer cells while
    /// protecting red cells.  Normalised so that a strong cavitation design
    /// (cancer_targeted_cav = 0.5) with clinical-target lysis (index = 0.001)
    /// yields a score of ≈ 1.0.
    #[serde(default)]
    pub therapeutic_window_score: f64,

    /// Oncology selectivity index ∈ [0, 1].
    ///
    /// Captures how strongly cavitation is concentrated onto cancer-enriched
    /// flow while shielding RBCs from venturi exposure:
    ///
    /// `oncology_selectivity_index = cancer_targeted_cavitation
    ///   × (1 − rbc_venturi_exposure_fraction)
    ///   × (0.5 + 0.5 × cancer_center_fraction)`
    ///
    /// Higher is better. This metric is used by SDT ranking modes that target
    /// preferential cancer-cell treatment.
    #[serde(default)]
    pub oncology_selectivity_index: f64,

    /// Cancer-vs-RBC cavitation bias index ∈ [0, 1].
    ///
    /// Compares cavitation load routed to cancer-enriched flow against
    /// cavitation load seen by RBC-rich venturi flow:
    ///
    /// `bias = cancer_load / (cancer_load + rbc_load)`
    ///
    /// with:
    /// - `cancer_load = cancer_targeted_cavitation`
    /// - `rbc_load = rbc_venturi_exposure_fraction × local_hematocrit_venturi × cavitation_intensity`
    ///
    /// Higher is better (more cavitation concentrated in cancer-targeted stream).
    #[serde(default)]
    pub cancer_rbc_cavitation_bias_index: f64,

    /// CIF remerge-proximity score ∈ [0, 1].
    ///
    /// Quantifies "remerge near outlet" behavior for staged CIF layouts using
    /// `cif_outlet_tail_length_mm`:
    /// - short post-remerge tails → score near `1`
    /// - long post-remerge tails → score near `0`
    ///
    /// Non-CIF topologies keep this at `0.0`.
    #[serde(default)]
    pub cif_remerge_proximity_score: f64,

    /// Composite selective cavitation delivery index ∈ [0, 1].
    ///
    /// Blends oncology selectivity, cancer-vs-RBC cavitation bias, and (for CIF)
    /// near-outlet remerge preference:
    ///
    /// `selective_delivery = oncology_selectivity_index × cancer_rbc_cavitation_bias_index × remerge_bonus`
    ///
    /// where `remerge_bonus = 0.5 + 0.5 × cif_remerge_proximity_score` for CIF,
    /// and `1.0` for non-CIF topologies.
    #[serde(default)]
    pub selective_cavitation_delivery_index: f64,

    /// Projected RBC haemolysis rate [% Hb release per hour] at operating flow.
    ///
    /// `rbc_lysis_rate_pct_per_h
    ///   = hemolysis_index_per_pass × flow_rate_ml_min × 60 / 5000 × 100`
    ///
    /// Assumes a 5 L adult blood volume.  Values < 0.1 %/h are considered
    /// acceptable for extended clinical operation (continuous-flow devices).
    #[serde(default)]
    pub rbc_lysis_rate_pct_per_h: f64,

    /// Estimated number of blood-volume passes in 15 minutes for a 3 kg patient.
    ///
    /// `passes_15m = (flow_rate_ml_min × 15) / (3 × 85)`
    #[serde(default)]
    pub projected_passes_15min_pediatric_3kg: f64,

    /// Projected cumulative hemolysis fraction over 15 minutes for a 3 kg patient.
    ///
    /// `HI_15m = 1 − (1 − HI_pass)^passes_15m`
    #[serde(default)]
    pub projected_hemolysis_15min_pediatric_3kg: f64,

    /// Projected cumulative hemolysis fraction over 15 minutes for a 5 L adult blood volume.
    ///
    /// `HI_15m = 1 − (1 − HI_pass)^passes_15m`
    #[serde(default)]
    pub projected_hemolysis_15min_adult: f64,

    /// Fraction of cancer cells in the active therapy zone at steady state.
    ///
    /// `cancer_therapy_zone_fraction = cancer_center_fraction × venturi_flow_fraction`
    ///
    /// Represents the instantaneous fraction of cancer cells receiving SDT dose
    /// in the cavitating venturi stream.  Maximised when cancer cells are strongly
    /// enriched in the center arm AND that arm carries a large flow fraction.
    #[serde(default)]
    pub cancer_therapy_zone_fraction: f64,

    /// Effective optical path length for 405 nm light transport through blood [m].
    ///
    /// Approximated from treatment-channel depth in the current 1D model.
    #[serde(default)]
    pub optical_path_length_405_m: f64,

    /// Relative 405 nm light-delivery index ∈ [0, 1].
    ///
    /// Beer-Lambert-style proxy:
    /// `exp(-μ_eff,405 · optical_path_length_405_m · hematocrit_scale)`.
    ///
    /// Higher values indicate stronger blue-light delivery to circulating cells.
    #[serde(default)]
    pub blue_light_delivery_index_405nm: f64,

    /// FDA main-channel safety margin [Pa].
    ///
    /// `safety_margin_pa = 150 Pa − max_main_channel_shear_pa`
    ///
    /// Positive values indicate compliant designs with headroom to the 150 Pa
    /// limit.  Negative values indicate a main-channel shear violation.
    /// Used in the milestone report to rank designs by how far below the limit
    /// they operate.
    #[serde(default)]
    pub safety_margin_pa: f64,

    /// Fraction of total channel path length in the active therapy zone ∈ [0, 1].
    ///
    /// Model-based approximation of the fraction of the chip's channel network
    /// that is doing therapeutic work (cancer-targeted or venturi-throat channels)
    /// vs. bypass plumbing or healthy-cell routing arms.
    ///
    /// - `1.0` — all-flow treatment topologies (single/double venturi, etc.).
    /// - `< 1.0` — selective-routing designs (primitive split trees, asymmetric bifurcation, etc.)
    ///   where a fraction of flow bypasses the therapy zone to protect RBCs.
    /// - `0.0` — non-venturi or leukapheresis-only designs.
    ///
    /// Larger values indicate more of the 45 × 45 mm treatment zone is used for
    /// direct SDT therapy.
    #[serde(default)]
    pub therapy_channel_fraction: f64,

    // ── Wall shear percentile statistics (FDA spatial assessment) ────────────
    /// 95th-percentile wall shear stress across all non-throat channels [Pa].
    ///
    /// Per ASTM F1841-20: spatial distribution of shear must be reported.
    /// P95 captures the predominant shear environment excluding brief peaks.
    #[serde(default)]
    pub wall_shear_p95_pa: f64,

    /// 99th-percentile wall shear stress across all non-throat channels [Pa].
    #[serde(default)]
    pub wall_shear_p99_pa: f64,

    /// Mean wall shear stress across all non-throat channels [Pa].
    #[serde(default)]
    pub wall_shear_mean_pa: f64,

    /// Coefficient of variation of wall shear stress across non-throat channels.
    ///
    /// `CV = std_dev / mean`.  Lower CV indicates more uniform shear exposure.
    #[serde(default)]
    pub wall_shear_cv: f64,

    /// FDA percentile-compliant: P95 ≤ 150 Pa AND P99 ≤ 300 Pa.
    #[serde(default)]
    pub fda_shear_percentile_compliant: bool,

    // ── Diffuser pressure recovery ───────────────────────────────────────────
    /// Diffuser pressure recovery downstream of venturi throat [Pa].
    ///
    /// `ΔP_recovery = C_D × ½ρ(v_throat² − v_outlet²)` (Idelchik Diagram 6-21).
    /// A non-zero value indicates the abrupt-expansion assumption has been
    /// replaced with a gradual-expansion diffuser model.
    #[serde(default)]
    pub diffuser_recovery_pa: f64,

    // ── Acoustic energy budget ───────────────────────────────────────────────
    /// Total hydraulic pumping power [W].
    ///
    /// `mechanical_power_w = total_pressure_drop_pa × Q_inlet`
    ///
    /// Represents the baseline mechanical energy supplied by the syringe pump.
    #[serde(default)]
    pub mechanical_power_w: f64,

    /// Fraction of pumping power directed into cavitation zone (dimensionless proxy).
    ///
    /// `acoustic_capture_efficiency = cavitation_potential × (venturi_shear_pa / P_atm)`
    ///
    /// Capped at 1.0.  Zero for non-venturi topologies.
    #[serde(default)]
    pub acoustic_capture_efficiency: f64,

    /// Specific cavitation energy per mL of blood processed [mJ/mL].
    ///
    /// `specific_cavitation_energy_j_ml =
    ///     acoustic_capture_efficiency × mechanical_power_w
    ///     × mean_residence_time_s / (flow_rate_ml_min / 60_000) × 1000`
    #[serde(default)]
    pub specific_cavitation_energy_j_ml: f64,

    /// Cavitation-amplified haemolysis index per pass (dimensionless).
    ///
    /// Applies `cavitation_amplified_hi` from `cfd_1d::hemolysis`:
    /// `HI_amplified = HI_base × (1 + 3 × cav_potential.clamp(0, 1))`
    ///
    /// For non-venturi topologies equals `hemolysis_index_per_pass`.
    #[serde(default)]
    pub hemolysis_index_per_pass_cavitation_amplified: f64,

    // ── Per-channel hemolysis decomposition ──────────────────────────────────
    /// HI contribution from venturi/treatment channel(s) only.
    ///
    /// After local-hematocrit correction for selective split-sequence topologies.
    /// Zero for non-venturi topologies.
    #[serde(default)]
    pub treatment_channel_hi: f64,

    /// Mean HI contribution across bypass (non-venturi) channels.
    ///
    /// Captures the average hemolytic load in peripheral/skimming arms.
    #[serde(default)]
    pub bypass_channel_hi_mean: f64,

    /// Maximum HI contribution among bypass (non-venturi) channels.
    ///
    /// Identifies the worst-case bypass arm for safety assessment.
    #[serde(default)]
    pub bypass_channel_hi_max: f64,

    /// Per-channel hemolysis breakdown for all channels in the solved network.
    ///
    /// Enables detailed reporting of which channel segments contribute most
    /// to cumulative hemolysis.  Sorted by descending HI contribution.
    #[serde(default)]
    pub per_channel_hemolysis: Vec<ChannelHemolysis>,

    #[serde(default)]
    pub acoustic_resonance_factor: f64,

    #[serde(default)]
    pub channel_resonance_score: f64,

    #[serde(default)]
    pub serial_cavitation_dose_fraction: f64,

    #[serde(default)]
    pub treatment_zone_dwell_time_s: f64,

    #[serde(default)]
    pub throat_temperature_rise_k: f64,

    #[serde(default)]
    pub fda_thermal_compliant: bool,
}

/// Deserialize an `f64` that may be JSON `null` (produced when `f64::INFINITY`
/// is serialized by serde_json).  Maps `null` → `f64::INFINITY`.
fn deserialize_f64_null_as_infinity<'de, D: Deserializer<'de>>(d: D) -> Result<f64, D::Error> {
    Ok(Option::<f64>::deserialize(d)?.unwrap_or(f64::INFINITY))
}

fn default_one() -> f64 {
    1.0
}
