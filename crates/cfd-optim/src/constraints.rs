//! Physical constants and design constraints for millifluidic SDT optimization.
//!
//! # Geometry reference
//! 96-well plates follow ANSI/SLAS 1-2004.
//! The 6 × 6 center treatment zone (rows B–G, columns D–I) is used as the
//! target area for uniform light / ultrasound exposure.
//!
//! # FDA hemolysis guidance
//! Maximum sustained wall shear stress for blood-contacting devices: **150 Pa**
//! (FDA Guidance, ASTM F1841).  Venturi throats may briefly exceed this limit
//! (microsecond transit time); the Giersiepen (1990) model is used to assess
//! cumulative hemolysis index rather than applying a blanket 150 Pa cutoff at
//! the throat.

// ── 96-well plate geometry (ANSI/SLAS 1-2004) ────────────────────────────────

/// Plate outer width [mm]
pub const PLATE_WIDTH_MM: f64 = 127.76;
/// Plate outer height [mm]
pub const PLATE_HEIGHT_MM: f64 = 85.47;
/// Well centre-to-centre pitch in both H and V directions [mm]
pub const WELL_PITCH_MM: f64 = 9.00;
/// A1 well centre, distance from left edge [mm]
pub const A1_X_MM: f64 = 14.38;
/// A1 well centre, distance from top edge [mm]
pub const A1_Y_MM: f64 = 11.24;

// ── 6 × 6 centre treatment zone ────────────────────────────────────────────

/// Number of rows in the treatment zone (rows B–G, 0-indexed 1–6)
pub const TREATMENT_ROWS: usize = 6;
/// Number of columns in the treatment zone (cols D–I, 0-indexed 3–8)
pub const TREATMENT_COLS: usize = 6;
/// Total well count in the treatment zone
pub const TREATMENT_WELL_COUNT: usize = TREATMENT_ROWS * TREATMENT_COLS; // 36

/// X coordinate of first treatment-zone well centre [mm]  (col D = index 3)
pub const TREATMENT_X_MIN_MM: f64 = 41.38; // A1_X_MM + 3 * WELL_PITCH_MM
/// X coordinate of last treatment-zone well centre [mm]   (col I = index 8)
pub const TREATMENT_X_MAX_MM: f64 = 86.38; // TREATMENT_X_MIN_MM + 5 * WELL_PITCH_MM
/// Y coordinate of first treatment-zone well centre [mm]  (row B = index 1)
pub const TREATMENT_Y_MIN_MM: f64 = 20.24; // A1_Y_MM + 1 * WELL_PITCH_MM
/// Y coordinate of last treatment-zone well centre [mm]   (row G = index 6)
pub const TREATMENT_Y_MAX_MM: f64 = 65.24; // TREATMENT_Y_MIN_MM + 5 * WELL_PITCH_MM
/// Treatment zone span in X [mm]
pub const TREATMENT_WIDTH_MM: f64 = 45.00;
/// Treatment zone span in Y [mm]
pub const TREATMENT_HEIGHT_MM: f64 = 45.00;
/// Treatment zone total area [mm²]
pub const TREATMENT_AREA_MM2: f64 = TREATMENT_WIDTH_MM * TREATMENT_HEIGHT_MM; // 2025

// ── Blood physical properties (37 °C) ─────────────────────────────────────

/// Blood density [kg / m³]   (Fung 1993)
pub const BLOOD_DENSITY_KG_M3: f64 = 1060.0;
/// Blood dynamic viscosity, high-shear Newtonian approximation [Pa · s]
/// Used for wall-shear quick-estimates; the non-Newtonian Casson model is used
/// for full physics analysis via `CassonBlood::normal_blood()`.
pub const BLOOD_VISCOSITY_PA_S: f64 = 3.45e-3;
/// Vapour pressure of blood at 37 °C [Pa]
pub const BLOOD_VAPOR_PRESSURE_PA: f64 = 6_280.0;

// ── FDA hemolysis limits ───────────────────────────────────────────────────

/// Maximum sustained wall shear stress: FDA guidance for blood-contacting devices [Pa]
pub const FDA_MAX_WALL_SHEAR_PA: f64 = 150.0;
/// Acceptable hemolysis index per device pass (Giersiepen 1990 model)
/// < 0.1 % fractional plasma haemoglobin increase
pub const HI_PASS_LIMIT: f64 = 0.001;

/// Extended shear stress limit for venturi throats with brief transit time [Pa].
///
/// FDA guidance allows higher peak shear when the exposure duration is below
/// [`FDA_TRANSIENT_TIME_S`].  A threshold of 300 Pa is used based on published
/// device-exemption practice for sharp-orifice venturi throats where red cell
/// transit is < 1–5 ms.
pub const FDA_TRANSIENT_SHEAR_PA: f64 = 300.0;

/// Maximum throat transit time to qualify for the FDA transient shear exception [s].
///
/// If the blood transit through the venturi throat is shorter than this value,
/// [`FDA_TRANSIENT_SHEAR_PA`] applies instead of [`FDA_MAX_WALL_SHEAR_PA`].
/// Based on FDA Guidance for VADs and published ASTM F1841 commentary.
pub const FDA_TRANSIENT_TIME_S: f64 = 5e-3; // 5 ms

/// Reference adult patient blood volume for clinical lysis-rate projection [mL].
///
/// Used to convert per-pass haemolysis index to a projected steady-state lysis
/// rate in % Hb release per hour at the device operating flow rate.
pub const PATIENT_BLOOD_VOLUME_ML: f64 = 5_000.0;

/// Blood-volume estimate [mL/kg] for pediatric clinical projections.
///
/// Used for 3 kg neonatal safety projections in milestone reporting.
pub const PEDIATRIC_BLOOD_VOLUME_ML_PER_KG: f64 = 85.0;

/// Reference pediatric weight [kg] used in milestone-12 cumulative-safety projections.
pub const PEDIATRIC_REFERENCE_WEIGHT_KG: f64 = 3.0;

/// Milestone-12 therapy window for cumulative hemolysis projection [min].
pub const MILESTONE_TREATMENT_DURATION_MIN: f64 = 15.0;

/// Normalisation reference for [`SdtMetrics::therapeutic_window_score`].
///
/// Chosen so that `cancer_targeted_cavitation = 0.5` with
/// `lysis_risk_index = 0.001` yields a score of 1.0:
/// `0.5 / (1e-6 + 0.001) / 500 ≈ 1.0`.
pub const THERAPEUTIC_WINDOW_REF: f64 = 500.0;

// ── Cavitation criterion ───────────────────────────────────────────────────

/// Standard atmospheric pressure [Pa]
pub const P_ATM_PA: f64 = 101_325.0;
/// Cavitation inception threshold (σ < `σ_crit` → cavitation)
pub const SIGMA_CRIT: f64 = 1.0;

// ── Giersiepen (1990) haemolysis model coefficients ───────────────────────
//   ΔHb/Hb = C_G · t^α · τ^β
//   Reference: Giersiepen et al., Int. J. Artif. Organs 13 (1990) 300–306

/// Pre-multiplied constant `C_G`
pub const GIERSIEPEN_C: f64 = 3.62e-5;
/// Time exponent α
pub const GIERSIEPEN_ALPHA: f64 = 0.765;
/// Shear stress exponent β
pub const GIERSIEPEN_BETA: f64 = 1.991;

// ── Design parameter sweep ranges ─────────────────────────────────────────

/// Volumetric flow rates to sweep [m³/s]  (100 – 600 mL/min — paediatric low-flow to high-throughput)
pub const FLOW_RATES_M3_S: [f64; 7] = [
    1.667e-6,  // 100 mL/min — paediatric / low-flow
    2.500e-6,  // 150 mL/min — light perfusion
    3.333e-6,  // 200 mL/min — dialysis low-end
    5.000e-6,  // 300 mL/min — dialysis mid
    6.667e-6,  // 400 mL/min — dialysis high-end
    8.333e-6,  // 500 mL/min — high-flow leukapheresis
    10.000e-6, // 600 mL/min — high-throughput SDT regime (NEW)
];

/// Gauge pressures at inlet [Pa]  (1, 2, 3, 4, 5 bar above atmospheric)
pub const INLET_GAUGES_PA: [f64; 5] = [100_000.0, 200_000.0, 300_000.0, 400_000.0, 500_000.0];

/// Venturi throat diameters [m]  (30 – 200 μm — extreme cavitation to moderate constriction)
pub const THROAT_DIAMETERS_M: [f64; 6] = [
    30e-6,  // 30 μm — extreme cavitation; σ < 1 at moderate flow (NEW)
    50e-6,  // 50 μm — strong cavitation
    75e-6,  // 75 μm — moderate-strong cavitation
    100e-6, // 100 μm — moderate cavitation
    150e-6, // 150 μm — mild cavitation
    200e-6, // 200 μm — gentle cavitation
];

/// Main channel widths for rectangular serpentine channels [m]  (0.5 – 8 mm).
///
/// Includes a thin 0.5 mm branch to capture 405 nm light-delivery scenarios
/// where optical penetration in whole blood is limited.
pub const CHANNEL_WIDTHS_M: [f64; 6] = [0.5e-3, 1.0e-3, 2.0e-3, 4.0e-3, 6.0e-3, 8.0e-3];

/// Fixed channel height used for all rectangular channels [m]
pub const CHANNEL_HEIGHT_M: f64 = 1.0e-3; // 1 mm — dialysis-compatible millifluidic

/// Fixed inlet/outlet port diameter for all venturi stages [m]  (4 mm port as specified)
pub const VENTURI_INLET_DIAM_M: f64 = 4.0e-3; // 4 mm

/// Serpentine segment lengths to sweep [m]  (one pass = full treatment width; 45 mm)
pub const SERPENTINE_SEG_LENGTHS_M: [f64; 2] = [
    TREATMENT_WIDTH_MM * 1e-3, // 0.045 m — one row sweep
    TREATMENT_WIDTH_MM * 1e-3, // same length; pairs with different segment counts
];

/// Serpentine segment counts to sweep (6 = single pass; 12 = double; 18 = triple — maximises dwell)
pub const SERPENTINE_SEGMENT_COUNTS: [usize; 3] = [6, 12, 18];

/// Fixed bend radius for serpentine turns [m]  (half well pitch)
pub const SERPENTINE_BEND_RADIUS_M: f64 = WELL_PITCH_MM * 0.5e-3; // 4.5 mm

/// Minimum straight-line serpentine segment length to cover ALL 36 treatment wells [m]
/// = `TREATMENT_ROWS` × `segment_length` = 6 × 45 mm = 270 mm
pub const FULL_GRID_SERPENTINE_LENGTH_M: f64 = 6.0 * TREATMENT_WIDTH_MM * 1e-3; // 0.27 m

// ── Leukapheresis micro-scale parameters ─────────────────────────────────────

/// Micro-scale channel widths for leukapheresis topologies [m]  (100–400 µm).
///
/// At 400 µm wide × 60 µm tall: `D_h` = 96 µm → `κ_WBC` = 12/96 = 0.125 > 0.07 ✓
/// At 200 µm wide × 60 µm tall: `D_h` = 92 µm → `κ_WBC` = 0.130 > 0.07 ✓
/// At 100 µm wide × 60 µm tall: `D_h` = 75 µm → `κ_WBC` = 0.160 > 0.07 ✓  (strong focusing)
pub const LEUKA_CHANNEL_WIDTHS_M: [f64; 3] = [100e-6, 200e-6, 400e-6];

/// Micro-scale channel height for leukapheresis [m]  (60 µm — matches Nivedita 2017, Wu 2019).
pub const LEUKA_CHANNEL_HEIGHT_M: f64 = 60e-6;

/// Leukapheresis flow rates per single chip [m³/s]  (150 µL/min – 1 mL/min).
///
/// Six chips in parallel at top rate = 9 mL/min; 10 chips = 15 mL/min — both
/// exceed the 10 mL/min clinical leukapheresis target.
pub const LEUKA_FLOW_RATES_M3_S: [f64; 3] = [2.5e-9, 5.0e-9, 1.5e-8];

// ── Trifurcation width-fraction sweep ────────────────────────────────────────

/// Center-arm width fraction for asymmetric trifurcation topologies.
///
/// Fraction of the parent channel width allocated to the center arm at each
/// trifurcation split; each peripheral arm receives `(1.0 − center_frac) / 2.0`
/// of the parent width.
///
/// - `0.250`: aggressive peripheral bias (strong RBC wall-skimming).
/// - `0.333` (= 1/3): symmetric split — all three arms equal width.
/// - `0.450 / 0.550`: center-biased — larger flow fraction to center arm,
///   stronger Zweifach-Fung bias toward large/stiff cells (cancer, WBC).
/// - `0.650`: highly center-biased — maximum cancer/WBC enrichment at venturi.
pub const TRIFURCATION_CENTER_FRACS: [f64; 5] = [0.250, 0.333, 0.450, 0.550, 0.650];

// ── Asymmetric bifurcation arm width ratio sweep ──────────────────────────

/// Narrow-arm width as a fraction of the wide arm for `AsymmetricBifurcationSerpentine`.
///
/// - `0.25`: 4:1 ratio — extreme enrichment; RBC bypass very low resistance
/// - `0.40`: 2.5:1 ratio — strong Zweifach-Fung bias
/// - `0.50`: 2:1 ratio (original default)
/// - `0.65`: 1.54:1 ratio — moderate bias; lower haemolysis risk in narrow arm
pub const BIFURCATION_ARM_RATIOS: [f64; 4] = [0.25, 0.40, 0.50, 0.65];

/// Independent left-arm width fraction for `AsymmetricTrifurcationVenturi`.
///
/// The right arm fraction is derived: `right_frac = 1 − center_frac − left_frac`.
/// Constraint: `center_frac + left_frac ≤ 0.85` to ensure right arm width ≥ 0.15.
pub const TRIFURCATION_LEFT_FRACS: [f64; 4] = [0.10, 0.20, 0.33, 0.45];

// ── Venturi throat geometry sweep ─────────────────────────────────────────

/// Throat length as a multiple of throat width.
///
/// Longer throats extend residence time in the cavitation zone:
/// - `2.0 × throat_width`: compact throat (original default)
/// - `5.0 × throat_width`: moderate SDT dose extension
/// - `10.0 × throat_width`: long throat — maximises cavitation dose at cost of pressure drop
pub const THROAT_LENGTH_FACTORS: [f64; 3] = [2.0, 5.0, 10.0];

// ── Channel height (aspect ratio) sweep ───────────────────────────────────

/// Channel height [m] for millifluidic topologies.
///
/// Aspect ratio `w/h` governs Dean-flow vortex strength and inertial focusing:
/// - `0.5 mm`: high aspect ratio (w/h up to 16) — strong shear gradients
/// - `1.0 mm`: nominal (original default)
/// - `2.0 mm`: low aspect ratio (w/h down to 0.5) — gentle shear, low PAI
pub const CHANNEL_HEIGHTS_M: [f64; 3] = [0.5e-3, 1.0e-3, 2.0e-3];

// ── Coagulation / platelet activation ─────────────────────────────────────

/// Maximum acceptable platelet activation index (PAI) per single device pass.
///
/// Empirical upper limit from Hellums (1994) for shear-induced platelet
/// activation in blood-contacting devices.  Designs whose PAI exceeds this
/// threshold are penalised in the scoring functions.
///
/// `PAI = 1.8 × 10⁻⁸ × τ^1.325 × t^0.462`
/// where τ = peak shear stress [Pa] and t = exposure duration [s].
pub const PAI_PASS_LIMIT: f64 = 5e-4;

/// Blood-flow target [mL/min] above which low-flow stasis clotting risk is treated as minimal.
///
/// Literature in extracorporeal blood purification reports reduced circuit-failure
/// risk as blood flow increases from ~200 toward ~250 mL/min under no-anticoagulation
/// protocols; this threshold is used as the "low-flow risk ≈ 0" anchor.
pub const CLOTTING_BFR_LOW_RISK_ML_MIN: f64 = 250.0;

/// Blood-flow level [mL/min] treated as high low-flow stasis risk.
///
/// This is intentionally conservative for millifluidic screening and reflects
/// a low-flow boundary where circuit stasis risk is expected to rise sharply.
pub const CLOTTING_BFR_HIGH_RISK_ML_MIN: f64 = 100.0;

/// Clinical caution threshold [mL/min] for low-flow clotting risk reporting.
///
/// Designs below this flow are flagged in reports as elevated stasis-clotting risk.
pub const CLOTTING_BFR_CAUTION_ML_MIN: f64 = 200.0;

/// Conservative sensitivity threshold [mL/min] for the "10 mL/s" clotting claim.
///
/// This is an exploratory reporting threshold (10 mL/s = 600 mL/min) used to
/// quantify how many candidates would remain if that strict assumption were
/// enforced. It is not used as the default clinical gate in scoring.
pub const CLOTTING_BFR_STRICT_10MLS_ML_MIN: f64 = 600.0;

/// Main-channel shear-rate threshold [1/s] above which low-shear stasis risk is minimal.
pub const CLOTTING_SHEAR_LOW_RISK_INV_S: f64 = 1_000.0;

/// Main-channel shear-rate threshold [1/s] below which low-shear stasis risk is high.
pub const CLOTTING_SHEAR_HIGH_RISK_INV_S: f64 = 100.0;

/// Residence-time threshold [s] below which residence-driven stasis risk is minimal.
pub const CLOTTING_RESIDENCE_LOW_RISK_S: f64 = 0.5;

/// Residence-time threshold [s] above which residence-driven stasis risk is high.
pub const CLOTTING_RESIDENCE_HIGH_RISK_S: f64 = 3.0;

// ── Bubble dynamics / Rayleigh–Plesset sonoluminescence parameters ────────────

/// Equilibrium cavitation bubble radius [m] — typical haematogenous nucleus in blood.
///
/// Used for Rayleigh-Plesset bubble-collapse energy estimates.
pub const R_BUBBLE_EQ_M: f64 = 20.0e-6; // 20 µm

/// Blood temperature [K] at 37 °C — reference for adiabatic collapse temperature.
pub const T_BLOOD_K: f64 = 310.15;

/// Adiabatic heat-capacity ratio γ for cavitation bubble contents (air–vapour mixture).
///
/// Used in the Rayleigh-Plesset adiabatic collapse formula:
/// `T_collapse / T_0 = (P_abs / P_vapor)^((γ−1)/γ)`.
pub const BUBBLE_GAMMA: f64 = 1.40;

/// Polytropic exponent k for steam-dominant bubble collapse (k ≈ 1.25 from Storey & Szeri 2000).
///
/// More conservative than fully adiabatic γ = 1.4; accounts for partial heat transfer
/// from the hot bubble core to the surrounding liquid during collapse.
pub const BUBBLE_POLYTROPIC_K: f64 = 1.25;

/// Reference inlet pressure for sonoluminescence normalisation [Pa] (3 bar gauge + atm).
///
/// At this pressure and full cavitation the sonoluminescence proxy equals 1.0.
pub const SONO_REF_P_ABS_PA: f64 = 401_325.0; // 3 bar gauge + 101 325 Pa

// ── 405 nm optical-delivery proxy (Lumidox II) ──────────────────────────────

/// Lumidox II illumination wavelength [nm] used in the optical-delivery proxy.
pub const LUMIDOX_BLUE_WAVELENGTH_NM: f64 = 405.0;

/// Effective attenuation coefficient for 405 nm propagation through blood [1/m].
///
/// Used in a Beer-Lambert proxy:
/// `I/I0 = exp(-μ_eff · path_length · hematocrit_scale)`.
pub const BLOOD_ATTENUATION_405NM_INV_M: f64 = 4_000.0;

// ── Cancer-targeting / inertial-focusing quality thresholds ──────────────────

/// Minimum confinement ratio κ = `a_cell` / `D_h` for effective inertial focusing.
///
/// Di Carlo (2009): cells focus deterministically when κ > 0.07.
pub const KAPPA_FOCUS_MIN: f64 = 0.07;

/// Reference venturi constriction velocity ratio for cavitation-intensity normalisation.
///
/// 80× corresponds to a 50 µm throat in a 4 mm channel (area ratio ≈ 6400 : 1,
/// velocity ratio ≈ 80 for rectangular cross-sections of equal height).
pub const VENTURI_VEL_RATIO_REF: f64 = 80.0;

/// Vena contracta coefficient for millifluidic Venturi throats.
///
/// The flow contracts beyond the geometric throat area; the effective jet
/// cross-section is `A_eff = Cc × A_geom` and the effective velocity is
/// `v_eff = v_geom / Cc`.  For sharp-edged orifices Cc ≈ 0.61 (Borda 1766),
/// for rounded contractions Cc → 0.95–1.0.  In millifluidic Venturi designs
/// the contraction is typically "smooth but short", giving Cc ≈ 0.85
/// (Idelchik, Handbook of Hydraulic Resistance, 7th ed., Diagram 4-9).
pub const VENTURI_CC: f64 = 0.85;

/// Diffuser discharge coefficient for gradual expansion downstream of venturi.
///
/// Models the pressure recovery when the flow expands from the venturi throat
/// back to the main channel cross-section.  An ideal diffuser recovers 100%
/// of the dynamic-pressure difference; real millifluidic diffusers achieve
/// ≈ 80% due to viscous losses and mild separation (Idelchik 1994,
/// Handbook of Hydraulic Resistance, Diagram 6-21).
///
/// `ΔP_recovery = C_D × ½ρ(v_throat² − v_outlet²)`
pub const DIFFUSER_DISCHARGE_COEFF: f64 = 0.80;

/// Rayleigh collapse time scale factor.
///
/// The inertial growth / collapse time for a vapour bubble of radius R₀
/// in liquid of density ρ at driving pressure ΔP = P_∞ − P_v is:
///
/// ```text
/// t_Rayleigh = 0.915 · R₀ · √(ρ / ΔP)
/// ```
///
/// (Lord Rayleigh 1917; Brennen 1995 _Cavitation and Bubble Dynamics_ §2.3).
pub const RAYLEIGH_COLLAPSE_FACTOR: f64 = 0.915;
