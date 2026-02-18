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

// ── Cavitation criterion ───────────────────────────────────────────────────

/// Standard atmospheric pressure [Pa]
pub const P_ATM_PA: f64 = 101_325.0;
/// Cavitation inception threshold (σ < σ_crit → cavitation)
pub const SIGMA_CRIT: f64 = 1.0;

// ── Giersiepen (1990) haemolysis model coefficients ───────────────────────
//   ΔHb/Hb = C_G · t^α · τ^β
//   Reference: Giersiepen et al., Int. J. Artif. Organs 13 (1990) 300–306

/// Pre-multiplied constant C_G
pub const GIERSIEPEN_C: f64 = 3.62e-5;
/// Time exponent α
pub const GIERSIEPEN_ALPHA: f64 = 0.765;
/// Shear stress exponent β
pub const GIERSIEPEN_BETA: f64 = 1.991;

// ── Design parameter sweep ranges ─────────────────────────────────────────

/// Volumetric flow rates to sweep [m³/s]  (1, 5, 10 mL/min)
pub const FLOW_RATES_M3_S: [f64; 3] = [1.667e-8, 8.333e-8, 1.667e-7];

/// Gauge pressures at inlet [Pa]  (1, 2, 3 bar above atmospheric)
pub const INLET_GAUGES_PA: [f64; 3] = [100_000.0, 200_000.0, 300_000.0];

/// Venturi throat diameters [m]  (50, 100, 150 μm)
pub const THROAT_DIAMETERS_M: [f64; 3] = [50e-6, 100e-6, 150e-6];

/// Main channel widths for rectangular serpentine channels [m]  (0.3, 0.5, 0.8 mm)
pub const CHANNEL_WIDTHS_M: [f64; 3] = [0.3e-3, 0.5e-3, 0.8e-3];

/// Fixed channel height used for all rectangular channels [m]
pub const CHANNEL_HEIGHT_M: f64 = 0.2e-3; // 200 μm

/// Fixed inlet diameter for all venturi stages [m]
pub const VENTURI_INLET_DIAM_M: f64 = 0.5e-3; // 500 μm

/// Serpentine segment lengths to sweep [m]  (one pass = full treatment width; 45 mm)
pub const SERPENTINE_SEG_LENGTHS_M: [f64; 2] = [
    TREATMENT_WIDTH_MM * 1e-3,           // 0.045 m — one row sweep
    TREATMENT_WIDTH_MM * 1e-3,           // same length; pairs with different segment counts
];

/// Serpentine segment counts to sweep (6 = single pass over all rows; 12 = double)
pub const SERPENTINE_SEGMENT_COUNTS: [usize; 2] = [6, 12];

/// Fixed bend radius for serpentine turns [m]  (half well pitch)
pub const SERPENTINE_BEND_RADIUS_M: f64 = WELL_PITCH_MM * 0.5e-3; // 4.5 mm

/// Minimum straight-line serpentine segment length to cover ALL 36 treatment wells [m]
/// = TREATMENT_ROWS × segment_length = 6 × 45 mm = 270 mm
pub const FULL_GRID_SERPENTINE_LENGTH_M: f64 = 6.0 * TREATMENT_WIDTH_MM * 1e-3; // 0.27 m
