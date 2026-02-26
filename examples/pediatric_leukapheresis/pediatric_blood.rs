//! Paediatric and neonatal blood parameters for leukapheresis modelling.
//!
//! All values are from clinical literature; references are given inline.
//! Every constant uses SI units (kg, m, s, Pa) unless a suffix indicates otherwise.

// ── Patient volumes ───────────────────────────────────────────────────────────

/// Total blood volume per kg body weight for neonates (Linderkamp 1977).
/// Units: mL / kg.
pub const TBV_PER_KG_ML: f64 = 85.0;

/// Maximum allowable extracorporeal volume as a fraction of total blood volume.
/// Paediatric haematology guideline: ECV ≤ 10 % TBV.
pub const MAX_ECV_FRACTION: f64 = 0.10;

/// Clinical target total throughput for leukapheresis. Units: mL / min.
pub const CLINICAL_FLOW_ML_MIN: f64 = 10.0;

/// Reference patient: 3 kg neonate.
pub const NEONATE_3KG_WEIGHT_KG: f64 = 3.0;

/// Total blood volume of a 3 kg neonate. Units: mL.
pub const NEONATE_3KG_TBV_ML: f64 = NEONATE_3KG_WEIGHT_KG * TBV_PER_KG_ML; // 255 mL

/// ECV budget for the 3 kg neonate patient. Units: mL.
pub const NEONATE_3KG_ECV_BUDGET_ML: f64 = NEONATE_3KG_TBV_ML * MAX_ECV_FRACTION; // 25.5 mL

// ── Blood properties ──────────────────────────────────────────────────────────

/// Diluted blood dynamic viscosity at HCT ≈ 4 % and 37 °C (Chien 1966). Units: Pa·s.
pub const DILUTED_BLOOD_VISCOSITY_PAS: f64 = 1.8e-3;

/// Typical feed hematocrit for post-dilution leukapheresis (1:4 dilution ratio).
/// Dimensionless.
pub const DILUTED_HCT: f64 = 0.04;

/// Whole blood haematocrit (normal adult / paediatric). Dimensionless.
pub const WHOLE_BLOOD_HCT: f64 = 0.45;

// ── Cell diameters ────────────────────────────────────────────────────────────

/// Neonatal erythrocyte (RBC) mean diameter — slightly larger than adult (Bertles 1970).
/// Units: m.
pub const NEONATAL_RBC_DIAMETER_M: f64 = 8.5e-6;

/// Average white blood cell (mixed WBC population) diameter.  Units: m.
pub const WBC_DIAMETER_M: f64 = 12.0e-6;

// ── Paper-reported reference performance (PMC10645346 review) ─────────────────

/// WBC recovery reported by Nivedita et al. (2017) — inertial spiral design.
pub const NIVEDITA_WBC_RECOVERY: f64 = 0.950;

/// RBC removal reported by Nivedita et al. (2017).
/// Fraction of RBCs directed to the peripheral outlet.
pub const NIVEDITA_RBC_REMOVAL: f64 = 0.940;

/// WBC separation efficiency reported by Wu Z. et al. (2019) — microstructure constriction.
pub const WU_WBC_RECOVERY: f64 = 0.897;

/// WBC purity in the centre outlet reported by Wu Z. et al. (2019).
pub const WU_WBC_PURITY: f64 = 0.910;
