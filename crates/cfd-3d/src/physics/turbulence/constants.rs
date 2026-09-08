//! Turbulence model constants for k-ω SST and Spalart-Allmaras models.
//!
//! All values are from primary literature:
//! - Menter (1994) for k-ω SST constants
//! - Spalart & Allmaras (1992) / Allmaras et al. (2012) for SA constants

// ── k-ω SST (Menter 1994) ────────────────────────────────────────────────────

/// SST eddy-viscosity limiter constant a₁ = 0.31.
pub const SST_A1: f64 = 0.31;

/// SST β* = 0.09 (k-equation destruction constant, both zones).
pub const SST_BETA_STAR: f64 = 0.09;

// ── Spalart-Allmaras (Spalart & Allmaras 1992; Allmaras et al. 2012) ─────────

/// SA wall-damping constant C_v1 = 7.1.
pub const SA_CV1: f64 = 7.1;

/// SA destruction shape constant C_w2 = 0.3.
pub const SA_CW2: f64 = 0.3;

/// SA destruction shape constant C_w3 = 2.0.
pub const SA_CW3: f64 = 2.0;

// ── DES (Detached Eddy Simulation) ───────────────────────────────────────────

/// DES model constant C_DES = 0.65 (Spalart et al. 1997).
///
/// Controls the transition between RANS and LES regions.  The DES length
/// scale is l_DES = C_DES · Δ where Δ is the local grid spacing.
pub const DES_C_DES: f64 = 0.65;

// ── Smagorinsky LES ─────────────────────────────────────────────────────────

/// Default Smagorinsky constant C_S = 0.1 (Smagorinsky 1963).
///
/// The standard value for channel and free-shear flows.  Range in practice:
/// 0.065 (channel) to 0.18 (isotropic turbulence).
pub const SMAGORINSKY_CS_DEFAULT: f64 = 0.1;

// ── Anisotropic minimum dissipation (AMD) ───────────────────────────────────

/// AMD coefficient C_A = 1/3 for second-order central differencing.
///
/// Rozema et al. (2015) show that the continuous derivation yields 1/12, and
/// the second-order discrete correction for centered finite differences is 1/3.
pub const AMD_C_A_SECOND_ORDER: f64 = 1.0 / 3.0;

/// AMD coefficient C_A = 1/12 for the continuous gradient-model derivation.
pub const AMD_C_A_CONTINUOUS: f64 = 1.0 / 12.0;

// ── Deardorff sub-grid scale ─────────────────────────────────────────────────

/// Exponent 1/3 used in the Deardorff (1980) sub-grid eddy viscosity:
/// ν_sgs = C_k · Δ · k_sgs^{1/3}.
pub const DEARDORFF_ONE_THIRD: f64 = 1.0 / 3.0;

// ── Sigma model (Nicoud et al. 2011) ─────────────────────────────────────────

/// Sigma model constant C_σ ≈ 1.35 (Nicoud et al. 2011).
///
/// Provides the correct asymptotic behaviour near walls without
/// requiring a damping function.
pub const SIGMA_C: f64 = 1.35;

// ── Vreman model (Vreman 2004) ───────────────────────────────────────────────

/// Vreman model constant C_v ≈ 0.07 (Vreman 2004).
///
/// The Vreman model uses this constant in ν_sgs = C_v · √(Bβ / αij αij).
/// This value is derived from C_s = 0.1: C_v ≈ 2.5 · C_s² ≈ 0.025;
/// in practice a calibrated value of 0.07 is common.
pub const VREMAN_CV: f64 = 0.07;

// ── WALE model (Nicoud & Ducros 1999) ───────────────────────────────────────

/// WALE model constant C_w = 0.325 for wall-bounded LES.
///
/// This is the standard calibration used in many wall-resolved LES codes.
pub const WALE_CW: f64 = 0.325;

// ── Wall functions (Spalding 1961; Launder & Spalding 1974) ──────────────────

/// Von Kármán constant κ = 0.41 for the log-law of the wall.
///
/// u⁺ = (1/κ) ln(y⁺) + B  for y⁺ > y⁺_transition.
pub const WALL_KAPPA: f64 = 0.41;

/// Additive constant B ≈ 5.2 in the log-law of the wall.
///
/// Standard value for smooth walls (Coles 1956).
pub const WALL_B: f64 = 5.2;

/// Transition y⁺ ≈ 11.06 between viscous sub-layer and log-law region.
///
/// Below this value the linear relation u⁺ = y⁺ holds; above it the
/// logarithmic profile applies.  The exact value satisfies
/// y⁺ = (1/κ) ln(y⁺) + B.
pub const WALL_Y_PLUS_TRANSITION: f64 = 11.06;
