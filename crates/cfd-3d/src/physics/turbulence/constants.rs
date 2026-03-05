//! Turbulence model constants for k-ω SST and Spalart-Allmaras models.
//!
//! All values are from primary literature:
//! - Menter (1994) for k-ω SST constants
//! - Spalart & Allmaras (1992) / Allmaras et al. (2012) for SA constants

// ── k-ω SST (Menter 1994) ────────────────────────────────────────────────────

/// SST eddy-viscosity limiter constant a₁ = 0.31.
#[allow(dead_code)]
pub const SST_A1: f64 = 0.31;

/// SST k-ω (inner) closure constant α₁ = 5/9.
#[allow(dead_code)]
pub const SST_ALPHA1: f64 = 5.0 / 9.0;

/// SST k-ε (outer) closure constant α₂ = 0.44.
#[allow(dead_code)]
pub const SST_ALPHA2: f64 = 0.44;

/// SST k-ω β₁ = 3/40 (inner destruction constant).
#[allow(dead_code)]
pub const SST_BETA1: f64 = 3.0 / 40.0;

/// SST k-ε β₂ = 0.0828 (outer destruction constant).
#[allow(dead_code)]
pub const SST_BETA2: f64 = 0.0828;

/// SST β* = 0.09 (k-equation destruction constant, both zones).
#[allow(dead_code)]
pub const SST_BETA_STAR: f64 = 0.09;

/// SST turbulent Prandtl number for k (inner zone) σ_k1 = 0.85.
#[allow(dead_code)]
pub const SST_SIGMA_K1: f64 = 0.85;

/// SST turbulent Prandtl number for k (outer zone) σ_k2 = 1.0.
#[allow(dead_code)]
pub const SST_SIGMA_K2: f64 = 1.0;

/// SST turbulent Prandtl number for ω (inner zone) σ_ω1 = 0.5.
#[allow(dead_code)]
pub const SST_SIGMA_OMEGA1: f64 = 0.5;

/// SST turbulent Prandtl number for ω (outer zone) σ_ω2 = 0.856.
#[allow(dead_code)]
pub const SST_SIGMA_OMEGA2: f64 = 0.856;

// ── Spalart-Allmaras (Spalart & Allmaras 1992; Allmaras et al. 2012) ─────────

/// SA production constant C_b1 = 0.1355.
#[allow(dead_code)]
pub const SA_CB1: f64 = 0.1355;

/// SA diffusion constant C_b2 = 0.622.
#[allow(dead_code)]
pub const SA_CB2: f64 = 0.622;

/// SA wall-damping constant C_v1 = 7.1.
#[allow(dead_code)]
pub const SA_CV1: f64 = 7.1;

/// SA destruction constant C_w1 = C_b1/κ² + (1 + C_b2)/σ.
/// Pre-computed: 0.1355/0.1681 + 1.622/0.6667 ≈ 3.239.
#[allow(dead_code)]
pub const SA_CW1: f64 = 3.239_067_817_268_056;

/// SA destruction shape constant C_w2 = 0.3.
#[allow(dead_code)]
pub const SA_CW2: f64 = 0.3;

/// SA destruction shape constant C_w3 = 2.0.
#[allow(dead_code)]
pub const SA_CW3: f64 = 2.0;

/// SA von Kármán constant κ = 0.41.
#[allow(dead_code)]
pub const SA_KAPPA: f64 = 0.41;

/// SA diffusion coefficient σ = 2/3.
#[allow(dead_code)]
pub const SA_SIGMA: f64 = 2.0 / 3.0;

// ── DES (Detached Eddy Simulation) ───────────────────────────────────────────

/// DES model constant C_DES = 0.65 (Spalart et al. 1997).
///
/// Controls the transition between RANS and LES regions.  The DES length
/// scale is l_DES = C_DES · Δ where Δ is the local grid spacing.
#[allow(dead_code)]
pub const DES_C_DES: f64 = 0.65;

// ── Smagorinsky LES ─────────────────────────────────────────────────────────

/// Default Smagorinsky constant C_S = 0.1 (Smagorinsky 1963).
///
/// The standard value for channel and free-shear flows.  Range in practice:
/// 0.065 (channel) to 0.18 (isotropic turbulence).
#[allow(dead_code)]
pub const SMAGORINSKY_CS_DEFAULT: f64 = 0.1;

// ── Deardorff sub-grid scale ─────────────────────────────────────────────────

/// Exponent 1/3 used in the Deardorff (1980) sub-grid eddy viscosity:
/// ν_sgs = C_k · Δ · k_sgs^{1/3}.
#[allow(dead_code)]
pub const DEARDORFF_ONE_THIRD: f64 = 1.0 / 3.0;

// ── Sigma model (Nicoud et al. 2011) ─────────────────────────────────────────

/// Sigma model constant C_σ ≈ 1.35 (Nicoud et al. 2011).
///
/// Provides the correct asymptotic behaviour near walls without
/// requiring a damping function.
#[allow(dead_code)]
pub const SIGMA_C: f64 = 1.35;

// ── Vreman model (Vreman 2004) ───────────────────────────────────────────────

/// Vreman model constant C_v ≈ 0.07 (Vreman 2004).
///
/// The Vreman model uses this constant in ν_sgs = C_v · √(Bβ / αij αij).
/// This value is derived from C_s = 0.1: C_v ≈ 2.5 · C_s² ≈ 0.025;
/// in practice a calibrated value of 0.07 is common.
#[allow(dead_code)]
pub const VREMAN_CV: f64 = 0.07;

// ── Wall functions (Spalding 1961; Launder & Spalding 1974) ──────────────────

/// Von Kármán constant κ = 0.41 for the log-law of the wall.
///
/// u⁺ = (1/κ) ln(y⁺) + B  for y⁺ > y⁺_transition.
#[allow(dead_code)]
pub const WALL_KAPPA: f64 = 0.41;

/// Additive constant B ≈ 5.2 in the log-law of the wall.
///
/// Standard value for smooth walls (Coles 1956).
#[allow(dead_code)]
pub const WALL_B: f64 = 5.2;

/// Transition y⁺ ≈ 11.06 between viscous sub-layer and log-law region.
///
/// Below this value the linear relation u⁺ = y⁺ holds; above it the
/// logarithmic profile applies.  The exact value satisfies
/// y⁺ = (1/κ) ln(y⁺) + B.
#[allow(dead_code)]
pub const WALL_Y_PLUS_TRANSITION: f64 = 11.06;
