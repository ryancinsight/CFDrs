//! Haemolysis index models for millifluidic and microfluidic flow.
//!
//! Provides the Giersiepen (1990) power-law haemolysis index model and a
//! conservative cavitation-amplification correction for SDT (Sonodynamic
//! Therapy) millifluidic devices.
//!
//! # Giersiepen (1990) Model
//!
//! The empirical power-law model relates haemoglobin release to shear stress
//! and exposure duration:
//!
//! ```text
//!   HI = C · t^α · τ^β
//! ```
//!
//! where:
//! - `C  = 3.62 × 10⁻⁵`  (fit constant)
//! - `α  = 0.765`         (time exponent)
//! - `β  = 1.991`         (shear exponent)
//! - `t` = exposure duration \[s]
//! - `τ` = wall shear stress \[Pa]
//!
//! Reference: Giersiepen M. et al. (1990) *Estimation of shear stress-related
//! blood damage in heart valve prostheses*. Int. J. Artif. Organs 13(5):300–306.
//!
//! # Cavitation Amplification
//!
//! In SDT millifluidic devices, acoustic bubble collapse generates micro-jets
//! and shockwaves that cause RBC membrane damage independently of the
//! macroscopic shear stress.  A conservative 3× amplification factor at full
//! cavitation potential is applied:
//!
//! ```text
//!   HI_amplified = HI_base × (1 + 3 × cav_potential)
//! ```

pub mod acoustic_radiation;
mod dynamics;
mod models;

// ── Giersiepen model constants (re-exported from cfd-core SSOT) ───────────────
//
// These are the single source of truth from cfd-core, re-exported here under
// their canonical names.  All three fidelity levels (1D/2D/3D) share identical
// constants, ensuring cross-fidelity HI comparisons are valid.

/// Giersiepen (1990) fit constant C — from cfd-core SSOT.
///
/// See [`cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_C`].
pub use cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_C;

/// Giersiepen (1990) time exponent α — from cfd-core SSOT.
///
/// See [`cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_TIME`].
pub use cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_TIME;

/// Giersiepen (1990) shear exponent β — from cfd-core SSOT.
///
/// See [`cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_STRESS`].
pub use cfd_core::physics::hemolysis::GIERSIEPEN_MILLIFLUIDIC_STRESS;

/// Conservative cavitation amplification slope — from cfd-core SSOT.
///
/// See [`cfd_core::physics::hemolysis::CAVITATION_HI_SLOPE`].
pub use cfd_core::physics::hemolysis::CAVITATION_HI_SLOPE;

// ── Public API ────────────────────────────────────────────────────────────────
pub use dynamics::{
    HemolysisExposure, P_REF_ATMOSPHERIC, RAYLEIGH_ALPHA, cavitation_hemolysis_amplification,
    collapse_jet_velocity, rayleigh_collapse_time,
};
pub use models::{
    SENSITIZER_K_ACT_CHLORIN_E6, SENSITIZER_K_ACT_HEMATOPORPHYRIN, TASKIN_BETA, TASKIN_C,
    cavitation_amplified_hi, giersiepen_hi, sonosensitizer_activation_efficiency, taskin_hi,
};

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests;
