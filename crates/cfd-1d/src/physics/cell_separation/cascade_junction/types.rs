//! Shared descriptors for the cascade junction model.
//!
//! These types are the module's data contract: they carry the per-stage arm
//! flow fractions, the treatment-arm hydraulic diameter, and the peripheral
//! recovery sub-splits that the routing leaves consume, plus the result
//! structs those leaves produce.

use aequitas::systems::si::quantities::{Length, Velocity};

/// Peripheral recovery sub-split descriptor.
///
/// When a non-treatment arm is further sub-split, the wider sub-arm can
/// feed recovered cells back to the treatment path.  This struct
/// parameterises that recovery routing for a single source arm.
#[derive(Debug, Clone, Copy)]
pub struct PeripheralRecovery {
    /// Index of the source arm in the parent stage's `arm_q_fracs`.
    pub source_arm_idx: usize,
    /// Sub-arm flow fractions (up to 5); only `[..n_sub_arms]` are used.
    pub sub_arm_q_fracs: [f64; 5],
    /// Number of active sub-arms (2–5).
    pub n_sub_arms: u8,
    /// Index of the sub-arm that feeds back to the treatment path.
    pub recovery_arm_idx: usize,
    /// Hydraulic diameter of the recovery sub-arm.
    pub recovery_dh_m: Length,
}

/// Per-stage descriptor for a mixed Bi/Tri selective-routing cascade.
///
/// Carries all arm flow fractions (enabling asymmetric splits) and the
/// treatment-arm hydraulic diameter (enabling κ-dependent β correction).
///
/// Index 0 of `arm_q_fracs` is always the **treatment arm** (center at
/// trifurcations, treatment branch at bifurcations).
///
/// For a bifurcation, `arm_q_fracs = [q_treat, q_bypass, 0.0]` and `n_arms = 2`.
/// For a trifurcation, `arm_q_fracs = [q_center, q_left, q_right]` and `n_arms = 3`.
#[derive(Debug, Clone, Copy)]
pub struct CascadeStage {
    /// Flow fractions for all arms; first entry is the treatment arm.
    pub arm_q_fracs: [f64; 5],
    /// Number of active arms (2 = bifurcation, 3 = trifurcation, 4 = quad, 5 = penta).
    pub n_arms: u8,
    /// Hydraulic diameter of the treatment arm at this stage.
    /// Used to compute κ = cell_diameter / Dh for β amplification.
    pub treatment_dh_m: Length,
    /// Inflow velocity into the junction.
    /// Used for PMC5114676 Zweifach-Fung high-velocity inversion mechanics.
    pub parent_v_in_m_s: Velocity,
    /// Optional peripheral recovery sub-splits (up to 4 per stage).
    pub peripheral_recoveries: [Option<PeripheralRecovery>; 4],
    /// Number of active peripheral recoveries.
    pub n_recoveries: u8,
}

/// Cell fractions at the center (venturi) arm and peripheral bypass arms after
/// N cascade trifurcation levels.
#[derive(Debug, Clone, Copy)]
pub struct CascadeJunctionResult {
    /// Fraction of input cancer cells that reach the deepest center arm (venturi).
    pub cancer_center_fraction: f64,
    /// Fraction of input WBCs that reach the deepest center arm.
    pub wbc_center_fraction: f64,
    /// Fraction of input RBCs that are diverted to peripheral bypass arms
    /// at any cascade level (`= 1 − rbc_center_fraction`).
    pub rbc_peripheral_fraction: f64,
    /// Separation efficiency = `|f_cancer_center − f_rbc_center|` ∈ [0, 1].
    ///
    /// High values indicate strong enrichment of cancer cells at the venturi
    /// and strong depletion of RBCs (which are protected in bypass channels).
    pub separation_efficiency: f64,
    /// Local hematocrit in the center arm at the deepest level, relative to
    /// the feed hematocrit. Computed from the flow fraction and RBC routing:
    /// `HCT_local = HCT_feed × (rbc_center / q_center_frac^n_levels)`.
    pub center_hematocrit_ratio: f64,
}

/// Cell fractions after a staged selective-routing path:
///
/// 1. `n_pretri` center-only cascade trifurcation stages (incremental skimming),
/// 2. one terminal trifurcation skimming stage,
/// 3. one terminal asymmetric bifurcation selecting the treatment arm.
#[derive(Debug, Clone, Copy)]
pub struct IncrementalFiltrationResult {
    /// Fraction of input cancer cells reaching the treatment venturi arm.
    pub cancer_center_fraction: f64,
    /// Fraction of input WBCs reaching the treatment venturi arm.
    pub wbc_center_fraction: f64,
    /// Fraction of input RBCs that still reach the treatment venturi arm.
    pub rbc_center_fraction: f64,
    /// Fraction of input RBCs skimmed to peripheral bypass streams.
    pub rbc_peripheral_fraction: f64,
    /// Separation efficiency = `|f_cancer_center − f_rbc_center|` ∈ [0, 1].
    pub separation_efficiency: f64,
    /// Local hematocrit ratio (HCT_local / HCT_feed) at the treatment arm.
    pub center_hematocrit_ratio: f64,
}