//! High-level cell separation model for millifluidic device design.
//!
//! [`CellSeparationModel`] combines the inertial lift and Dean drag force
//! models from [`crate::cell_separation::margination`] to predict the
//! **separation efficiency** of a rectangular microchannel for two cell
//! populations: a target (cancer) cell and a background (healthy) cell.
//!
//! # Separation efficiency
//!
//! Separation efficiency `η` is defined as the fraction of the channel
//! cross-section that separates the two cell populations:
//!
//! ```text
//! η = |x̃_cancer − x̃_healthy| / 1.0
//! ```
//!
//! where `x̃` is the normalised equilibrium lateral position (0 = center,
//! 1 = wall).  `η = 1` means the two populations are at opposite extremes
//! of the channel; `η = 0` means they co-focus at the same position.
//!
//! # Center-channel fraction
//!
//! For a venturi-based cell separation device, the center channel collects
//! cells within `|x̃| < x̃_split` where `x̃_split` is the split position
//! (typically 0.3 for a center channel that is 60% of the channel width).
//!
//! The **cancer cell center fraction** is the probability that a cancer cell
//! is within the center channel:
//!
//! ```text
//! P_cancer_center = 1 − x̃_cancer / 1.0   (if x̃_cancer < x̃_split)
//!                 = 0                       (if x̃_cancer ≥ x̃_split)
//! ```
//!
//! The **healthy cell peripheral fraction** is the probability that a healthy
//! cell is in the peripheral (wall) channels:
//!
//! ```text
//! P_healthy_wall = x̃_healthy / 1.0   (if x̃_healthy > x̃_split)
//!               = 0                   (if x̃_healthy ≤ x̃_split)
//! ```
//!
//! # References
//! - Di Carlo, D. (2009). Inertial microfluidics. *Lab Chip*, 9, 3038–3046.
//! - Gossett, D. R. & Di Carlo, D. (2009). *Anal. Chem.*, 81, 8459–8465.
//! - Hur, S. C. et al. (2011). *Lab Chip*, 11, 912–920.

use crate::cell_separation::margination::{lateral_equilibrium, EquilibriumResult};
use crate::cell_separation::properties::CellProperties;
use serde::{Deserialize, Serialize};

// ── Analysis result ───────────────────────────────────────────────────────────

/// Full separation analysis for a pair of cell populations in a microchannel.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct CellSeparationAnalysis {
    /// Equilibrium position of the target (cancer) cell population.
    pub target_equilibrium: EquilibriumResult,

    /// Equilibrium position of the background (healthy) cell population.
    pub background_equilibrium: EquilibriumResult,

    /// Separation efficiency `η = |x̃_target − x̃_background|` ∈ [0, 1].
    ///
    /// 0 = no separation; 1 = maximum separation (opposite walls).
    pub separation_efficiency: f64,

    /// Fraction of target (cancer) cells collected in the center channel.
    ///
    /// Computed for a center channel spanning `|x̃| < x̃_split`.
    pub target_center_fraction: f64,

    /// Fraction of background (healthy) cells in the peripheral channels.
    ///
    /// Computed for peripheral channels spanning `|x̃| > x̃_split`.
    pub background_peripheral_fraction: f64,

    /// Combined purity metric: geometric mean of center fraction and
    /// peripheral fraction.  `purity = √(target_center × background_peripheral)`.
    ///
    /// Ranges from 0 (no separation) to 1 (perfect separation).
    pub purity: f64,

    /// Channel split position `x̃_split` used for fraction calculations.
    pub x_tilde_split: f64,
}

// ── Model ─────────────────────────────────────────────────────────────────────

/// Cell separation model for a rectangular microchannel.
///
/// Predicts the lateral equilibrium positions of two cell populations and
/// computes separation efficiency and center-channel collection fractions.
///
/// # Usage
///
/// ```rust,no_run
/// use cfd_1d::cell_separation::{CellProperties, CellSeparationModel};
///
/// let cancer = CellProperties::mcf7_breast_cancer();
/// let healthy = CellProperties::red_blood_cell();
///
/// let model = CellSeparationModel::new(
///     500e-6,  // channel width [m]
///     200e-6,  // channel height [m]
///     None,    // straight channel (no curvature)
/// );
///
/// let analysis = model.analyze(
///     &cancer,
///     &healthy,
///     1060.0,  // blood density [kg/m³]
///     3.5e-3,  // blood viscosity [Pa·s]
///     0.05,    // mean velocity [m/s]
/// ).expect("cell focusing requires κ > 0.07");
///
/// println!("Separation efficiency: {:.2}", analysis.separation_efficiency);
/// ```
#[derive(Debug, Clone)]
pub struct CellSeparationModel {
    /// Channel width [m].
    pub channel_width_m: f64,
    /// Channel height [m] (shorter dimension for rectangular channels).
    pub channel_height_m: f64,
    /// Radius of curvature [m] for curved channels, or `None` for straight.
    pub bend_radius_m: Option<f64>,
    /// Center-channel split position `x̃_split` ∈ (0, 1).
    ///
    /// Cells with `x̃ < x̃_split` are collected in the center channel.
    /// Default: 0.3 (center channel spans 60% of channel half-width).
    pub x_tilde_split: f64,
}

impl CellSeparationModel {
    /// Construct a new model for a rectangular channel.
    ///
    /// # Arguments
    /// - `channel_width_m` — channel width [m]
    /// - `channel_height_m` — channel height [m] (shorter dimension)
    /// - `bend_radius_m` — radius of curvature [m], or `None` for straight
    #[must_use]
    pub fn new(
        channel_width_m: f64,
        channel_height_m: f64,
        bend_radius_m: Option<f64>,
    ) -> Self {
        Self {
            channel_width_m,
            channel_height_m,
            bend_radius_m,
            x_tilde_split: 0.3,
        }
    }

    /// Set the center-channel split position `x̃_split`.
    ///
    /// Cells with `x̃ < x̃_split` are collected in the center channel.
    /// Must be in (0, 1).
    #[must_use]
    pub fn with_split(mut self, x_tilde_split: f64) -> Self {
        self.x_tilde_split = x_tilde_split.clamp(0.05, 0.95);
        self
    }

    /// Hydraulic diameter `D_h = 2wh / (w + h)` [m].
    #[inline]
    #[must_use]
    pub fn hydraulic_diameter_m(&self) -> f64 {
        let w = self.channel_width_m;
        let h = self.channel_height_m;
        2.0 * w * h / (w + h)
    }

    /// Analyse the separation of two cell populations in this channel.
    ///
    /// # Arguments
    /// - `target` — target (cancer) cell properties
    /// - `background` — background (healthy) cell properties
    /// - `fluid_density_kg_m3` — fluid density [kg/m³]
    /// - `dynamic_viscosity_pa_s` — fluid dynamic viscosity [Pa·s]
    /// - `mean_velocity_m_s` — mean channel velocity [m/s]
    ///
    /// # Returns
    /// `Some(CellSeparationAnalysis)` if at least one cell type focuses
    /// (κ > 0.07), or `None` if neither cell type will focus in this channel.
    #[must_use]
    pub fn analyze(
        &self,
        target: &CellProperties,
        background: &CellProperties,
        fluid_density_kg_m3: f64,
        dynamic_viscosity_pa_s: f64,
        mean_velocity_m_s: f64,
    ) -> Option<CellSeparationAnalysis> {
        let make_dispersed = || EquilibriumResult {
            x_tilde_eq: 0.5,
            lateral_position_m: 0.0,
            residual_force_n: 0.0,
            dean_drag_n: 0.0,
            reynolds_number: 0.0, // not computed if skipped, or could duplicate
            dean_number: 0.0,
            will_focus: false,
        };

        let target_eq = lateral_equilibrium(
            target,
            fluid_density_kg_m3,
            dynamic_viscosity_pa_s,
            mean_velocity_m_s,
            self.channel_width_m,
            self.channel_height_m,
            self.bend_radius_m,
        )
        .unwrap_or_else(make_dispersed);

        let background_eq = lateral_equilibrium(
            background,
            fluid_density_kg_m3,
            dynamic_viscosity_pa_s,
            mean_velocity_m_s,
            self.channel_width_m,
            self.channel_height_m,
            self.bend_radius_m,
        )
        .unwrap_or_else(make_dispersed);

        // If neither focuses, return None (device does nothing)
        if !target_eq.will_focus && !background_eq.will_focus {
            return None;
        }

        // Separation efficiency:
        // If both focus: |x_t - x_b|
        // If one focuses and other is dispersed (x=0.5): |x_t - 0.5|
        let sep_eff = (target_eq.x_tilde_eq - background_eq.x_tilde_eq).abs();

        let split = self.x_tilde_split;

        // Target (cancer) center fraction
        let target_center = if target_eq.will_focus {
            if target_eq.x_tilde_eq < split {
                // Focuses near center
                (1.0 - target_eq.x_tilde_eq / split).clamp(0.0, 1.0)
            } else {
                // Focuses near wall
                0.0
            }
        } else {
            // Dispersed: fraction in center channel = split width / total width
            split
        };

        // Background (healthy) peripheral fraction
        let background_peripheral = if background_eq.will_focus {
            if background_eq.x_tilde_eq > split {
                // Focuses near wall
                ((background_eq.x_tilde_eq - split) / (1.0 - split)).clamp(0.0, 1.0)
            } else {
                // Focuses near center
                0.0
            }
        } else {
            // Dispersed: fraction in peripheral channel = (1 - split)
            1.0 - split
        };

        // Purity: geometric mean of target center and background peripheral fractions
        let purity = (target_center * background_peripheral).sqrt();

        Some(CellSeparationAnalysis {
            target_equilibrium: target_eq,
            background_equilibrium: background_eq,
            separation_efficiency: sep_eff,
            target_center_fraction: target_center,
            background_peripheral_fraction: background_peripheral,
            purity,
            x_tilde_split: split,
        })
    }
}
