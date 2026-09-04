//! Pries phase-separation model for bifurcation hematocrit partitioning.
//!
//! ## Theorem — Plasma Skimming (Pries et al. 1989; Pries & Secomb 2005)
//!
//! At a bifurcation, the hematocrit in each daughter branch differs from the
//! feed hematocrit due to the cell-free layer effect:
//!
//! ```text
//! logit(FQ_E) = A + B · logit((FQ_B − X₀) / (1 − 2X₀))
//! ```
//!
//! where:
//! - FQ_B = Q_daughter / Q_total (fractional blood flow to daughter)
//! - FQ_E = (H_daughter · Q_daughter) / (H_feed · Q_total) (fractional RBC flux)
//! - X₀ = minimum flow fraction required for RBCs to enter the daughter
//! - A depends on daughter-diameter asymmetry
//! - B = 1 + 6.98 · (1 − H_feed) / D_feed
//! - logit(x) = ln(x / (1−x))
//!
//! **Physical basis**: The cell-free layer (thickness δ ≈ 1-3 µm) near the
//! wall is depleted of RBCs. At a bifurcation, the daughter branch that
//! receives more wall-adjacent plasma (typically the smaller one) gets
//! a disproportionately low hematocrit.
//!
//! **Reference**: Pries, A.R., Ley, K., Claassen, M. & Gaehtgens, P. (1989).
//! "Red cell distribution at microvascular bifurcations",
//! *Microvasc. Res.* 38:81-101.
//! Also: Pries, A.R. & Secomb, T.W. (2005). "Microvascular blood viscosity
//! in vivo and the endothelial surface layer", *Am. J. Physiol.* 289:H2657-H2664.

use super::fahraeus_lindqvist::secomb_phase_separation_x0;
use aequitas::systems::si::quantities::Length;
use cfd_core::error::{Error, Result};

/// Result of the Pries phase-separation model for one daughter branch.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PhaseSeparationResult {
    /// Fractional blood flow into the daughter branch.
    pub flow_fraction: f64,
    /// Fractional erythrocyte flux into the daughter branch.
    pub cell_fraction: f64,
    /// Hematocrit in the daughter branch.
    pub daughter_hematocrit: f64,
    /// Minimum flow fraction required for cells to enter the daughter branch.
    pub x0: f64,
}

/// Compute the threshold-aware Pries phase-separation model for one daughter.
pub fn pries_phase_separation(
    feed_hematocrit: f64,
    flow_fraction: f64,
    diameter_daughter_alpha: Length,
    diameter_daughter_beta: Length,
    diameter_feed: Length,
) -> Result<PhaseSeparationResult> {
    let checked = checked_pries_phase_separation(
        feed_hematocrit,
        flow_fraction,
        diameter_daughter_alpha,
        diameter_daughter_beta,
        diameter_feed,
    )?;

    Ok(PhaseSeparationResult {
        daughter_hematocrit: checked
            .daughter_hematocrit
            .clamp(0.0, (2.0 * feed_hematocrit).min(1.0)),
        cell_fraction: checked.cell_fraction.clamp(0.0, 1.0),
        ..checked
    })
}

/// Checked threshold-aware Pries phase-separation model for one daughter.
pub fn checked_pries_phase_separation(
    feed_hematocrit: f64,
    flow_fraction: f64,
    diameter_daughter_alpha: Length,
    diameter_daughter_beta: Length,
    diameter_feed: Length,
) -> Result<PhaseSeparationResult> {
    let diameter_daughter_alpha_m = diameter_daughter_alpha.into_base();
    let diameter_daughter_beta_m = diameter_daughter_beta.into_base();
    let diameter_feed_m = diameter_feed.into_base();
    if !feed_hematocrit.is_finite()
        || !flow_fraction.is_finite()
        || !diameter_daughter_alpha_m.is_finite()
        || !diameter_daughter_beta_m.is_finite()
        || !diameter_feed_m.is_finite()
    {
        return Err(Error::InvalidConfiguration(
            "Pries phase-separation inputs must be finite".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&feed_hematocrit) {
        return Err(Error::InvalidConfiguration(
            "Pries phase-separation feed hematocrit must lie in [0, 1]".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&flow_fraction) {
        return Err(Error::InvalidConfiguration(
            "Pries phase-separation flow fraction must lie in [0, 1]".to_string(),
        ));
    }
    if diameter_daughter_alpha_m <= 0.0 || diameter_daughter_beta_m <= 0.0 || diameter_feed_m <= 0.0
    {
        return Err(Error::InvalidConfiguration(
            "Pries phase-separation diameters must be positive".to_string(),
        ));
    }

    let h_feed = feed_hematocrit;
    let d_feed = diameter_feed_m * 1.0e6;
    let d_alpha = diameter_daughter_alpha_m * 1.0e6;
    let d_beta = diameter_daughter_beta_m * 1.0e6;

    if flow_fraction <= 0.0 {
        return Ok(PhaseSeparationResult {
            flow_fraction,
            cell_fraction: 0.0,
            daughter_hematocrit: 0.0,
            x0: 0.0,
        });
    }
    if flow_fraction >= 1.0 {
        return Ok(PhaseSeparationResult {
            flow_fraction,
            cell_fraction: 1.0,
            daughter_hematocrit: h_feed,
            x0: 0.0,
        });
    }

    if h_feed < 1e-15 {
        return Ok(PhaseSeparationResult {
            flow_fraction,
            cell_fraction: 0.0,
            daughter_hematocrit: 0.0,
            x0: 0.0,
        });
    }

    let x0 = secomb_phase_separation_x0(Length::from_base(d_feed * 1.0e-6), h_feed)?;
    if x0 >= 0.5 {
        return Err(Error::InvalidConfiguration(
            "Pries phase-separation x0 must be less than 0.5".to_string(),
        ));
    }

    let diameter_ratio_sq = d_alpha * d_alpha / (d_beta * d_beta);
    let a_param =
        -13.29 * (diameter_ratio_sq - 1.0) / (diameter_ratio_sq + 1.0) * (1.0 - h_feed) / d_feed;
    let b_param = 1.0 + 6.98 * (1.0 - h_feed) / d_feed;

    let cell_fraction = if flow_fraction <= x0 {
        0.0
    } else if flow_fraction >= 1.0 - x0 {
        1.0
    } else {
        let normalized = (flow_fraction - x0) / (1.0 - 2.0 * x0);
        let logit_input = logit(normalized);
        let logit_fqe = a_param + b_param * logit_input;
        inv_logit(logit_fqe)
    };

    if !(0.0..=1.0).contains(&cell_fraction) {
        return Err(Error::InvalidConfiguration(
            "Pries phase-separation produced cell fraction outside [0, 1]".to_string(),
        ));
    }

    let daughter_hematocrit = cell_fraction * h_feed / flow_fraction;
    let daughter_upper_bound = (2.0 * h_feed).min(1.0);
    if !(0.0..=daughter_upper_bound).contains(&daughter_hematocrit) {
        return Err(Error::InvalidConfiguration(
            "Pries phase-separation produced daughter hematocrit outside the physical compact-model bounds".to_string(),
        ));
    }

    Ok(PhaseSeparationResult {
        flow_fraction,
        cell_fraction,
        daughter_hematocrit,
        x0,
    })
}

/// Compute the daughter-branch hematocrit after plasma skimming at a
/// microvascular bifurcation.
///
/// This convenience wrapper keeps the historical compact 1D API and infers the
/// missing sibling branch from Murray cubic diameter conservation before
/// dispatching to [`checked_pries_phase_separation`].
///
/// # Theorem — Compact Pries Projection
///
/// If a parent diameter `D0` and one daughter diameter `D1` are known, Murray's
/// law gives the sibling diameter
///
/// ```text
/// D2 = (max(D0^3 - D1^3, D1^3))^(1/3).
/// ```
///
/// Passing `(D1, D2, D0)` into the Pries phase-separation law preserves the
/// explicit `X0` cell-entry threshold and the bounded erythrocyte-flux logit
/// map, so the compact wrapper remains a projection of the full bifurcation
/// model instead of an independent screening law.
///
/// **Proof.** Murray conservation closes the missing geometric degree of
/// freedom with a positive sibling diameter.  [`checked_pries_phase_separation`]
/// then maps valid flow fractions through the thresholded normalized interval
/// `(FQ_B - X0)/(1 - 2X0)`, assigns zero RBC flux below `X0`, and computes
/// daughter hematocrit from conserved RBC flux `H_i Q_i = FQ_E H_0 Q_0`.
/// Therefore the wrapper inherits boundedness and threshold behavior from the
/// full model. ∎
///
/// # Arguments
/// * `feed_hematocrit` - Feed (parent) hematocrit [0, 1]
/// * `flow_fraction` - Fractional volumetric flow to this daughter branch,
///   Q_daughter / Q_total [0, 1]
/// * `diameter_daughter` - Daughter branch diameter as an SI `Length`.
/// * `diameter_feed` - Feed (parent) branch diameter as an SI `Length`.
///
/// # Returns
/// Daughter branch hematocrit, clamped to [0, min(1, 2 × H_feed)].
#[inline]
pub fn plasma_skimming_hematocrit(
    feed_hematocrit: f64,
    flow_fraction: f64,
    diameter_daughter: Length,
    diameter_feed: Length,
) -> Result<f64> {
    Ok(checked_plasma_skimming_hematocrit(
        feed_hematocrit,
        flow_fraction,
        diameter_daughter,
        diameter_feed,
    )?
    .clamp(0.0, (2.0 * feed_hematocrit).min(1.0)))
}

/// Checked compact plasma-skimming wrapper with explicit output-bound enforcement.
pub fn checked_plasma_skimming_hematocrit(
    feed_hematocrit: f64,
    flow_fraction: f64,
    diameter_daughter: Length,
    diameter_feed: Length,
) -> Result<f64> {
    if !feed_hematocrit.is_finite()
        || !flow_fraction.is_finite()
        || !diameter_daughter.into_base().is_finite()
        || !diameter_feed.into_base().is_finite()
    {
        return Err(Error::InvalidConfiguration(
            "Plasma-skimming inputs must be finite".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&feed_hematocrit) {
        return Err(Error::InvalidConfiguration(
            "Plasma-skimming feed hematocrit must lie in [0, 1]".to_string(),
        ));
    }
    if !(0.0..=1.0).contains(&flow_fraction) {
        return Err(Error::InvalidConfiguration(
            "Plasma-skimming flow fraction must lie in [0, 1]".to_string(),
        ));
    }
    if diameter_daughter.into_base() <= 0.0 || diameter_feed.into_base() <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Plasma-skimming diameters must be positive".to_string(),
        ));
    }

    let ht_feed = feed_hematocrit;
    let fq = flow_fraction;
    let d_daughter = diameter_daughter.into_base() * 1.0e6;
    let d_feed = diameter_feed.into_base() * 1.0e6;

    if ht_feed < 1e-15 || fq < 1e-15 {
        return Ok(0.0);
    }
    if fq >= 1.0 - 1e-15 {
        return Ok(ht_feed);
    }

    let sibling = inferred_murray_sibling_diameter(d_daughter, d_feed);
    Ok(checked_pries_phase_separation(
        ht_feed,
        fq,
        Length::from_base(d_daughter / 1.0e6),
        Length::from_base(sibling / 1.0e6),
        Length::from_base(d_feed / 1.0e6),
    )?
    .daughter_hematocrit)
}

#[inline]
fn inferred_murray_sibling_diameter(d_daughter: f64, d_feed: f64) -> f64 {
    if d_daughter >= d_feed {
        d_daughter
    } else {
        (d_feed.powi(3) - d_daughter.powi(3))
            .max(d_daughter.powi(3))
            .cbrt()
    }
}

#[inline]
fn logit(x: f64) -> f64 {
    (x / (1.0 - x)).ln()
}

#[inline]
fn inv_logit(x: f64) -> f64 {
    1.0 / (1.0 + (-x).exp())
}

#[cfg(test)]
#[path = "plasma_skimming_tests.rs"]
mod plasma_skimming_tests;
