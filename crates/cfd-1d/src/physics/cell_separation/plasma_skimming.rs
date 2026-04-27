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
    diameter_daughter_alpha: f64,
    diameter_daughter_beta: f64,
    diameter_feed: f64,
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
    diameter_daughter_alpha: f64,
    diameter_daughter_beta: f64,
    diameter_feed: f64,
) -> Result<PhaseSeparationResult> {
    if !feed_hematocrit.is_finite()
        || !flow_fraction.is_finite()
        || !diameter_daughter_alpha.is_finite()
        || !diameter_daughter_beta.is_finite()
        || !diameter_feed.is_finite()
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
    if diameter_daughter_alpha <= 0.0 || diameter_daughter_beta <= 0.0 || diameter_feed <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Pries phase-separation diameters must be positive".to_string(),
        ));
    }

    let h_feed = feed_hematocrit;
    let d_feed = diameter_feed;
    let d_alpha = diameter_daughter_alpha;
    let d_beta = diameter_daughter_beta;

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

    let x0 = secomb_phase_separation_x0(d_feed, h_feed)?;
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
/// * `diameter_daughter` - Daughter branch diameter [µm]
/// * `diameter_feed` - Feed (parent) branch diameter [µm]
///
/// # Returns
/// Daughter branch hematocrit, clamped to [0, min(1, 2 × H_feed)].
#[inline]
pub fn plasma_skimming_hematocrit(
    feed_hematocrit: f64,
    flow_fraction: f64,
    diameter_daughter: f64,
    diameter_feed: f64,
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
    diameter_daughter: f64,
    diameter_feed: f64,
) -> Result<f64> {
    if !feed_hematocrit.is_finite()
        || !flow_fraction.is_finite()
        || !diameter_daughter.is_finite()
        || !diameter_feed.is_finite()
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
    if diameter_daughter <= 0.0 || diameter_feed <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "Plasma-skimming diameters must be positive".to_string(),
        ));
    }

    let ht_feed = feed_hematocrit;
    let fq = flow_fraction;
    let d_daughter = diameter_daughter;
    let d_feed = diameter_feed;

    if ht_feed < 1e-15 || fq < 1e-15 {
        return Ok(0.0);
    }
    if fq >= 1.0 - 1e-15 {
        return Ok(ht_feed);
    }

    let sibling = inferred_murray_sibling_diameter(d_daughter, d_feed);
    Ok(
        checked_pries_phase_separation(ht_feed, fq, d_daughter, sibling, d_feed)?
            .daughter_hematocrit,
    )
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
mod tests {
    use super::*;

    const HT_NORMAL: f64 = 0.45;
    const D_FEED: f64 = 100.0; // 100 µm parent vessel

    /// Equal 50/50 flow split with equal diameters: daughter hematocrit
    /// should be approximately equal to feed hematocrit (symmetric case,
    /// no skimming bias).
    #[test]
    fn test_plasma_skimming_equal_split() -> cfd_core::error::Result<()> {
        let ht = plasma_skimming_hematocrit(HT_NORMAL, 0.5, D_FEED, D_FEED)?;
        let rel_err = (ht - HT_NORMAL).abs() / HT_NORMAL;
        assert!(
            rel_err < 0.05,
            "Equal split hematocrit ({:.4}) should be ~feed ({:.4}), rel_err = {:.4}",
            ht,
            HT_NORMAL,
            rel_err
        );
        Ok(())
    }

    /// The smaller daughter branch (receiving 20% of flow) should get
    /// less hematocrit than the feed due to plasma skimming.
    #[test]
    fn test_plasma_skimming_smaller_daughter_less_hematocrit() -> cfd_core::error::Result<()> {
        let ht_small = plasma_skimming_hematocrit(HT_NORMAL, 0.2, 60.0, D_FEED)?;
        assert!(
            ht_small < HT_NORMAL,
            "Smaller daughter Ht ({:.4}) should be < feed Ht ({:.4})",
            ht_small,
            HT_NORMAL
        );
        Ok(())
    }

    /// Daughter hematocrit should always remain in [0, 1] regardless of inputs.
    #[test]
    fn test_plasma_skimming_bounded() -> cfd_core::error::Result<()> {
        let test_cases = [
            (0.45, 0.1, 30.0, 100.0),
            (0.45, 0.9, 90.0, 100.0),
            (0.80, 0.3, 50.0, 100.0),
            (0.10, 0.5, 100.0, 100.0),
            (0.45, 0.01, 20.0, 200.0),
            (0.45, 0.99, 200.0, 200.0),
        ];

        for &(ht, fq, dd, df) in &test_cases {
            let result = plasma_skimming_hematocrit(ht, fq, dd, df)?;
            assert!(
                (0.0..=1.0).contains(&result),
                "Ht={}, fq={}, dd={}, df={} → result {} out of [0,1]",
                ht,
                fq,
                dd,
                df,
                result
            );
        }
        Ok(())
    }

    /// Zero flow fraction means no blood enters this daughter branch,
    /// so its hematocrit should be zero.
    #[test]
    fn test_plasma_skimming_zero_flow_zero_hematocrit() -> cfd_core::error::Result<()> {
        let ht = plasma_skimming_hematocrit(HT_NORMAL, 0.0, 50.0, D_FEED)?;
        assert!(
            ht.abs() < 1e-15,
            "Zero flow fraction should give zero hematocrit, got {:.10}",
            ht
        );
        Ok(())
    }

    /// Higher feed hematocrit should produce higher daughter hematocrit
    /// at the same flow fraction and geometry.
    #[test]
    fn test_plasma_skimming_increases_with_feed_hematocrit() -> cfd_core::error::Result<()> {
        let fq = 0.4;
        let dd = 70.0;

        let ht_low = plasma_skimming_hematocrit(0.20, fq, dd, D_FEED)?;
        let ht_mid = plasma_skimming_hematocrit(0.35, fq, dd, D_FEED)?;
        let ht_high = plasma_skimming_hematocrit(0.50, fq, dd, D_FEED)?;

        assert!(
            ht_low < ht_mid && ht_mid < ht_high,
            "Daughter Ht should increase with feed Ht: {:.4} < {:.4} < {:.4}",
            ht_low,
            ht_mid,
            ht_high
        );
        Ok(())
    }

    #[test]
    fn test_pries_phase_separation_enforces_x0_threshold() -> cfd_core::error::Result<()> {
        let result = pries_phase_separation(HT_NORMAL, 0.01, 30.0, 90.0, 20.0)?;
        assert!(
            result.x0 > 0.01,
            "Expected meaningful X0, got {:.4}",
            result.x0
        );
        assert_eq!(result.cell_fraction, 0.0);
        assert_eq!(result.daughter_hematocrit, 0.0);
        Ok(())
    }

    #[test]
    fn test_pries_phase_separation_biases_wider_daughter() -> cfd_core::error::Result<()> {
        let wide = pries_phase_separation(HT_NORMAL, 0.55, 90.0, 40.0, 100.0)?;
        let narrow = pries_phase_separation(HT_NORMAL, 0.45, 40.0, 90.0, 100.0)?;
        assert!(
            wide.cell_fraction > narrow.cell_fraction,
            "Wider daughter should receive more RBC flux than the narrower sibling: wide={:.4}, narrow={:.4}",
            wide.cell_fraction,
            narrow.cell_fraction
        );
        Ok(())
    }

    #[test]
    fn test_checked_pries_matches_legacy_nominal_case() -> cfd_core::error::Result<()> {
        let legacy = pries_phase_separation(HT_NORMAL, 0.55, 90.0, 40.0, 100.0)?;
        let checked = checked_pries_phase_separation(HT_NORMAL, 0.55, 90.0, 40.0, 100.0)?;

        assert!((legacy.cell_fraction - checked.cell_fraction).abs() < 1e-12);
        assert!((legacy.daughter_hematocrit - checked.daughter_hematocrit).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn test_checked_compact_wrapper_rejects_nonphysical_diameter() {
        let err = checked_plasma_skimming_hematocrit(HT_NORMAL, 0.5, 0.0, D_FEED)
            .expect_err("checked plasma-skimming wrapper must reject zero daughter diameter");
        assert!(err.to_string().contains("diameters"));
    }

    #[test]
    fn test_checked_compact_wrapper_matches_legacy_nominal_case() -> cfd_core::error::Result<()> {
        let legacy = plasma_skimming_hematocrit(HT_NORMAL, 0.4, 60.0, D_FEED)?;
        let checked = checked_plasma_skimming_hematocrit(HT_NORMAL, 0.4, 60.0, D_FEED)?;

        assert!((legacy - checked).abs() < 1e-12);
        Ok(())
    }

    #[test]
    fn test_compact_wrapper_preserves_pries_x0_threshold() -> cfd_core::error::Result<()> {
        let ht = plasma_skimming_hematocrit(HT_NORMAL, 0.01, 30.0, 20.0)?;
        assert_eq!(
            ht, 0.0,
            "Compact wrapper must inherit the Pries cell-entry threshold"
        );
        Ok(())
    }

    #[test]
    fn test_legacy_wrapper_monotone_for_validation_triplet() -> cfd_core::error::Result<()> {
        let small = plasma_skimming_hematocrit(HT_NORMAL, 0.2, 30.0, D_FEED)?;
        let medium = plasma_skimming_hematocrit(HT_NORMAL, 0.4, 60.0, D_FEED)?;
        let large = plasma_skimming_hematocrit(HT_NORMAL, 0.6, 90.0, D_FEED)?;
        assert!(
            small < medium && medium < large,
            "Legacy wrapper should preserve monotone hematocrit ranking for typical design-space scans: small={small:.4}, medium={medium:.4}, large={large:.4}"
        );
        Ok(())
    }
}
