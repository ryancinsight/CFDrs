//! Pries et al. (1990) plasma skimming model for bifurcation hematocrit partitioning.
//!
//! ## Theorem — Plasma Skimming (Pries et al. 1990)
//!
//! At a bifurcation, the hematocrit in each daughter branch differs from the
//! feed hematocrit due to the cell-free layer effect:
//!
//! ```text
//! H_daughter / H_feed = FQE^(1/X₀) · (1 − FQE)^(−1/X₀) · B / (B + 1)
//! ```
//!
//! Simplified Pries (1990) logit model:
//!
//! ```text
//! logit(FQ_E) = A + B · logit(FQ_B)
//! ```
//!
//! where:
//! - FQ_B = Q_daughter / Q_total (fractional blood flow to daughter)
//! - FQ_E = (H_daughter · Q_daughter) / (H_feed · Q_total) (fractional RBC flux)
//! - A depends on diameter ratio: A = −13.29 · ((D_α/D_β)² − 1) / ((D_α/D_β)² + 1)
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
//! Also: Pries, A.R. et al. (1990). "Blood flow in microvascular networks",
//! *Circ. Res.* 67:826-834.

/// Compute the daughter-branch hematocrit after plasma skimming at a
/// microvascular bifurcation.
///
/// At bifurcations the cell-free layer near channel walls preferentially
/// enters the smaller daughter branch, reducing its hematocrit below the
/// feed. This simplified empirical model from Pries et al. (1990) captures
/// the essential physics: the daughter receiving less flow gets
/// disproportionately less RBCs.
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
#[must_use]
pub fn plasma_skimming_hematocrit(
    feed_hematocrit: f64,
    flow_fraction: f64,
    diameter_daughter: f64,
    diameter_feed: f64,
) -> f64 {
    let ht_feed = feed_hematocrit.clamp(0.0, 1.0);
    let fq = flow_fraction.clamp(0.0, 1.0);
    let d_daughter = diameter_daughter.max(1.0);
    let d_feed = diameter_feed.max(1.0);

    // Trivial cases
    if ht_feed < 1e-15 || fq < 1e-15 {
        return 0.0;
    }
    if fq > 1.0 - 1e-15 {
        // This branch gets all the flow — hematocrit equals feed
        return ht_feed;
    }

    // Pries (1990) simplified logit model parameters
    // A depends on diameter ratio (asymmetry parameter).
    // In the original Pries formulation, α and β are the two daughter branches.
    // We approximate the other daughter's diameter from Murray's law:
    //   D_other³ = D_feed³ − D_daughter³  →  D_other = (D_feed³ − D_daughter³)^(1/3)
    // When D_daughter ≥ D_feed (e.g. symmetric bifurcation where both daughters
    // have the same diameter as the feed), we assume symmetric split: D_other = D_daughter.
    let a = if d_daughter >= d_feed {
        // Symmetric or larger daughter: A = 0 (no asymmetry bias)
        0.0
    } else {
        // Estimate the other daughter's diameter via Murray's law:
        //   D_other³ = D_feed³ − D_daughter³
        let d_other_cubed = d_feed.powi(3) - d_daughter.powi(3);
        let d_other = d_other_cubed.max(1.0).cbrt();
        // Pries convention: ratio D_other/D_daughter. When the other branch
        // is larger (D_other > D_daughter), A < 0 → less RBC flux to this
        // (smaller) branch, consistent with plasma skimming.
        let d_ratio = d_other / d_daughter;
        let d_ratio_sq = d_ratio * d_ratio;
        -13.29 * (d_ratio_sq - 1.0) / (d_ratio_sq + 1.0)
    };

    // B depends on feed hematocrit and diameter
    let b = 1.0 + 6.98 * (1.0 - ht_feed) / d_feed;

    // logit(FQ_B)
    let logit_fqb = (fq / (1.0 - fq)).ln();

    // logit(FQ_E) = A + B * logit(FQ_B)
    let logit_fqe = a + b * logit_fqb;

    // Invert logit: FQ_E = 1 / (1 + exp(-logit_fqe))
    let fqe = 1.0 / (1.0 + (-logit_fqe).exp());

    // FQ_E is fractional erythrocyte flux to daughter:
    //   FQ_E = (H_daughter * Q_daughter) / (H_feed * Q_total)
    //        = (H_daughter / H_feed) * fq
    // Therefore: H_daughter = H_feed * FQ_E / fq
    let ht_daughter = ht_feed * fqe / fq;

    // Clamp to physically meaningful range
    ht_daughter.clamp(0.0, (2.0 * ht_feed).min(1.0))
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
    fn test_plasma_skimming_equal_split() {
        let ht = plasma_skimming_hematocrit(HT_NORMAL, 0.5, D_FEED, D_FEED);
        let rel_err = (ht - HT_NORMAL).abs() / HT_NORMAL;
        assert!(
            rel_err < 0.05,
            "Equal split hematocrit ({:.4}) should be ~feed ({:.4}), rel_err = {:.4}",
            ht,
            HT_NORMAL,
            rel_err
        );
    }

    /// The smaller daughter branch (receiving 20% of flow) should get
    /// less hematocrit than the feed due to plasma skimming.
    #[test]
    fn test_plasma_skimming_smaller_daughter_less_hematocrit() {
        let ht_small = plasma_skimming_hematocrit(HT_NORMAL, 0.2, 60.0, D_FEED);
        assert!(
            ht_small < HT_NORMAL,
            "Smaller daughter Ht ({:.4}) should be < feed Ht ({:.4})",
            ht_small,
            HT_NORMAL
        );
    }

    /// Daughter hematocrit should always remain in [0, 1] regardless of inputs.
    #[test]
    fn test_plasma_skimming_bounded() {
        let test_cases = [
            (0.45, 0.1, 30.0, 100.0),
            (0.45, 0.9, 90.0, 100.0),
            (0.80, 0.3, 50.0, 100.0),
            (0.10, 0.5, 100.0, 100.0),
            (0.45, 0.01, 20.0, 200.0),
            (0.45, 0.99, 200.0, 200.0),
        ];

        for &(ht, fq, dd, df) in &test_cases {
            let result = plasma_skimming_hematocrit(ht, fq, dd, df);
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
    }

    /// Zero flow fraction means no blood enters this daughter branch,
    /// so its hematocrit should be zero.
    #[test]
    fn test_plasma_skimming_zero_flow_zero_hematocrit() {
        let ht = plasma_skimming_hematocrit(HT_NORMAL, 0.0, 50.0, D_FEED);
        assert!(
            ht.abs() < 1e-15,
            "Zero flow fraction should give zero hematocrit, got {:.10}",
            ht
        );
    }

    /// Higher feed hematocrit should produce higher daughter hematocrit
    /// at the same flow fraction and geometry.
    #[test]
    fn test_plasma_skimming_increases_with_feed_hematocrit() {
        let fq = 0.4;
        let dd = 70.0;

        let ht_low = plasma_skimming_hematocrit(0.20, fq, dd, D_FEED);
        let ht_mid = plasma_skimming_hematocrit(0.35, fq, dd, D_FEED);
        let ht_high = plasma_skimming_hematocrit(0.50, fq, dd, D_FEED);

        assert!(
            ht_low < ht_mid && ht_mid < ht_high,
            "Daughter Ht should increase with feed Ht: {:.4} < {:.4} < {:.4}",
            ht_low,
            ht_mid,
            ht_high
        );
    }
}
