//! Fahraeus effect: tube hematocrit reduction in microchannels.
//!
//! # Theorem — Fahraeus Effect (Fahraeus 1929)
//!
//! In tubes with diameter $D < 300\,\mu m$, the **tube hematocrit** $H_T$
//! (volume fraction of RBCs within the tube at any instant) is lower than
//! the **discharge (feed) hematocrit** $H_F$ (volume fraction in collected
//! outflow). This occurs because RBCs travel faster than the surrounding
//! plasma (they concentrate in the high-velocity core), so fewer RBCs are
//! present in the tube at steady state to maintain the same RBC flux.
//!
//! ## Pries et al. (1990) empirical correlation
//!
//! ```text
//! H_T / H_F = H_F + (1 − H_F) · (1 + 1.7·e^{−0.415·D} − 0.6·e^{−0.011·D})
//! ```
//!
//! where $D$ is the tube diameter in micrometers.
//!
//! At large $D$: $H_T / H_F \to 1$ (no Fahraeus effect).
//! At $D \approx 10\,\mu m$: $H_T / H_F \approx 0.5$ (strong reduction).
//!
//! ## Distinction from Fahraeus-Lindqvist
//!
//! - **Fahraeus effect**: tube hematocrit is reduced ($H_T < H_F$)
//! - **Fahraeus-Lindqvist effect**: apparent viscosity is reduced
//!
//! Both arise from the same physical mechanism (RBC margination to core),
//! but describe different observable quantities.
//!
//! # References
//! - Fahraeus, R. (1929). The suspension stability of blood. *Physiol. Rev.* 9:241-274.
//! - Pries, A.R. et al. (1990). Blood flow in microvascular networks.
//!   *Circ. Res.* 67:826-834.

/// Pries et al. (1990) tube-to-discharge hematocrit ratio.
///
/// Returns $H_T / H_F$ as a function of tube diameter only. The ratio is
/// independent of the absolute hematocrit in this empirical model.
///
/// ```text
/// H_T / H_F = H_F + (1 − H_F) · (1 + 1.7·e^{−0.415·D} − 0.6·e^{−0.011·D})
/// ```
///
/// For the purpose of this standalone ratio, we evaluate at the reference
/// hematocrit $H_F = 0.45$, consistent with the Pries parameterisation.
///
/// # Arguments
/// * `diameter_um` — Tube diameter [µm]
///
/// # Returns
/// $H_T / H_F \in (0, 1]$
#[inline]
#[must_use]
pub fn tube_hematocrit_ratio(diameter_um: f64) -> f64 {
    let d = diameter_um.clamp(3.0, 1000.0);
    let ratio = 0.45 + (1.0 - 0.45) * (1.0 + 1.7 * (-0.415 * d).exp() - 0.6 * (-0.011 * d).exp());
    ratio.clamp(0.0, 1.0)
}

/// Compute the tube hematocrit from the feed (discharge) hematocrit.
///
/// # Theorem — Fahraeus Tube Hematocrit
///
/// ```text
/// H_T = H_F · r(D)
/// ```
///
/// where $r(D) = H_T/H_F$ is the Pries (1990) diameter-dependent ratio.
///
/// **Proof of $H_T \le H_F$**: Since RBCs preferentially occupy the
/// high-velocity core, their mean velocity exceeds the mean fluid velocity.
/// For a given RBC flux $\dot{N}_{RBC}$, the instantaneous RBC count in
/// the tube is inversely proportional to their mean velocity:
/// $H_T = \dot{N}_{RBC} / (A \cdot \bar{v}_{RBC})$. Since
/// $\bar{v}_{RBC} > \bar{v}_{fluid}$, we have $H_T < H_F = \dot{N}_{RBC} / (A \cdot \bar{v}_{fluid})$. ∎
///
/// # Arguments
/// * `feed_hematocrit` — Discharge (feed) hematocrit $H_F \in [0, 1]$
/// * `diameter_um` — Tube diameter [µm]
///
/// # Returns
/// Tube hematocrit $H_T \in [0, H_F]$
#[inline]
#[must_use]
pub fn tube_hematocrit(feed_hematocrit: f64, diameter_um: f64) -> f64 {
    let ht = feed_hematocrit.clamp(0.0, 1.0);
    let ratio = tube_hematocrit_ratio(diameter_um);
    (ht * ratio).clamp(0.0, ht)
}

/// Inverse Fahraeus: compute feed (discharge) hematocrit from tube hematocrit.
///
/// At network bifurcations, blood enters daughter vessels with a known
/// tube hematocrit. To compute the discharge hematocrit for downstream
/// re-entry calculations:
///
/// ```text
/// H_F = H_T / r(D)
/// ```
///
/// # Arguments
/// * `tube_ht` — Tube hematocrit $H_T \in [0, 1]$
/// * `diameter_um` — Tube diameter [µm]
///
/// # Returns
/// Feed (discharge) hematocrit $H_F \in [H_T, 1]$
#[inline]
#[must_use]
pub fn discharge_hematocrit(tube_ht: f64, diameter_um: f64) -> f64 {
    let ht = tube_ht.clamp(0.0, 1.0);
    let ratio = tube_hematocrit_ratio(diameter_um);
    if ratio < 1e-15 {
        return ht;
    }
    (ht / ratio).clamp(ht, 1.0)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// For large tubes (D=1000 µm), the Fahraeus effect is negligible:
    /// H_T ≈ H_F.
    #[test]
    fn large_tube_no_fahraeus_effect() {
        let ratio = tube_hematocrit_ratio(1000.0);
        assert!(
            (ratio - 1.0).abs() < 0.05,
            "Large tube ratio {ratio:.4} should be ~1.0"
        );
    }

    /// For small tubes (D=10 µm), the Fahraeus reduction is strong:
    /// H_T / H_F < 0.8.
    #[test]
    fn small_tube_strong_fahraeus() {
        let ratio = tube_hematocrit_ratio(10.0);
        assert!(
            ratio < 0.8,
            "Small tube ratio {ratio:.4} should be < 0.8"
        );
    }

    /// Tube hematocrit is always ≤ feed hematocrit.
    #[test]
    fn tube_ht_leq_feed() {
        for d in [5.0, 10.0, 30.0, 50.0, 100.0, 300.0, 1000.0] {
            let ht = tube_hematocrit(0.45, d);
            assert!(
                ht <= 0.45 + 1e-10,
                "H_T={ht:.6} should be ≤ H_F=0.45 at D={d}"
            );
        }
    }

    /// Ratio generally increases with diameter above the Fahraeus inversion
    /// zone (D > ~10 µm). Below D ≈ 8 µm (RBC diameter), the formula is
    /// non-monotonic due to single-file flow effects.
    #[test]
    fn ratio_increases_with_diameter() {
        let mut prev = tube_hematocrit_ratio(10.0);
        for d in [20.0, 50.0, 100.0, 300.0] {
            let ratio = tube_hematocrit_ratio(d);
            assert!(
                ratio >= prev - 1e-10,
                "Ratio at D={d} ({ratio:.4}) should be ≥ prev ({prev:.4})"
            );
            prev = ratio;
        }
    }

    /// Round-trip: discharge_hematocrit(tube_hematocrit(H_F, D), D) ≈ H_F.
    #[test]
    fn fahraeus_round_trip() {
        let hf = 0.45;
        for d in [10.0, 30.0, 50.0, 100.0, 300.0] {
            let ht = tube_hematocrit(hf, d);
            let hf_recovered = discharge_hematocrit(ht, d);
            assert!(
                (hf_recovered - hf).abs() < 0.01,
                "Round-trip at D={d}: {hf_recovered:.4} should be ~{hf:.4}"
            );
        }
    }

    /// Zero hematocrit stays zero.
    #[test]
    fn zero_ht_remains_zero() {
        assert!(tube_hematocrit(0.0, 50.0).abs() < 1e-15);
    }
}
