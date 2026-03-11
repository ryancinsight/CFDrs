//! Conductance-weighted flow-fraction computation for parallel channels.
//!
//! ## Theorem — Conductance-Weighted Flow Splitting (Hagen–Poiseuille)
//!
//! **Theorem**: For *N* equal-height, equal-length rectangular channels sharing a
//! common pressure drop ΔP, the volumetric flow rate through channel *i* is:
//!
//! ```text
//! Qᵢ = Gᵢ · ΔP     where  Gᵢ = wᵢ³ h / (12 μ L · 1/φ(αᵢ))
//! ```
//!
//! When height and length are identical across all branches the correction
//! factor φ and constant prefactors cancel, giving the flow fraction:
//!
//! ```text
//! fᵢ = Qᵢ / Q_total = wᵢ³ / Σⱼ wⱼ³
//! ```
//!
//! **Proof sketch**: The rectangular-duct resistance is R = 12μL / (wh³ φ(α)).
//! The conductance G = 1/R ∝ wh³φ for a single channel.  When h and L are
//! shared, G ∝ w³ (the φ(α) ratio also cancels for equal-height channels
//! because α = h/w varies, but in the common millifluidic case h ≪ w the
//! leading-order conductance is dominated by w³).  Kirchhoff's junction law
//! gives Qᵢ = Gᵢ ΔP, so fᵢ = Gᵢ / Σ Gⱼ = wᵢ³ / Σ wⱼ³.
//!
//! **Validity**: Laminar regime (Re < 2300), equal channel height and length,
//! Newtonian incompressible fluid, fully developed flow.

/// Per-branch flow fractions for a set of parallel rectangular channels of
/// equal height and length.
///
/// Given `n` branch widths, returns `n` fractions that sum to 1.0, weighted
/// by the cubic conductance `w³` (first-order Hagen–Poiseuille for equal-height
/// channels).
///
/// Falls back to linear width weighting when total cubic conductance is
/// negligibly small (< 10⁻³⁰).
///
/// # Returns
///
/// A `Vec<f64>` of per-branch flow fractions in the same order as `widths`.
#[must_use]
pub fn parallel_channel_flow_fractions(widths: &[f64]) -> Vec<f64> {
    if widths.is_empty() {
        return Vec::new();
    }
    let total_cond: f64 = widths.iter().map(|w| w.powi(3)).sum();
    if total_cond > 1.0e-30 {
        widths.iter().map(|w| w.powi(3) / total_cond).collect()
    } else {
        let total_width: f64 = widths.iter().sum();
        if total_width > 0.0 {
            widths.iter().map(|w| w / total_width).collect()
        } else {
            vec![0.0; widths.len()]
        }
    }
}

/// Cumulative treatment-path flow fraction through a cascade of split stages.
///
/// Each "stage" is represented by a slice of `(width, is_treatment_path)` pairs.
/// Returns `(per_stage_fracs, cumulative)` where:
/// - `per_stage_fracs[i]` = fraction of flow entering stage *i* that takes
///   treatment-path branches (via [`parallel_channel_flow_fractions`]).
/// - `cumulative` = product of all per-stage fractions = overall treatment-path
///   volume fraction from inlet to outlet.
///
/// Stages with zero total width are skipped (no contribution to the product).
#[must_use]
pub fn cascade_treatment_flow_fractions(
    stages: &[&[(f64, bool)]],
) -> (Vec<f64>, f64) {
    let mut cumulative = 1.0_f64;
    let mut per_stage = Vec::with_capacity(stages.len());
    for branches in stages {
        let widths: Vec<f64> = branches.iter().map(|(w, _)| *w).collect();
        let total_width: f64 = widths.iter().sum();
        if total_width <= 0.0 {
            continue;
        }
        let fracs = parallel_channel_flow_fractions(&widths);
        let treatment_frac: f64 = branches
            .iter()
            .zip(fracs.iter())
            .filter(|((_, is_treat), _)| *is_treat)
            .map(|(_, f)| f)
            .sum();
        cumulative *= treatment_frac;
        per_stage.push(treatment_frac);
    }
    (per_stage, cumulative)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn equal_widths_give_equal_fractions() {
        let fracs = parallel_channel_flow_fractions(&[1.0, 1.0, 1.0]);
        for f in &fracs {
            assert!((f - 1.0 / 3.0).abs() < 1e-12);
        }
    }

    #[test]
    fn single_branch_gives_unity() {
        let fracs = parallel_channel_flow_fractions(&[0.5]);
        assert_eq!(fracs.len(), 1);
        assert!((fracs[0] - 1.0).abs() < 1e-15);
    }

    #[test]
    fn cubic_scaling_two_to_one() {
        // Branch widths 2:1 → conductance ratio 8:1 → fractions 8/9, 1/9
        let fracs = parallel_channel_flow_fractions(&[2.0, 1.0]);
        assert!((fracs[0] - 8.0 / 9.0).abs() < 1e-12);
        assert!((fracs[1] - 1.0 / 9.0).abs() < 1e-12);
    }

    #[test]
    fn empty_widths() {
        assert!(parallel_channel_flow_fractions(&[]).is_empty());
    }

    #[test]
    fn cascade_two_stages() {
        // Stage 1: treatment branch 2mm, bypass 1mm → treat frac = 8/9
        // Stage 2: treatment branch 1mm, bypass 1mm → treat frac = 1/2
        // Cumulative = 8/9 * 1/2 = 4/9
        let s1: &[(f64, bool)] = &[(2.0, true), (1.0, false)];
        let s2: &[(f64, bool)] = &[(1.0, true), (1.0, false)];
        let (per_stage, cum) = cascade_treatment_flow_fractions(&[s1, s2]);
        assert_eq!(per_stage.len(), 2);
        assert!((per_stage[0] - 8.0 / 9.0).abs() < 1e-12);
        assert!((per_stage[1] - 0.5).abs() < 1e-12);
        assert!((cum - 4.0 / 9.0).abs() < 1e-12);
    }
}
