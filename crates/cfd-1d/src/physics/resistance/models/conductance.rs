//! Conductance-weighted flow-fraction computation for parallel channels.
//!
//! ## Theorem — Conductance-Weighted Flow Splitting (Hagen–Poiseuille)
//!
//! **Theorem**: For *N* rectangular channels sharing a common pressure drop ΔP
//! and the same axial length `L`, the volumetric flow rate through channel *i* is:
//!
//! ```text
//! Qᵢ = Gᵢ · ΔP     where  Gᵢ = 2 Aᵢ D_h,ᵢ² / (Poᵢ μ L)
//! ```
//!
//! with:
//!
//! ```text
//! Aᵢ      = wᵢ hᵢ
//! D_h,ᵢ   = 2 wᵢ hᵢ / (wᵢ + hᵢ)
//! Poᵢ     = Po(εᵢ)
//! εᵢ      = min(wᵢ, hᵢ) / max(wᵢ, hᵢ)
//! ```
//!
//! Therefore the flow fraction is:
//!
//! ```text
//! fᵢ = Qᵢ / Q_total = Gᵢ / Σⱼ Gⱼ
//! ```
//!
//! **Proof sketch**: The laminar rectangular-duct resistance can be written as
//! `R = Po μ L / (2 A D_h²)`, where `Po` is the aspect-ratio-dependent
//! Poiseuille number (Bahrami et al. 2006 / Shah-London limit). Since
//! `G = 1/R`, each branch flow under a shared `ΔP` is proportional to
//! `A D_h² / Po`. Kirchhoff's junction law then gives `Qᵢ = Gᵢ ΔP`, so the
//! normalized flow split is `fᵢ = Gᵢ / Σ Gⱼ`.
//!
//! **Validity**: Laminar regime (Re < 2300), fully developed flow, smooth walls,
//! Newtonian incompressible fluid. In the wide-channel limit `w ≫ h`, this
//! reduces to approximately linear width weighting (`G ∝ w`) rather than a
//! generic `w³` law.

use std::f64::consts::PI;

fn rectangular_poiseuille_number(width_m: f64, height_m: f64) -> f64 {
    let width = width_m.max(0.0);
    let height = height_m.max(0.0);
    if width <= 0.0 || height <= 0.0 {
        return 0.0;
    }
    let epsilon = width.min(height) / width.max(height);
    if epsilon <= 1.0e-18 {
        return 96.0;
    }

    let denominator = (1.0 + epsilon).powi(2)
        * (1.0 - (192.0 * epsilon / PI.powi(5)) * (PI / (2.0 * epsilon)).tanh());
    if denominator <= 1.0e-18 {
        0.0
    } else {
        96.0 / denominator
    }
}

fn rectangular_conductance_factor(width_m: f64, height_m: f64) -> f64 {
    let width = width_m.max(0.0);
    let height = height_m.max(0.0);
    if width <= 0.0 || height <= 0.0 {
        return 0.0;
    }

    let area = width * height;
    let hydraulic_diameter = 2.0 * width * height / (width + height).max(1.0e-18);
    let poiseuille_number = rectangular_poiseuille_number(width, height);
    if !poiseuille_number.is_finite() || poiseuille_number <= 1.0e-18 {
        0.0
    } else {
        2.0 * area * hydraulic_diameter.powi(2) / poiseuille_number
    }
}

/// Per-branch flow fractions for a set of parallel rectangular channels of
/// equal axial length.
///
/// Given `n` rectangular branch dimensions `(width, height)`, returns `n`
/// fractions that sum to 1.0, weighted by the laminar rectangular conductance
/// factor `A D_h² / Po(ε)`.
///
/// Falls back to area weighting when total conductance is negligibly small
/// (< 10⁻³⁰).
///
/// # Returns
///
/// A `Vec<f64>` of per-branch flow fractions in the same order as `branches`.
#[must_use]
pub fn parallel_channel_flow_fractions(branches: &[(f64, f64)]) -> Vec<f64> {
    if branches.is_empty() {
        return Vec::new();
    }

    let total_cond: f64 = branches
        .iter()
        .map(|(width, height)| rectangular_conductance_factor(*width, *height))
        .sum();
    if total_cond > 1.0e-30 {
        branches
            .iter()
            .map(|(width, height)| rectangular_conductance_factor(*width, *height) / total_cond)
            .collect()
    } else {
        let total_area: f64 = branches
            .iter()
            .map(|(width, height)| width.max(0.0) * height.max(0.0))
            .sum();
        if total_area > 0.0 {
            branches
                .iter()
                .map(|(width, height)| width.max(0.0) * height.max(0.0) / total_area)
                .collect()
        } else {
            vec![0.0; branches.len()]
        }
    }
}

/// Cumulative treatment-path flow fraction through a cascade of split stages.
///
/// Each "stage" is represented by a slice of `(width, height, is_treatment)`
/// triples.
/// Returns `(per_stage_fracs, cumulative)` where:
/// - `per_stage_fracs[i]` = fraction of flow entering stage *i* that takes
///   treatment-path branches (via [`parallel_channel_flow_fractions`]).
/// - `cumulative` = product of all per-stage fractions = overall treatment-path
///   volume fraction from inlet to outlet.
///
/// Stages with zero total area are skipped (no contribution to the product).
#[must_use]
pub fn cascade_treatment_flow_fractions(stages: &[&[(f64, f64, bool)]]) -> (Vec<f64>, f64) {
    let mut cumulative = 1.0_f64;
    let mut per_stage = Vec::with_capacity(stages.len());
    for branches in stages {
        let branch_dimensions: Vec<(f64, f64)> = branches
            .iter()
            .map(|(width, height, _)| (*width, *height))
            .collect();
        let total_area: f64 = branch_dimensions
            .iter()
            .map(|(width, height)| width.max(0.0) * height.max(0.0))
            .sum();
        if total_area <= 0.0 {
            continue;
        }
        let fracs = parallel_channel_flow_fractions(&branch_dimensions);
        let treatment_frac: f64 = branches
            .iter()
            .zip(fracs.iter())
            .filter(|((_, _, is_treat), _)| *is_treat)
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
    fn equal_geometries_give_equal_fractions() {
        let fracs = parallel_channel_flow_fractions(&[(1.0, 0.25), (1.0, 0.25), (1.0, 0.25)]);
        for f in &fracs {
            assert!((f - 1.0 / 3.0).abs() < 1e-12);
        }
    }

    #[test]
    fn single_branch_gives_unity() {
        let fracs = parallel_channel_flow_fractions(&[(0.5, 0.1)]);
        assert_eq!(fracs.len(), 1);
        assert!((fracs[0] - 1.0).abs() < 1.0e-15);
    }

    #[test]
    fn wide_duct_limit_tracks_linear_width_scaling() {
        let fracs = parallel_channel_flow_fractions(&[(4.0, 0.01), (2.0, 0.01), (1.0, 0.01)]);
        let expected = [4.0 / 7.0, 2.0 / 7.0, 1.0 / 7.0];
        for (actual, expected) in fracs.iter().zip(expected) {
            assert!((actual - expected).abs() < 1.0e-3);
        }
    }

    #[test]
    fn wider_branch_still_draws_more_flow_at_fixed_height() {
        let fracs = parallel_channel_flow_fractions(&[(2.0, 1.0), (1.0, 1.0)]);
        assert!(fracs[0] > fracs[1]);
        assert!((fracs.iter().sum::<f64>() - 1.0).abs() < 1.0e-12);
    }

    #[test]
    fn empty_branches() {
        assert!(parallel_channel_flow_fractions(&[]).is_empty());
    }

    #[test]
    fn cascade_two_stages() {
        let s1: &[(f64, f64, bool)] = &[(2.0, 0.001, true), (1.0, 0.001, false)];
        let s2: &[(f64, f64, bool)] = &[(1.0, 0.001, true), (1.0, 0.001, false)];
        let (per_stage, cum) = cascade_treatment_flow_fractions(&[s1, s2]);
        assert_eq!(per_stage.len(), 2);
        assert!((per_stage[0] - 2.0 / 3.0).abs() < 5.0e-4);
        assert!((per_stage[1] - 0.5).abs() < 1e-12);
        assert!((cum - 1.0 / 3.0).abs() < 5.0e-4);
    }

    #[test]
    fn negative_dimensions_clamped_to_zero() {
        let fracs = parallel_channel_flow_fractions(&[(-1.0, 1.0), (1.0, 1.0)]);
        assert_eq!(fracs[0], 0.0);
        assert_eq!(fracs[1], 1.0);

        let fracs2 = parallel_channel_flow_fractions(&[(-2.0, -1.0), (-1.0, -0.5)]);
        assert_eq!(fracs2[0], 0.0);
        assert_eq!(fracs2[1], 0.0);
    }
}
