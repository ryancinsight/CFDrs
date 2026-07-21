//! Latin Hypercube Sampling for Milestone 12 candidate generation.
//!
//! Replaces the exhaustive grid sweep (6.8M candidates) with a Tyche-owned
//! stratified design (~50K candidates) across the continuous parameter space
//! of each topology. Physics-invariant mirror variants are deduplicated to a
//! single `base` evaluation; top-K winners are postfix-replicated.
//!
//! Tyche is the single source of truth for the Latin-hypercube algorithm and
//! its stratification proof. This module owns only the mapping from unit-space
//! samples into CFD candidate parameters.

use core::num::NonZeroU32;
use std::sync::Arc;

use cfd_schematics::topology::presets::{
    enumerate_milestone12_topologies, Milestone12TopologyRequest,
};
use cfd_schematics::TreatmentActuationMode;
use tyche_core::{sampling::Counter, sampling::UserDomain, Design, LatinHypercube, Seed, SplitMix64};

use crate::design::space::sweep::milestone12::CandidateParams;
use crate::error::OptimError;

/// Default LHS seed for reproducible Milestone 12 runs.
const DEFAULT_LHS_SEED: u64 = 0x4346_4472_735F_4C48;

/// Default number of LHS samples per topology.
///
/// Release builds evaluate 5 000 points per acoustic or venturi tier. The
/// final candidate count is the validated product of this count and the
/// topology-tier count.
#[cfg(not(any(test, debug_assertions)))]
const SAMPLES_PER_TOPOLOGY: usize = 5_000;
#[cfg(any(test, debug_assertions))]
const SAMPLES_PER_TOPOLOGY: usize = 32;

/// Static dimension bounds drawn from `Milestone12Dimensions`.
const FLOW_RANGE: (f64, f64) = (1.333e-6, 3.333e-6);
const GAUGE_RANGE: (f64, f64) = (25_000.0, 200_000.0);
const WIDTH_RANGE: (f64, f64) = (4.0e-3, 8.0e-3);
const THROAT_RANGE: (f64, f64) = (35.0e-6, 120.0e-6);
const TL_FACTOR_RANGE: (f64, f64) = (2.0, 10.0);
const PRETRI_FRAC_RANGE: (f64, f64) = (0.33, 0.55);
const TRI_CENTER_FRAC_RANGE: (f64, f64) = (0.250, 0.650);
const BI_TREAT_FRAC_RANGE: (f64, f64) = (0.60, 0.76);

/// Categorical values for segment count.
const N_SEGS_OPTIONS: [usize; 3] = [1, 5, 7];
/// Categorical values for venturi count.
const VT_COUNT_OPTIONS: [u8; 4] = [1, 2, 3, 4];

/// Index map for the 8-dimensional continuous LHS sample vector.
///
/// Acoustic candidates use dimensions 0–5 (6D).
/// Venturi candidates use all 8 dimensions.
const DIM_Q: usize = 0;
const DIM_GAUGE: usize = 1;
const DIM_WIDTH: usize = 2;
const DIM_PRETRI: usize = 3;
const DIM_TRI: usize = 4;
const DIM_BI: usize = 5;
const DIM_THROAT: usize = 6;
const DIM_TL_FACTOR: usize = 7;

/// Number of continuous dimensions for acoustic (non-venturi) candidates.
const ACOUSTIC_DIMS: usize = 6;
/// Number of continuous dimensions for venturi candidates.
const VENTURI_DIMS: usize = 8;

/// Generate Milestone 12 candidate parameters via Latin Hypercube Sampling.
///
/// Only `base` mirror variants are evaluated (mirrors are physics-invariant).
/// The caller is responsible for postfix-replicating top-K winners into
/// all 4 mirror variants for report SVG/figure output.
///
/// # Examples
///
/// ```
/// use cfd_optim::generate_milestone12_lhs_params;
///
/// let candidates = generate_milestone12_lhs_params();
/// assert!(candidates
///     .iter()
///     .all(|candidate| !candidate.request.mirror_x && !candidate.request.mirror_y));
/// ```
#[must_use]
pub fn generate_milestone12_lhs_params() -> Vec<CandidateParams> {
    generate_milestone12_lhs_params_seeded(DEFAULT_LHS_SEED, SAMPLES_PER_TOPOLOGY)
        .expect("invariant: the compile-time default produces a valid candidate capacity")
}

/// Seeded variant for testing and custom sample counts.
///
/// Each topology tier receives a distinct counter-derived Tyche seed. The
/// generated unit-space point is held in one fixed array and mapped directly
/// into a candidate, without an intermediate design matrix.
///
/// # Errors
///
/// Returns [`OptimError::InvalidParameter`] when the requested sample count is
/// zero, exceeds Tyche's `u32` design bound, or would overflow candidate IDs.
///
/// # Examples
///
/// ```
/// use cfd_optim::generate_milestone12_lhs_params_seeded;
///
/// let first = generate_milestone12_lhs_params_seeded(42, 8)?;
/// let replay = generate_milestone12_lhs_params_seeded(42, 8)?;
/// assert_eq!(first[0].q, replay[0].q);
/// assert_eq!(first[0].gauge, replay[0].gauge);
/// # Ok::<(), cfd_optim::OptimError>(())
/// ```
pub fn generate_milestone12_lhs_params_seeded(
    seed: u64,
    samples_per_topology: usize,
) -> Result<Vec<CandidateParams>, OptimError> {
    let sample_count = validated_sample_count(samples_per_topology)?;
    let root_seed = Seed::new(seed);
    let mut design_ordinal = 0_u64;

    // Deduplicate mirrors: only keep `base` variants (no mirror_x, no
    // mirror_y) for physics evaluation.
    let topologies: Vec<Arc<_>> = enumerate_milestone12_topologies()
        .into_iter()
        .filter(|req| !req.mirror_x && !req.mirror_y)
        .map(Arc::new)
        .collect();

    let venturi_tiers = topologies
        .iter()
        .filter(|request| supports_venturi(request))
        .count();
    let tier_count = topologies
        .len()
        .checked_add(venturi_tiers)
        .ok_or_else(|| OptimError::InvalidParameter("candidate tier count overflowed".into()))?;
    let total_capacity = validated_candidate_capacity(tier_count, sample_count)?;
    let mut candidates = Vec::new();
    candidates.try_reserve_exact(total_capacity).map_err(|_| {
        OptimError::InvalidParameter(format!(
            "unable to reserve storage for {total_capacity} candidates"
        ))
    })?;
    let mut idx: u32 = 0;

    for request in &topologies {
        let supports_venturi = supports_venturi(request);

        // --- Acoustic tier (ultrasound-only) ---
        let acoustic_design = next_design::<ACOUSTIC_DIMS>(root_seed, design_ordinal, sample_count);
        design_ordinal = design_ordinal
            .checked_add(1)
            .expect("invariant: design tier count fits u64");
        for sample_index in 0..acoustic_design.sample_count() {
            let mut sample = [0.0; ACOUSTIC_DIMS];
            acoustic_design
                .sample_unit_into(sample_index, &mut sample)
                .expect("invariant: loop index is inside the Tyche design");
            idx = idx
                .checked_add(1)
                .expect("invariant: prevalidated candidate count fits u32");
            let n_segs = pick_categorical(&N_SEGS_OPTIONS, sample[DIM_Q]);
            candidates.push(CandidateParams {
                idx,
                request: Arc::clone(request),
                q: lerp(FLOW_RANGE, sample[DIM_Q]),
                gauge: lerp(GAUGE_RANGE, sample[DIM_GAUGE]),
                d_throat: 0.0,
                throat_len: 0.0,
                w_ch: lerp(WIDTH_RANGE, sample[DIM_WIDTH]),
                n_segs,
                pretri_center_frac: lerp(PRETRI_FRAC_RANGE, sample[DIM_PRETRI]),
                terminal_tri_center_frac: lerp(TRI_CENTER_FRAC_RANGE, sample[DIM_TRI]),
                bi_treat_frac: lerp(BI_TREAT_FRAC_RANGE, sample[DIM_BI]),
                treatment_actuation_mode: TreatmentActuationMode::UltrasoundOnly,
                vt_count: 0,
            });
        }

        // --- Venturi tier (cavitation) ---
        if supports_venturi {
            let venturi_design =
                next_design::<VENTURI_DIMS>(root_seed, design_ordinal, sample_count);
            design_ordinal = design_ordinal
                .checked_add(1)
                .expect("invariant: design tier count fits u64");
            for sample_index in 0..venturi_design.sample_count() {
                let mut sample = [0.0; VENTURI_DIMS];
                venturi_design
                    .sample_unit_into(sample_index, &mut sample)
                    .expect("invariant: loop index is inside the Tyche design");
                idx = idx
                    .checked_add(1)
                    .expect("invariant: prevalidated candidate count fits u32");
                let n_segs = pick_categorical(&N_SEGS_OPTIONS, sample[DIM_Q]);
                let vt_count = pick_categorical(&VT_COUNT_OPTIONS, sample[DIM_GAUGE]);
                let d_throat = lerp(THROAT_RANGE, sample[DIM_THROAT]);
                let tl_factor = lerp(TL_FACTOR_RANGE, sample[DIM_TL_FACTOR]);
                candidates.push(CandidateParams {
                    idx,
                    request: Arc::clone(request),
                    q: lerp(FLOW_RANGE, sample[DIM_Q]),
                    gauge: lerp(GAUGE_RANGE, sample[DIM_GAUGE]),
                    d_throat,
                    throat_len: d_throat * tl_factor,
                    w_ch: lerp(WIDTH_RANGE, sample[DIM_WIDTH]),
                    n_segs,
                    pretri_center_frac: lerp(PRETRI_FRAC_RANGE, sample[DIM_PRETRI]),
                    terminal_tri_center_frac: lerp(TRI_CENTER_FRAC_RANGE, sample[DIM_TRI]),
                    bi_treat_frac: lerp(BI_TREAT_FRAC_RANGE, sample[DIM_BI]),
                    treatment_actuation_mode: TreatmentActuationMode::VenturiCavitation,
                    vt_count,
                });
            }
        }
    }

    Ok(candidates)
}

fn validated_sample_count(samples: usize) -> Result<NonZeroU32, OptimError> {
    let bounded = u32::try_from(samples).map_err(|_| {
        OptimError::InvalidParameter(format!(
            "samples per topology {samples} exceeds Tyche's u32 design bound"
        ))
    })?;
    NonZeroU32::new(bounded)
        .ok_or_else(|| OptimError::InvalidParameter("samples per topology must be non-zero".into()))
}

fn validated_candidate_capacity(
    tier_count: usize,
    sample_count: NonZeroU32,
) -> Result<usize, OptimError> {
    let samples = usize::try_from(sample_count.get()).map_err(|_| {
        OptimError::InvalidParameter("target usize cannot represent the Tyche sample count".into())
    })?;
    let capacity = tier_count.checked_mul(samples).ok_or_else(|| {
        OptimError::InvalidParameter("candidate capacity overflowed usize".into())
    })?;
    u32::try_from(capacity).map_err(|_| {
        OptimError::InvalidParameter(format!(
            "candidate count {capacity} exceeds the u32 identifier space"
        ))
    })?;
    Ok(capacity)
}

fn next_design<const PARAMETERS: usize>(
    root_seed: Seed,
    ordinal: u64,
    sample_count: NonZeroU32,
) -> LatinHypercube<PARAMETERS, SplitMix64> {
    let seed = Seed::new(Counter::<UserDomain<0>, SplitMix64>::word(root_seed, ordinal, 0));
    LatinHypercube::new(seed, sample_count)
}

fn supports_venturi(request: &Milestone12TopologyRequest) -> bool {
    request
        .split_kinds
        .iter()
        .all(|kind| matches!(kind, cfd_schematics::SplitKind::NFurcation(2..=5)))
}

/// Linear interpolation: `lo + u * (hi - lo)` where `u ∈ [0, 1)`.
fn lerp((lo, hi): (f64, f64), u: f64) -> f64 {
    lo + u * (hi - lo)
}

/// Pick a categorical value by mapping a `[0, 1)` uniform to one of `N`
/// equally-probable bins.
fn pick_categorical<T: Copy, const OPTIONS: usize>(options: &[T; OPTIONS], unit: f64) -> T {
    assert!(OPTIONS > 0, "invariant: categorical options are non-empty");
    let count = u32::try_from(OPTIONS).expect("invariant: categorical option count fits in u32");
    for index in 1..OPTIONS {
        let boundary =
            u32::try_from(index).expect("invariant: categorical option index fits in u32");
        if unit < f64::from(boundary) / f64::from(count) {
            return options[index - 1];
        }
    }
    options[OPTIONS - 1]
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lhs_generator_produces_both_modes() {
        let params = generate_milestone12_lhs_params();
        assert!(!params.is_empty(), "LHS should produce candidates");
        let has_acoustic = params.iter().any(|p| !p.is_venturi());
        let has_venturi = params.iter().any(|p| p.is_venturi());
        assert!(has_acoustic, "LHS should include acoustic candidates");
        assert!(has_venturi, "LHS should include venturi candidates");
    }

    #[test]
    fn lhs_generator_produces_every_requested_tier_sample() {
        const SAMPLES: usize = 7;
        let topologies: Vec<_> = enumerate_milestone12_topologies()
            .into_iter()
            .filter(|request| !request.mirror_x && !request.mirror_y)
            .collect();
        let venturi_tiers = topologies
            .iter()
            .filter(|request| supports_venturi(request))
            .count();
        let expected = (topologies.len() + venturi_tiers) * SAMPLES;

        let params = generate_milestone12_lhs_params_seeded(42, SAMPLES).unwrap();

        assert_eq!(params.len(), expected);
        assert_eq!(
            params.last().map(|candidate| candidate.idx),
            u32::try_from(expected).ok()
        );
    }

    #[test]
    fn lhs_generator_only_produces_base_mirrors() {
        let params = generate_milestone12_lhs_params();
        for p in &params {
            assert!(
                !p.request.mirror_x && !p.request.mirror_y,
                "LHS should only produce base mirror variants, got {:?}",
                p.request.design_name,
            );
        }
    }

    #[test]
    fn lhs_generator_is_reproducible() {
        let a = generate_milestone12_lhs_params_seeded(42, 16).unwrap();
        let b = generate_milestone12_lhs_params_seeded(42, 16).unwrap();
        assert_eq!(a.len(), b.len());
        for (pa, pb) in a.iter().zip(b.iter()) {
            assert_eq!(pa.q, pb.q);
            assert_eq!(pa.gauge, pb.gauge);
            assert_eq!(pa.w_ch, pb.w_ch);
        }
    }

    #[test]
    fn lhs_seed_changes_continuous_parameters() {
        let first = generate_milestone12_lhs_params_seeded(1, 16).unwrap();
        let second = generate_milestone12_lhs_params_seeded(2, 16).unwrap();

        assert!(
            first
                .iter()
                .zip(&second)
                .any(|(left, right)| left.q != right.q || left.gauge != right.gauge),
            "distinct study seeds must change at least one continuous coordinate"
        );
    }

    #[test]
    fn lhs_rejects_zero_samples() {
        let error = generate_milestone12_lhs_params_seeded(42, 0).unwrap_err();

        assert_eq!(
            error.to_string(),
            "invalid parameter: samples per topology must be non-zero"
        );
    }

    #[test]
    fn lhs_parameters_within_bounds() {
        let params = generate_milestone12_lhs_params();
        for p in &params {
            assert!(
                p.q >= FLOW_RANGE.0 && p.q < FLOW_RANGE.1,
                "flow out of range"
            );
            assert!(
                p.gauge >= GAUGE_RANGE.0 && p.gauge < GAUGE_RANGE.1,
                "gauge out of range"
            );
            assert!(
                p.w_ch >= WIDTH_RANGE.0 && p.w_ch < WIDTH_RANGE.1,
                "width out of range"
            );
            if p.is_venturi() {
                assert!(p.d_throat >= THROAT_RANGE.0 && p.d_throat < THROAT_RANGE.1);
            }
        }
    }
}
