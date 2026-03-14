//! Latin Hypercube Sampling for Milestone 12 candidate generation.
//!
//! Replaces the exhaustive grid sweep (6.8M candidates) with stratified
//! random sampling (~50K candidates) across the continuous parameter space
//! of each topology.  Physics-invariant mirror variants are deduplicated
//! to a single `base` evaluation; top-K winners are postfix-replicated.
//!
//! ## Algorithm
//!
//! **Theorem** (McKay, Beckman & Conover, 1979): For a monotone function of
//! any single dimension, the variance of the LHS estimator is never larger
//! than that of simple random sampling, and is usually smaller for smooth
//! responses.  The stratification guarantees every "slice" of each dimension
//! is represented exactly once.

mod sampler;

use std::sync::Arc;
use rand::SeedableRng;
use rand::rngs::StdRng;

use cfd_schematics::topology::presets::enumerate_milestone12_topologies;
use cfd_schematics::TreatmentActuationMode;

use crate::design::space::sweep::milestone12::CandidateParams;
use self::sampler::LatinHypercubeSampler;

/// Default LHS seed for reproducible Milestone 12 runs.
const DEFAULT_LHS_SEED: u64 = 0x4346_4472_735F_4C48;

/// Default number of LHS samples per topology.
///
/// 5 000 provides excellent coverage for the 6–8 continuous dimensions of
/// the selective sweep while reducing the total candidate count from 6.8M
/// to ~50K.
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
/// Acoustic candidates use dimensions 0–4 (5D).
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
#[must_use]
pub fn generate_milestone12_lhs_params() -> Vec<CandidateParams> {
    generate_milestone12_lhs_params_seeded(DEFAULT_LHS_SEED, SAMPLES_PER_TOPOLOGY)
}

/// Seeded variant for testing and custom sample counts.
#[must_use]
pub fn generate_milestone12_lhs_params_seeded(
    seed: u64,
    samples_per_topology: usize,
) -> Vec<CandidateParams> {
    let mut rng = StdRng::seed_from_u64(seed);

    // Deduplicate mirrors: only keep `base` variants (no mirror_x, no
    // mirror_y) for physics evaluation.
    let topologies: Vec<Arc<_>> = enumerate_milestone12_topologies()
        .into_iter()
        .filter(|req| !req.mirror_x && !req.mirror_y)
        .map(Arc::new)
        .collect();

    let total_capacity = topologies.len() * samples_per_topology * 2;
    let mut candidates = Vec::with_capacity(total_capacity);
    let mut idx: u32 = 0;

    for request in &topologies {
        let split_kinds = &request.split_kinds;
        let supports_venturi = split_kinds.iter().all(|sk| {
            matches!(sk, cfd_schematics::SplitKind::NFurcation(2..=5))
        });

        // --- Acoustic tier (ultrasound-only) ---
        let acoustic_samples =
            LatinHypercubeSampler::generate(ACOUSTIC_DIMS, samples_per_topology, &mut rng);
        for sample in &acoustic_samples {
            idx += 1;
            let n_segs = pick_categorical(&N_SEGS_OPTIONS, sample[DIM_Q], &mut rng);
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
            let venturi_samples =
                LatinHypercubeSampler::generate(VENTURI_DIMS, samples_per_topology, &mut rng);
            for sample in &venturi_samples {
                idx += 1;
                let n_segs = pick_categorical(&N_SEGS_OPTIONS, sample[DIM_Q], &mut rng);
                let vt_count = pick_categorical(&VT_COUNT_OPTIONS, sample[DIM_GAUGE], &mut rng);
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

    candidates
}

/// Compute a lightweight config hash for pipeline staleness detection.
///
/// The hash captures the LHS seed, samples-per-topology, and dimension bounds
/// so that stage summaries from a different configuration are detectable at
/// report-assembly time.
#[must_use]
pub fn pipeline_config_hash() -> String {
    use std::hash::{Hash, Hasher};
    let mut hasher = std::collections::hash_map::DefaultHasher::new();
    DEFAULT_LHS_SEED.hash(&mut hasher);
    SAMPLES_PER_TOPOLOGY.hash(&mut hasher);
    FLOW_RANGE.0.to_bits().hash(&mut hasher);
    FLOW_RANGE.1.to_bits().hash(&mut hasher);
    GAUGE_RANGE.0.to_bits().hash(&mut hasher);
    GAUGE_RANGE.1.to_bits().hash(&mut hasher);
    WIDTH_RANGE.0.to_bits().hash(&mut hasher);
    WIDTH_RANGE.1.to_bits().hash(&mut hasher);
    THROAT_RANGE.0.to_bits().hash(&mut hasher);
    THROAT_RANGE.1.to_bits().hash(&mut hasher);
    format!("{:016x}", hasher.finish())
}

/// Linear interpolation: `lo + u * (hi - lo)` where `u ∈ [0, 1]`.
fn lerp((lo, hi): (f64, f64), u: f64) -> f64 {
    lo + u * (hi - lo)
}

/// Pick a categorical value by mapping a `[0, 1)` uniform to one of `N`
/// equally-probable bins.
fn pick_categorical<T: Copy>(options: &[T], u: f64, _rng: &mut StdRng) -> T {
    let index = (u * options.len() as f64).floor() as usize;
    options[index.min(options.len() - 1)]
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
        let a = generate_milestone12_lhs_params_seeded(42, 16);
        let b = generate_milestone12_lhs_params_seeded(42, 16);
        assert_eq!(a.len(), b.len());
        for (pa, pb) in a.iter().zip(b.iter()) {
            assert_eq!(pa.q, pb.q);
            assert_eq!(pa.gauge, pb.gauge);
            assert_eq!(pa.w_ch, pb.w_ch);
        }
    }

    #[test]
    fn lhs_parameters_within_bounds() {
        let params = generate_milestone12_lhs_params();
        for p in &params {
            assert!(p.q >= FLOW_RANGE.0 && p.q <= FLOW_RANGE.1, "flow out of range");
            assert!(p.gauge >= GAUGE_RANGE.0 && p.gauge <= GAUGE_RANGE.1, "gauge out of range");
            assert!(p.w_ch >= WIDTH_RANGE.0 && p.w_ch <= WIDTH_RANGE.1, "width out of range");
            if p.is_venturi() {
                assert!(p.d_throat >= THROAT_RANGE.0 && p.d_throat <= THROAT_RANGE.1);
            }
        }
    }
}
