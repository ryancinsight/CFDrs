//! Core Latin Hypercube Sampler.
//!
//! ## Algorithm
//!
//! **Theorem** (McKay, Beckman & Conover, 1979): Latin Hypercube Sampling
//! partitions each of `D` dimensions into `N` equal-probability strata and
//! places exactly one sample per stratum.  For monotone response surfaces,
//! the LHS estimator variance is bounded above by simple random sampling.
//!
//! **Implementation**: For each dimension `d ∈ 0..D`, generate a random
//! permutation `π_d` of `{0, 1, …, N-1}`.  Sample `i` gets value
//! `(π_d[i] + u_d) / N` where `u_d ~ Uniform(0, 1)`.  This produces an
//! `N × D` matrix of values in `[0, 1)^D` with stratification constraints.

use rand::Rng;
use rand::seq::SliceRandom;

/// Latin Hypercube Sampler generating `N` stratified samples in `[0, 1)^D`.
pub struct LatinHypercubeSampler;

impl LatinHypercubeSampler {
    /// Generate `n_samples` points in `[0, 1)^dims`.
    ///
    /// Each returned `Vec<f64>` has length `dims`.  The samples are
    /// stratified: for every dimension, each of the `n_samples` equal
    /// strata contains exactly one sample.
    ///
    /// # Panics
    ///
    /// Panics if `dims == 0` or `n_samples == 0`.
    pub fn generate<R: Rng>(dims: usize, n_samples: usize, rng: &mut R) -> Vec<Vec<f64>> {
        assert!(dims > 0, "LHS requires at least 1 dimension");
        assert!(n_samples > 0, "LHS requires at least 1 sample");

        // Pre-allocate output matrix.
        let mut samples: Vec<Vec<f64>> = (0..n_samples)
            .map(|_| vec![0.0; dims])
            .collect();

        let inv_n = 1.0 / n_samples as f64;

        for d in 0..dims {
            // Build permutation [0, 1, …, n_samples-1].
            let mut perm: Vec<usize> = (0..n_samples).collect();
            perm.shuffle(rng);

            for (i, &stratum) in perm.iter().enumerate() {
                let jitter: f64 = rng.gen();
                samples[i][d] = (stratum as f64 + jitter) * inv_n;
            }
        }

        samples
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use rand::SeedableRng;
    use rand::rngs::StdRng;

    #[test]
    fn unit_interval_bounds() {
        let mut rng = StdRng::seed_from_u64(0xCAFE);
        let samples = LatinHypercubeSampler::generate(4, 100, &mut rng);
        assert_eq!(samples.len(), 100);
        for sample in &samples {
            assert_eq!(sample.len(), 4);
            for &val in sample {
                assert!(val >= 0.0 && val < 1.0, "value {val} out of [0,1)");
            }
        }
    }

    #[test]
    fn stratification_property() {
        // For each dimension, each stratum [k/N, (k+1)/N) should contain
        // exactly one sample.
        let mut rng = StdRng::seed_from_u64(0xBEEF);
        let n = 50;
        let dims = 3;
        let samples = LatinHypercubeSampler::generate(dims, n, &mut rng);

        for d in 0..dims {
            let mut strata = vec![0u32; n];
            for sample in &samples {
                let stratum = (sample[d] * n as f64).floor() as usize;
                let stratum = stratum.min(n - 1);
                strata[stratum] += 1;
            }
            for (k, &count) in strata.iter().enumerate() {
                assert_eq!(
                    count, 1,
                    "dimension {d}, stratum {k}: expected 1 sample, got {count}"
                );
            }
        }
    }

    #[test]
    fn reproducibility() {
        let a = LatinHypercubeSampler::generate(3, 20, &mut StdRng::seed_from_u64(42));
        let b = LatinHypercubeSampler::generate(3, 20, &mut StdRng::seed_from_u64(42));
        for (sa, sb) in a.iter().zip(b.iter()) {
            for (va, vb) in sa.iter().zip(sb.iter()) {
                assert_eq!(va, vb);
            }
        }
    }
}
