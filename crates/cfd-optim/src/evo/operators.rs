//! Genetic operators: selection, crossover, and mutation.

use rand::Rng;

/// Tournament selection: return the index of the winner among `k` random participants.
pub(super) fn tournament_select<R: Rng>(fitnesses: &[f64], k: usize, rng: &mut R) -> usize {
    let n = fitnesses.len();
    let mut best = rng.gen_range(0..n);
    for _ in 1..k {
        let challenger = rng.gen_range(0..n);
        if fitnesses[challenger] > fitnesses[best] {
            best = challenger;
        }
    }
    best
}

/// Simulated Binary Crossover (SBX, Deb & Agrawal 1995).
///
/// Produces two offspring from two parents.  `eta` is the distribution index
/// (larger = children closer to parents).
pub(super) fn sbx_crossover<R: Rng>(
    p1: &[f64],
    p2: &[f64],
    eta: f64,
    rng: &mut R,
) -> (Vec<f64>, Vec<f64>) {
    let mut c1 = p1.to_vec();
    let mut c2 = p2.to_vec();

    for i in 0..p1.len() {
        if rng.gen::<f64>() < 0.5 {
            let u = rng.gen::<f64>();
            let beta = if u <= 0.5 {
                (2.0 * u).powf(1.0 / (eta + 1.0))
            } else {
                (1.0 / (2.0 * (1.0 - u))).powf(1.0 / (eta + 1.0))
            };
            c1[i] = (0.5 * ((1.0 + beta) * p1[i] + (1.0 - beta) * p2[i])).clamp(0.0, 1.0);
            c2[i] = (0.5 * ((1.0 - beta) * p1[i] + (1.0 + beta) * p2[i])).clamp(0.0, 1.0);
        }
    }
    (c1, c2)
}

/// Probability that a single offspring's gene[0] (topology selector) is
/// completely re-randomised to a uniform draw from [0, 1].
pub(super) const TOPO_JUMP_RATE: f64 = 0.25;

/// Apply a topology jump to gene[0]: with probability [`TOPO_JUMP_RATE`],
/// replace it with a uniform draw from [0, 1].
pub(super) fn apply_topology_jump<R: Rng>(genes: &mut Vec<f64>, rng: &mut R) {
    if rng.gen::<f64>() < TOPO_JUMP_RATE {
        genes[0] = rng.gen::<f64>();
    }
}

/// Polynomial mutation (Deb & Deb 2012).
///
/// Mutates each gene with probability `p_m`.  `eta_m` is the mutation
/// distribution index (larger = smaller perturbations).
pub(super) fn polynomial_mutation<R: Rng>(genes: &mut [f64], eta_m: f64, p_m: f64, rng: &mut R) {
    for g in genes.iter_mut() {
        if rng.gen::<f64>() < p_m {
            let u = rng.gen::<f64>();
            let delta = if u < 0.5 {
                (2.0 * u).powf(1.0 / (eta_m + 1.0)) - 1.0
            } else {
                1.0 - (2.0 * (1.0 - u)).powf(1.0 / (eta_m + 1.0))
            };
            *g = (*g + delta).clamp(0.0, 1.0);
        }
    }
}
