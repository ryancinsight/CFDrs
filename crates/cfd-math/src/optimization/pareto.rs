//! Pareto-dominance filtering and NSGA-II crowding distance.
//!
//! Provides non-dominated sorting and diversity-preserving crowding distance
//! as used in multi-objective evolutionary algorithms (NSGA-II, MOEA/D).
//!
//! # Complexity
//!
//! [`pareto_front`] uses an O(N² × M) naive non-dominated sort that is
//! practical for N ≤ 50,000 candidates and M ≤ 10 objectives.  For the
//! typical millifluidic design space (≤ 20,000 feasible candidates, 3
//! objectives) this completes in well under one second.
//!
//! # Layout
//!
//! The objective count `M` is const-generic so the candidate layout is a
//! flat `[T; M]` array — no per-candidate heap allocation, cache-coherent
//! traversal. Contrast with the jagged `Vec<Vec<T>>` form this replaced.

use eunomia::RealField;

/// Sense (direction) for a single optimization objective.
///
/// Replaces the `&[bool]` parallel-array encoding, making objective
/// direction explicit and exhaustiveness-checked.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum ObjectiveSense {
    /// Higher values are preferred (maximise).
    Maximize,
    /// Lower values are preferred (minimise).
    Minimize,
}

/// Return the indices of non-dominated (Pareto-optimal) candidates.
///
/// `candidates[i]` is the objective vector for candidate `i`, with `M`
/// objectives laid out contiguously. `senses[j]` specifies whether objective
/// `j` should be maximised or minimised.
///
/// # Panics
///
/// Does not panic. An empty candidate set returns an empty index vector.
///
/// # Examples
///
/// ```
/// use cfd_math::optimization::pareto::{pareto_front, ObjectiveSense};
///
/// // Maximise both objectives: [0.5, 0.5] dominates [0.3, 0.3].
/// let front = pareto_front(
///     &[[0.5_f64, 0.5], [0.3, 0.3]],
///     &[ObjectiveSense::Maximize, ObjectiveSense::Maximize],
/// );
/// assert_eq!(front, vec![0]);
/// ```
#[must_use]
pub fn pareto_front<T, const M: usize>(
    candidates: &[[T; M]],
    senses: &[ObjectiveSense; M],
) -> Vec<usize>
where
    T: PartialOrd + Copy,
{
    let n = candidates.len();
    if n == 0 {
        return Vec::new();
    }

    let mut dominated = vec![false; n];

    for i in 0..n {
        if dominated[i] {
            continue;
        }
        for j in (i + 1)..n {
            if dominated[j] {
                continue;
            }
            match dominance_relation(&candidates[i], &candidates[j], senses) {
                Dominance::IOverJ => dominated[j] = true,
                Dominance::JOverI => {
                    dominated[i] = true;
                    break;
                }
                Dominance::Neither => {}
            }
        }
    }

    (0..n).filter(|&i| !dominated[i]).collect()
}

/// Pareto-dominance relationship between two objective vectors.
#[derive(Debug, PartialEq, Eq)]
enum Dominance {
    /// Candidate `i` dominates candidate `j`.
    IOverJ,
    /// Candidate `j` dominates candidate `i`.
    JOverI,
    /// Neither dominates the other.
    Neither,
}

fn dominance_relation<T, const M: usize>(
    a: &[T; M],
    b: &[T; M],
    senses: &[ObjectiveSense; M],
) -> Dominance
where
    T: PartialOrd + Copy,
{
    let mut a_better = false;
    let mut b_better = false;

    for ((&ai, &bi), sense) in a.iter().zip(b.iter()).zip(senses.iter()) {
        let (a_wins, b_wins) = match sense {
            ObjectiveSense::Maximize => (ai > bi, bi > ai),
            ObjectiveSense::Minimize => (ai < bi, bi < ai),
        };
        if a_wins {
            a_better = true;
        }
        if b_wins {
            b_better = true;
        }
        if a_better && b_better {
            return Dominance::Neither;
        }
    }

    match (a_better, b_better) {
        (true, false) => Dominance::IOverJ,
        (false, true) => Dominance::JOverI,
        _ => Dominance::Neither,
    }
}

/// Compute NSGA-II crowding distances for a Pareto-front member set.
///
/// `front[i]` is the objective vector of the i-th front member with `M`
/// objectives laid out contiguously (same shape as [`pareto_front`] input).
/// Boundary members on each objective axis receive `T::INFINITY`.
/// Higher crowding distance indicates a more isolated (diverse) solution and
/// should be preserved during selection.
#[must_use]
pub fn crowding_distances<T, const M: usize>(front: &[[T; M]]) -> Vec<T>
where
    T: RealField,
{
    let n = front.len();
    if n == 0 {
        return Vec::new();
    }
    let inf = T::infinity();
    if n <= 2 {
        return (0..n).map(|_| inf).collect();
    }

    let mut dist = vec![T::ZERO; n];

    for obj_idx in 0..M {
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            front[a][obj_idx]
                .partial_cmp(&front[b][obj_idx])
                .unwrap_or(core::cmp::Ordering::Equal)
        });

        dist[order[0]] = inf;
        dist[order[n - 1]] = inf;

        let range = front[order[n - 1]][obj_idx] - front[order[0]][obj_idx];
        let epsilon = T::from_f64(1.0e-12);
        if range < epsilon {
            continue;
        }

        for k in 1..(n - 1) {
            let delta = front[order[k + 1]][obj_idx] - front[order[k - 1]][obj_idx];
            dist[order[k]] += delta / range;
        }
    }

    dist
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── pareto_front ──────────────────────────────────────────────────────────

    #[test]
    fn empty_input_returns_empty() {
        let front = pareto_front::<f64, 1>(&[], &[ObjectiveSense::Maximize]);
        assert!(front.is_empty());
    }

    #[test]
    fn single_candidate_is_always_pareto_optimal() {
        let candidates = [[0.5_f64, 0.5]];
        let front = pareto_front(&candidates, &[ObjectiveSense::Maximize; 2]);
        assert_eq!(front, vec![0]);
    }

    #[test]
    fn dominated_point_excluded_maximise() {
        // [0.5, 0.5] dominates [0.3, 0.3] when both are maximised.
        let candidates = [[0.5_f64, 0.5], [0.3, 0.3]];
        let front = pareto_front(&candidates, &[ObjectiveSense::Maximize; 2]);
        assert_eq!(front, vec![0]);
    }

    #[test]
    fn dominated_point_excluded_minimise() {
        // [0.3, 0.3] dominates [0.5, 0.5] when both are minimised.
        let candidates = [[0.5_f64, 0.5], [0.3, 0.3]];
        let front = pareto_front(&candidates, &[ObjectiveSense::Minimize; 2]);
        assert_eq!(front, vec![1]);
    }

    #[test]
    fn two_incomparable_points_both_in_front() {
        // [1.0, 0.0] and [0.0, 1.0] are incomparable (maximise both).
        let candidates = [[1.0_f64, 0.0], [0.0, 1.0]];
        let front = pareto_front(&candidates, &[ObjectiveSense::Maximize; 2]);
        assert_eq!(front.len(), 2);
    }

    #[test]
    fn three_objectives_pareto_front() {
        // [1,1,1] dominates all others; [1,0,1] is non-dominated w.r.t. [0.5,0.5,0.5].
        let candidates = [
            [1.0_f64, 1.0, 1.0], // non-dominated
            [0.5_f64, 0.5, 0.5], // dominated by first
            [1.0_f64, 0.0, 1.0], // dominated by first (obj2: 0 < 1)
        ];
        let front = pareto_front(&candidates, &[ObjectiveSense::Maximize; 3]);
        assert!(front.contains(&0));
        assert!(!front.contains(&1));
    }

    #[test]
    fn pareto_front_is_generic_over_scalar() {
        let candidates_f32 = [[0.5_f32, 0.5], [0.3_f32, 0.3]];
        let front_f32 = pareto_front(&candidates_f32, &[ObjectiveSense::Maximize; 2]);
        assert_eq!(front_f32, vec![0]);

        let candidates_f64 = [[0.5_f64, 0.5], [0.3_f64, 0.3]];
        let front_f64 = pareto_front(&candidates_f64, &[ObjectiveSense::Maximize; 2]);
        assert_eq!(front_f64, vec![0]);
    }

    // ── crowding_distances ────────────────────────────────────────────────────

    #[test]
    fn crowding_empty_returns_empty() {
        let d = crowding_distances::<f64, 2>(&[]);
        assert!(d.is_empty());
    }

    #[test]
    fn crowding_single_is_infinite() {
        let d = crowding_distances(&[[0.5_f64, 0.5]]);
        assert_eq!(d.len(), 1);
        assert!(d[0].is_infinite());
    }

    #[test]
    fn crowding_two_are_both_infinite() {
        let d = crowding_distances(&[[0.0_f64, 1.0], [1.0, 0.0]]);
        assert!(d.iter().all(|v| v.is_infinite()));
    }

    #[test]
    fn crowding_boundary_is_infinite_middle_is_finite() {
        // Three evenly spaced points on the front.
        let front = [[0.0_f64, 1.0], [0.5, 0.5], [1.0, 0.0]];
        let d = crowding_distances(&front);
        assert!(d[0].is_infinite());
        assert!(d[2].is_infinite());
        assert!(d[1].is_finite() && d[1] > 0.0);
    }

    #[test]
    fn crowding_uniform_spacing_gives_equal_interior_distances() {
        // Five evenly spaced points; interior distances should all be equal.
        let front: Vec<[f64; 2]> = (0..5)
            .map(|i| {
                let t = i as f64 / 4.0;
                [t, 1.0 - t]
            })
            .collect();
        let d = crowding_distances(&front);
        assert!(d[0].is_infinite());
        assert!(d[4].is_infinite());
        let first = d[1];
        for &v in &d[1..4] {
            assert!((v - first).abs() < 1e-10, "expected uniform: {d:?}");
        }
    }

    #[test]
    fn crowding_distances_is_generic_over_scalar() {
        let front_f32: Vec<[f32; 2]> = (0..3)
            .map(|i| {
                let t = i as f32 / 2.0;
                [t, 1.0 - t]
            })
            .collect();
        let d_f32 = crowding_distances(&front_f32);
        assert!(d_f32[0].is_infinite());
        assert!(d_f32[2].is_infinite());
        assert!(d_f32[1].is_finite());
    }
}
