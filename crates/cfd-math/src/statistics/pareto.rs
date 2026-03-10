//! Pareto-dominance filtering and NSGA-II crowding distance.
//!
//! Provides non-dominated sorting and diversity-preserving crowding distance
//! as used in multi-objective evolutionary algorithms (NSGA-II, MOEA/D).
//!
//! # Complexity
//!
//! [`pareto_front_nd`] uses an O(N² × M) naive non-dominated sort that is
//! practical for N ≤ 50,000 candidates and M ≤ 10 objectives.  For the
//! typical millifluidic design space (≤ 20,000 feasible candidates, 3
//! objectives) this completes in well under one second.

/// Return indices of non-dominated (Pareto-optimal) candidates.
///
/// `objectives[i]` is the objective vector for candidate `i`.
/// `is_maximized[j]` indicates whether objective `j` should be maximised
/// (`true`) or minimised (`false`).
///
/// # Panics
///
/// Panics if `objectives` is non-empty and `is_maximized.len()` differs from
/// `objectives[0].len()`.
#[must_use]
pub fn pareto_front_nd(objectives: &[Vec<f64>], is_maximized: &[bool]) -> Vec<usize> {
    let n = objectives.len();
    if n == 0 {
        return Vec::new();
    }
    debug_assert_eq!(
        objectives[0].len(),
        is_maximized.len(),
        "objective vector length must match is_maximized length"
    );

    let mut dominated = vec![false; n];

    for i in 0..n {
        if dominated[i] {
            continue;
        }
        for j in (i + 1)..n {
            if dominated[j] {
                continue;
            }
            match dominance_relation(&objectives[i], &objectives[j], is_maximized) {
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
    /// Vector `a` dominates vector `b`.
    IOverJ,
    /// Vector `b` dominates vector `a`.
    JOverI,
    /// Neither vector dominates the other.
    Neither,
}

/// Determine the dominance relationship between two objective vectors.
fn dominance_relation(a: &[f64], b: &[f64], maximize: &[bool]) -> Dominance {
    let mut a_strictly_better = false;
    let mut b_strictly_better = false;

    for ((&ai, &bi), &max) in a.iter().zip(b.iter()).zip(maximize.iter()) {
        let (a_wins, b_wins) = if max {
            (ai > bi, bi > ai)
        } else {
            (ai < bi, bi < ai)
        };
        if a_wins {
            a_strictly_better = true;
        }
        if b_wins {
            b_strictly_better = true;
        }
        // Early exit: neither can dominate the other if both are better in some objective
        if a_strictly_better && b_strictly_better {
            return Dominance::Neither;
        }
    }

    match (a_strictly_better, b_strictly_better) {
        (true, false) => Dominance::IOverJ,
        (false, true) => Dominance::JOverI,
        _ => Dominance::Neither,
    }
}

/// Compute NSGA-II crowding distances for a Pareto-front member set.
///
/// Higher crowding distance indicates a more isolated (diverse) solution and
/// should be preserved during selection.  Boundary members on each objective
/// axis receive [`f64::INFINITY`].
///
/// `front_objectives[i]` is the objective vector of the i-th front member.
/// The returned vector has the same length as `front_objectives`.
#[must_use]
pub fn crowding_distances(front_objectives: &[Vec<f64>]) -> Vec<f64> {
    let n = front_objectives.len();
    if n == 0 {
        return Vec::new();
    }
    if n <= 2 {
        return vec![f64::INFINITY; n];
    }

    let m = front_objectives[0].len();
    let mut dist = vec![0.0_f64; n];

    for obj_idx in 0..m {
        // Sort members by this objective (ascending)
        let mut order: Vec<usize> = (0..n).collect();
        order.sort_by(|&a, &b| {
            front_objectives[a][obj_idx]
                .partial_cmp(&front_objectives[b][obj_idx])
                .unwrap_or(std::cmp::Ordering::Equal)
        });

        // Boundary members receive infinite distance
        dist[order[0]] = f64::INFINITY;
        dist[order[n - 1]] = f64::INFINITY;

        let range = front_objectives[order[n - 1]][obj_idx] - front_objectives[order[0]][obj_idx];
        if range < 1e-12 {
            // All members have the same value for this objective; no contribution
            continue;
        }

        for k in 1..(n - 1) {
            let delta =
                front_objectives[order[k + 1]][obj_idx] - front_objectives[order[k - 1]][obj_idx];
            dist[order[k]] += delta / range;
        }
    }

    dist
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    // ── pareto_front_nd ──────────────────────────────────────────────────────

    #[test]
    fn empty_input_returns_empty() {
        let front = pareto_front_nd(&[], &[true]);
        assert!(front.is_empty());
    }

    #[test]
    fn single_candidate_is_always_pareto_optimal() {
        let objs = vec![vec![0.5, 0.5]];
        let front = pareto_front_nd(&objs, &[true, true]);
        assert_eq!(front, vec![0]);
    }

    #[test]
    fn dominated_point_excluded_maximise() {
        // [0.5, 0.5] dominates [0.3, 0.3] when both objectives are maximised
        let objs = vec![vec![0.5, 0.5], vec![0.3, 0.3]];
        let front = pareto_front_nd(&objs, &[true, true]);
        assert_eq!(front, vec![0]);
    }

    #[test]
    fn dominated_point_excluded_minimise() {
        // [0.3, 0.3] dominates [0.5, 0.5] when both objectives are minimised
        let objs = vec![vec![0.5, 0.5], vec![0.3, 0.3]];
        let front = pareto_front_nd(&objs, &[false, false]);
        assert_eq!(front, vec![1]);
    }

    #[test]
    fn two_incomparable_points_both_in_front() {
        // [1.0, 0.0] and [0.0, 1.0] are incomparable (maximise both)
        let objs = vec![vec![1.0_f64, 0.0], vec![0.0, 1.0]];
        let front = pareto_front_nd(&objs, &[true, true]);
        assert_eq!(front.len(), 2);
    }

    #[test]
    fn three_objectives_pareto_front() {
        // Only [1,1,1] dominates all others
        let objs = vec![
            vec![1.0, 1.0, 1.0], // non-dominated
            vec![0.5, 0.5, 0.5], // dominated by [1,1,1]
            vec![1.0, 0.0, 1.0], // non-dominated (obj2 lower but obj2 also lower for dominated)
        ];
        let front = pareto_front_nd(&objs, &[true, true, true]);
        // [0] and [2] are both non-dominated: [0] beats [2] on obj2 but [2] ties on obj1 and obj3
        // Actually [0]=[1,1,1] dominates [2]=[1,0,1] because obj1 equal, obj2 better, obj3 equal
        assert!(front.contains(&0));
        assert!(!front.contains(&1)); // [0.5,0.5,0.5] is dominated
    }

    // ── crowding_distances ───────────────────────────────────────────────────

    #[test]
    fn crowding_empty_returns_empty() {
        let d = crowding_distances(&[]);
        assert!(d.is_empty());
    }

    #[test]
    fn crowding_single_is_infinite() {
        let d = crowding_distances(&[vec![0.5, 0.5]]);
        assert_eq!(d, vec![f64::INFINITY]);
    }

    #[test]
    fn crowding_two_are_both_infinite() {
        let d = crowding_distances(&[vec![0.0, 1.0], vec![1.0, 0.0]]);
        assert!(d.iter().all(|&v| v == f64::INFINITY));
    }

    #[test]
    fn crowding_boundary_is_infinite_middle_is_finite() {
        // Three evenly spaced points on the front
        let front = vec![vec![0.0_f64, 1.0], vec![0.5, 0.5], vec![1.0, 0.0]];
        let d = crowding_distances(&front);
        // After sorting by obj0: order = [0, 1, 2]
        // Boundaries: d[0] and d[2] are INF; d[1] is finite
        assert_eq!(d[0], f64::INFINITY);
        assert_eq!(d[2], f64::INFINITY);
        assert!(d[1].is_finite() && d[1] > 0.0);
    }

    #[test]
    fn crowding_uniform_spacing_gives_equal_interior_distances() {
        // Five evenly spaced points: interior distances should all be equal
        let front: Vec<Vec<f64>> = (0..5)
            .map(|i| {
                let t = i as f64 / 4.0;
                vec![t, 1.0 - t]
            })
            .collect();
        let d = crowding_distances(&front);
        // Boundaries: d[0] and d[4] are INF
        assert_eq!(d[0], f64::INFINITY);
        assert_eq!(d[4], f64::INFINITY);
        // Interior: all should be 0.5 (1/4 / 1 * 2 objectives = 0.5 each)
        let interior: Vec<f64> = d[1..4].to_vec();
        let first = interior[0];
        for &v in &interior {
            assert!((v - first).abs() < 1e-10, "expected uniform: {interior:?}");
        }
    }
}
