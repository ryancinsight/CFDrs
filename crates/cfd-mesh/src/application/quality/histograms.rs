//! Fixed-width histogram for mesh quality distributions.
//!
//! ## Complexity
//!
//! | Operation | Complexity | Notes |
//! |-----------|------------|-------|
//! | `compute` | O(n)       | n = value count |
//! | `percentile` | O(n log n) | sort-based; result is approximate (bin midpoint) |
//!
//! ## References
//!
//! Sturges, H.A. (1926). "The Choice of a Class Interval". *JASA* 21(153): 65-66.

use crate::domain::core::scalar::Real;

// ── Histogram ─────────────────────────────────────────────────────────────────

/// Fixed-width histogram over a set of scalar values.
///
/// `edges` has `n_bins + 1` elements: `edges[i]` is the left boundary of bin `i`,
/// and `edges[n_bins]` is the right boundary of the last bin.
/// `bins[i]` is the count of values in `[edges[i], edges[i+1])`.
/// Values exactly equal to `max` fall into the last bin.
#[derive(Clone, Debug)]
pub struct Histogram {
    /// Per-bin counts.
    pub bins: Vec<usize>,
    /// Bin boundary values; length = `bins.len() + 1`.
    pub edges: Vec<Real>,
    /// Minimum value observed.
    pub min: Real,
    /// Maximum value observed.
    pub max: Real,
}

impl Histogram {
    /// Build a fixed-width histogram from `values` with `n_bins` bins.
    ///
    /// Returns `None` when `values` is empty or `n_bins` is zero.
    ///
    /// All finite values are bucketed.  Non-finite values (NaN, ±∞) are silently
    /// skipped so that degenerate triangle metrics do not corrupt the histogram.
    #[must_use]
    pub fn compute(values: &[Real], n_bins: usize) -> Option<Self> {
        if values.is_empty() || n_bins == 0 {
            return None;
        }
        let finite: Vec<Real> = values.iter().copied().filter(|v| v.is_finite()).collect();
        if finite.is_empty() {
            return None;
        }
        let min = finite.iter().copied().fold(Real::INFINITY, Real::min);
        let max = finite.iter().copied().fold(Real::NEG_INFINITY, Real::max);

        let mut edges = Vec::with_capacity(n_bins + 1);
        let range = max - min;
        let bin_w = if range < Real::EPSILON {
            1.0
        } else {
            range / n_bins as Real
        };
        for i in 0..=n_bins {
            edges.push(min + i as Real * bin_w);
        }
        // Ensure the last edge is exactly max (floating-point safety).
        *edges.last_mut().unwrap() = max + Real::EPSILON;

        let mut bins = vec![0usize; n_bins];
        let inv_w = 1.0 / bin_w;
        for &v in &finite {
            let idx = ((v - min) * inv_w).floor() as usize;
            let idx = idx.min(n_bins - 1);
            bins[idx] += 1;
        }

        Some(Histogram {
            bins,
            edges,
            min,
            max,
        })
    }

    /// Number of bins.
    #[must_use]
    pub fn n_bins(&self) -> usize {
        self.bins.len()
    }

    /// Bin width (uniform).
    #[must_use]
    pub fn bin_width(&self) -> Real {
        if self.edges.len() < 2 {
            return 0.0;
        }
        self.edges[1] - self.edges[0]
    }

    /// Bin midpoint for bin index `i`.
    #[must_use]
    pub fn midpoint(&self, i: usize) -> Real {
        0.5 * (self.edges[i] + self.edges[i + 1])
    }

    /// Approximate p-th percentile (0.0 = min, 1.0 = max) using bin midpoints.
    ///
    /// Returns `None` if the histogram is empty.
    #[must_use]
    pub fn percentile(&self, p: f64) -> Option<Real> {
        let total: usize = self.bins.iter().sum();
        if total == 0 {
            return None;
        }
        let target = (p.clamp(0.0, 1.0) * total as f64) as usize;
        let mut cumulative = 0usize;
        for (i, &count) in self.bins.iter().enumerate() {
            cumulative += count;
            if cumulative > target {
                return Some(self.midpoint(i));
            }
        }
        Some(self.midpoint(self.bins.len() - 1))
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn histogram_uniform_values_fills_all_bins() {
        let values: Vec<Real> = (0..100).map(|i| i as Real).collect();
        let h = Histogram::compute(&values, 10).unwrap();
        assert_eq!(h.bins.len(), 10);
        // Each bin should have 10 values.
        for &count in &h.bins {
            assert!(count > 0, "all bins should be non-empty for uniform data");
        }
        assert_eq!(h.bins.iter().sum::<usize>(), 100);
    }

    #[test]
    fn histogram_single_value_all_in_one_bin() {
        let values = vec![42.0_f64; 50];
        let h = Histogram::compute(&values, 5).unwrap();
        assert_eq!(h.bins.iter().sum::<usize>(), 50);
    }

    #[test]
    fn histogram_nan_is_skipped() {
        let values = vec![1.0, f64::NAN, 2.0, 3.0];
        let h = Histogram::compute(&values, 3).unwrap();
        assert_eq!(h.bins.iter().sum::<usize>(), 3);
    }

    #[test]
    fn histogram_empty_returns_none() {
        assert!(Histogram::compute(&[], 5).is_none());
    }

    #[test]
    fn histogram_zero_bins_returns_none() {
        assert!(Histogram::compute(&[1.0, 2.0], 0).is_none());
    }

    #[test]
    fn histogram_percentile_median_is_central() {
        // 0..100 uniform, 10 bins → p50 should be near 50.
        let values: Vec<Real> = (0..100).map(|i| i as Real).collect();
        let h = Histogram::compute(&values, 10).unwrap();
        let p50 = h.percentile(0.5).unwrap();
        assert!(
            (p50 - 50.0).abs() < 15.0,
            "median should be near 50, got {p50}"
        );
    }
}
