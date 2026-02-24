//! # Knot Vector
//!
//! A knot vector `Ξ = {ξ₀, ξ₁, …, ξₘ}` is a non-decreasing sequence of
//! parameter values that partitions the parameter domain.  For a B-spline of
//! degree `p` with `n+1` control points, the knot vector has `m+1 = n+p+2`
//! entries.
//!
//! ## Valid knot vectors
//!
//! A knot vector is *valid* for degree `p` and `n+1` control points iff:
//! 1. It is non-decreasing: `ξᵢ ≤ ξᵢ₊₁` for all `i`.
//! 2. It has exactly `n + p + 2` entries.
//!
//! A *clamped* (open) knot vector has the first and last knot repeated `p+1`
//! times, so the curve interpolates its first and last control points.
//!
//! ## Span search
//!
//! Given a parameter value `ξ`, the **active knot span** is the index `i`
//! such that `ξᵢ ≤ ξ < ξᵢ₊₁` and `ξᵢ < ξᵢ₊₁`.
//! (At the right boundary `ξ = ξₙ₊₁` we return `i = n`.)
//!
//! Binary search gives O(log m) lookup.

use crate::domain::core::scalar::Real;

/// A validated knot vector.
///
/// Cloned cheaply — internally backed by a `Vec<Real>`.
#[derive(Clone, Debug, PartialEq)]
pub struct KnotVector {
    knots: Vec<Real>,
}

impl KnotVector {
    /// Create from a raw vector, validating non-decreasingness.
    ///
    /// # Errors
    /// Returns `Err` if the vector is empty or has a decreasing pair.
    pub fn try_new(knots: Vec<Real>) -> Result<Self, KnotError> {
        if knots.is_empty() {
            return Err(KnotError::Empty);
        }
        for i in 0..knots.len() - 1 {
            if knots[i] > knots[i + 1] {
                return Err(KnotError::Decreasing { index: i });
            }
        }
        Ok(Self { knots })
    }

    /// Create without validation (caller guarantees non-decreasingness).
    ///
    /// # Safety (logical, not memory)
    /// Providing a decreasing knot vector leads to incorrect basis evaluation
    /// but not undefined memory behavior.
    pub fn new_unchecked(knots: Vec<Real>) -> Self {
        Self { knots }
    }

    /// Construct a **clamped uniform** knot vector for degree `p` and
    /// `n+1` control points.
    ///
    /// The result is `{0,…,0, 1/(m-2p), 2/(m-2p), …, 1,…,1}` where the
    /// first and last knot are repeated `p+1` times.
    ///
    /// # Panics
    /// Panics if `n < p` (not enough control points for the degree).
    pub fn clamped_uniform(n: usize, p: usize) -> Self {
        assert!(n >= p, "need at least p+1 control points for degree p");
        let m = n + p + 1; // last knot index
        let interior = m - 2 * p; // number of interior knot spans
        let mut knots = Vec::with_capacity(m + 1);
        // p+1 zeros
        for _ in 0..=p {
            knots.push(0.0);
        }
        // interior knots
        for j in 1..interior {
            knots.push(j as Real / interior as Real);
        }
        // p+1 ones
        for _ in 0..=p {
            knots.push(1.0);
        }
        debug_assert_eq!(knots.len(), m + 1);
        Self { knots }
    }

    /// Number of knots (= `n + p + 2` for a valid B-spline).
    #[inline]
    pub fn len(&self) -> usize {
        self.knots.len()
    }

    /// Raw slice of knot values.
    #[inline]
    pub fn as_slice(&self) -> &[Real] {
        &self.knots
    }

    /// Get knot value by index.
    #[inline]
    pub fn get(&self, i: usize) -> Real {
        self.knots[i]
    }

    /// The parameter domain: `[knots[0], knots[last]]`.
    pub fn domain(&self) -> (Real, Real) {
        (*self.knots.first().unwrap(), *self.knots.last().unwrap())
    }

    /// Find the knot span index `i` such that `ξᵢ ≤ ξ < ξᵢ₊₁`.
    ///
    /// Uses binary search — O(log m).
    ///
    /// At the right boundary (`ξ == ξₘ`), returns the last non-zero-length
    /// span.  This makes the parameter domain a half-open interval
    /// `[ξ₀, ξₘ)` with a special-case for the right endpoint.
    ///
    /// # Panics
    /// Panics if `ξ` is outside `[ξ₀, ξₘ]`.
    pub fn find_span(&self, t: Real, n: usize) -> usize {
        let (lo, hi) = self.domain();
        assert!(
            t >= lo && t <= hi,
            "parameter {t} outside knot domain [{lo}, {hi}]"
        );
        // Special case: right boundary
        if t == self.knots[n + 1] {
            // Walk back to the last non-degenerate span
            let mut span = n;
            while span > 0 && self.knots[span] == self.knots[span + 1] {
                span -= 1;
            }
            return span;
        }
        // Binary search in knots[0..=n+1]
        let mut low = 0usize;
        let mut high = n + 1;
        let mut mid = (low + high) / 2;
        while t < self.knots[mid] || t >= self.knots[mid + 1] {
            if t < self.knots[mid] {
                high = mid;
            } else {
                low = mid;
            }
            mid = (low + high) / 2;
        }
        mid
    }

    /// Return the multiplicity of knot value `t` (number of times it appears).
    pub fn multiplicity(&self, t: Real) -> usize {
        self.knots
            .iter()
            .filter(|&&k| (k - t).abs() < 1e-14)
            .count()
    }

    /// Return `true` if this is a clamped (open) knot vector for degree `p`:
    /// the first `p+1` entries are equal and the last `p+1` entries are equal.
    pub fn is_clamped(&self, p: usize) -> bool {
        let m = self.knots.len() - 1;
        if m < 2 * p {
            return false;
        }
        let first = self.knots[0];
        let last = self.knots[m];
        self.knots[..=p].iter().all(|&k| (k - first).abs() < 1e-14)
            && self.knots[m - p..]
                .iter()
                .all(|&k| (k - last).abs() < 1e-14)
    }
}

/// Error produced when constructing an invalid [`KnotVector`].
#[derive(Clone, Debug, PartialEq)]
pub enum KnotError {
    /// Empty knot vector.
    Empty,
    /// Knot at `index+1` is less than knot at `index`.
    Decreasing {
        /// Index of the first offending pair.
        index: usize,
    },
    /// The number of knots is inconsistent with degree and control-point count.
    WrongLength {
        /// Number of knots provided.
        got: usize,
        /// Number of knots expected.
        expected: usize,
    },
}

impl std::fmt::Display for KnotError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            KnotError::Empty => write!(f, "knot vector is empty"),
            KnotError::Decreasing { index } => {
                write!(f, "knot vector is decreasing at index {}", index)
            }
            KnotError::WrongLength { got, expected } => {
                write!(f, "knot vector has {} entries; expected {}", got, expected)
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn empty_fails() {
        assert_eq!(KnotVector::try_new(vec![]), Err(KnotError::Empty));
    }

    #[test]
    fn decreasing_fails() {
        let r = KnotVector::try_new(vec![0.0, 1.0, 0.5]);
        assert!(matches!(r, Err(KnotError::Decreasing { index: 1 })));
    }

    #[test]
    fn clamped_uniform_cubic() {
        // Cubic (p=3) with 6 control points (n=5): m = 5+3+1 = 9, 10 knots
        let kv = KnotVector::clamped_uniform(5, 3);
        assert_eq!(kv.len(), 10);
        assert!(kv.is_clamped(3));
        assert_eq!(kv.get(0), 0.0);
        assert_eq!(kv.get(9), 1.0);
    }

    #[test]
    fn find_span_basic() {
        // Linear (p=1), 3 control points: knots = [0, 0, 0.5, 1, 1]
        let kv = KnotVector::try_new(vec![0.0, 0.0, 0.5, 1.0, 1.0]).unwrap();
        // n = 3 - 1 = 2 (number of control points - 1)
        let n = 2;
        assert_eq!(kv.find_span(0.0, n), 1);
        assert_eq!(kv.find_span(0.25, n), 1);
        assert_eq!(kv.find_span(0.5, n), 2);
        assert_eq!(kv.find_span(0.75, n), 2);
        assert_eq!(kv.find_span(1.0, n), 2); // right boundary
    }

    #[test]
    fn multiplicity() {
        let kv = KnotVector::try_new(vec![0.0, 0.0, 0.0, 0.5, 1.0, 1.0, 1.0]).unwrap();
        assert_eq!(kv.multiplicity(0.0), 3);
        assert_eq!(kv.multiplicity(0.5), 1);
        assert_eq!(kv.multiplicity(1.0), 3);
    }
}
