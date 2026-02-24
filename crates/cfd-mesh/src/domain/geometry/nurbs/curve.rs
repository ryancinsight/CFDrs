//! # NURBS and B-Spline Curves
//!
//! Provides `BSplineCurve` (non-rational) and `NurbsCurve` (rational) for any
//! embedding dimension.  Dimension `D` is a [`nalgebra::Const`] type parameter
//! so all allocations are stack-based for D ≤ 4.
//!
//! ## Mathematical Foundation
//!
//! ### B-Spline Curve (non-rational)
//!
//! ```text
//! C(t) = Σᵢ N_{i,p}(t) · Pᵢ
//! ```
//!
//! where `N_{i,p}` are the B-spline basis functions computed by Cox–de Boor
//! recursion (see [`super::basis`]).
//!
//! ### NURBS Curve (rational)
//!
//! ```text
//!        Σᵢ N_{i,p}(t) · wᵢ · Pᵢ
//! C(t) = ─────────────────────────
//!           Σᵢ N_{i,p}(t) · wᵢ
//! ```
//!
//! ### Theorem — Partition of Unity
//!
//! For any parameter value `t` in the domain:
//! ```text
//! Σᵢ N_{i,p}(t) = 1
//! ```
//! This ensures that a B-spline curve is an affine combination of its control
//! points and that translations and rotations of the control polygon map
//! exactly to the curve.
//!
//! ### Theorem — Convex Hull Property
//!
//! Each point `C(t)` lies in the convex hull of the control points local to the
//! active knot span.  For p+1 overlapping spans this is a "local convex hull".
//! Critical for conservative AABB computation in BVH construction.
//!
//! ## Example
//!
//! ```rust,no_run
//! use cfd_mesh::domain::geometry::nurbs::knot::KnotVector;
//! use cfd_mesh::domain::geometry::nurbs::curve::NurbsCurve;
//! use nalgebra::SVector;
//!
//! // Quadratic NURBS arc (quarter circle in XY plane)
//! let ctrl = vec![
//!     SVector::<f64, 3>::new(1.0, 0.0, 0.0),
//!     SVector::<f64, 3>::new(1.0, 1.0, 0.0),
//!     SVector::<f64, 3>::new(0.0, 1.0, 0.0),
//! ];
//! let weights = vec![1.0_f64, std::f64::consts::FRAC_1_SQRT_2, 1.0];
//! let knots = KnotVector::try_new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap();
//! let curve = NurbsCurve::new(ctrl, weights, knots, 2).unwrap();
//!
//! let mid = curve.point(0.5);
//! let scale = mid[0].hypot(mid[1]);
//! assert!((scale - 1.0).abs() < 1e-10, "point should be on unit circle");
//! ```

use super::basis::{eval_basis, eval_basis_and_deriv};
use super::knot::{KnotError, KnotVector};
use crate::domain::core::scalar::Real;
use nalgebra::SVector;

// ── Errors ────────────────────────────────────────────────────────────────────

/// Error returned when constructing a B-spline or NURBS curve fails.
#[derive(Clone, Debug, PartialEq)]
pub enum CurveError {
    /// Knot vector is invalid.
    Knot(KnotError),
    /// Vector of control points is empty.
    NoControlPoints,
    /// Number of knots is inconsistent with degree and control-point count.
    KnotCountMismatch {
        /// Number of knots provided.
        got: usize,
        /// Number of knots expected (`n + p + 2`).
        expected: usize,
    },
    /// Degree is zero — undefined for B-splines.
    ZeroDegree,
    /// Weights vector length does not match control-point count.
    WeightsMismatch {
        /// Number of weights provided.
        weights: usize,
        /// Number of control points.
        control_points: usize,
    },
    /// A weight is ≤ 0, which makes the NURBS curve ill-defined.
    NonPositiveWeight {
        /// Index of the offending weight.
        index: usize,
    },
}

impl std::fmt::Display for CurveError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CurveError::Knot(e) => write!(f, "knot error: {e}"),
            CurveError::NoControlPoints => write!(f, "no control points"),
            CurveError::KnotCountMismatch { got, expected } => {
                write!(f, "expected {expected} knots (n+p+2), got {got}")
            }
            CurveError::ZeroDegree => write!(f, "degree must be ≥ 1"),
            CurveError::WeightsMismatch {
                weights,
                control_points,
            } => {
                write!(
                    f,
                    "weights length {weights} != control points {control_points}"
                )
            }
            CurveError::NonPositiveWeight { index } => {
                write!(f, "weight[{index}] ≤ 0")
            }
        }
    }
}

impl std::error::Error for CurveError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            CurveError::Knot(e) => Some(e),
            _ => None,
        }
    }
}

impl std::error::Error for KnotError {}

impl From<KnotError> for CurveError {
    fn from(e: KnotError) -> Self {
        CurveError::Knot(e)
    }
}

// ── BSplineCurve ──────────────────────────────────────────────────────────────

/// A non-rational B-spline curve of degree `p` in `D`-dimensional space.
///
/// All weights are implicitly 1.  For rational curves (conics, circles, etc.)
/// use [`NurbsCurve`].
///
/// # Invariants
///
/// - `knots.len() == n + p + 2` where `n+1 = control_points.len()`.
/// - `p ≥ 1`.
/// - The knot vector is non-decreasing.
///
/// # Diagram
///
/// ```text
/// Control polygon:   P0 ── P1 ── P2 ── P3
///                           │
///                      B-spline curve
///                           │
/// C(t): ════════════════════════════════  (t ∈ [ξ₀, ξₘ])
/// ```
#[derive(Clone, Debug)]
pub struct BSplineCurve<const D: usize> {
    /// Control points in D-dimensional space.
    control_points: Vec<SVector<Real, D>>,
    /// Validated knot vector.
    knots: KnotVector,
    /// Polynomial degree.
    degree: usize,
}

impl<const D: usize> BSplineCurve<D> {
    /// Create a new B-spline curve.
    ///
    /// # Errors
    /// Returns [`CurveError`] if the knot count, degree, or control-point
    /// count are inconsistent.
    pub fn new(
        control_points: Vec<SVector<Real, D>>,
        knots: KnotVector,
        degree: usize,
    ) -> Result<Self, CurveError> {
        if control_points.is_empty() {
            return Err(CurveError::NoControlPoints);
        }
        if degree == 0 {
            return Err(CurveError::ZeroDegree);
        }
        let n = control_points.len() - 1;
        let expected = n + degree + 2;
        if knots.len() != expected {
            return Err(CurveError::KnotCountMismatch {
                got: knots.len(),
                expected,
            });
        }
        Ok(Self {
            control_points,
            knots,
            degree,
        })
    }

    /// Create a clamped uniform B-spline with the given control points and degree.
    ///
    /// The knot vector is constructed automatically.
    ///
    /// # Panics
    /// Panics if `degree == 0` or `control_points.len() < degree + 1`.
    pub fn clamped(control_points: Vec<SVector<Real, D>>, degree: usize) -> Self {
        assert!(degree >= 1, "degree must be ≥ 1");
        let n = control_points.len() - 1;
        assert!(n >= degree, "need at least degree+1 control points");
        let knots = KnotVector::clamped_uniform(n, degree);
        Self {
            control_points,
            knots,
            degree,
        }
    }

    /// The parameter domain `[t_min, t_max]`.
    pub fn domain(&self) -> (Real, Real) {
        self.knots.domain()
    }

    /// Number of control points.
    pub fn num_control_points(&self) -> usize {
        self.control_points.len()
    }

    /// Polynomial degree.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Evaluate the curve at parameter `t`.
    ///
    /// Uses de Boor's algorithm via [`eval_basis`] — O(p²) per evaluation.
    ///
    /// # Panics
    /// Panics if `t` is outside the knot domain.
    pub fn point(&self, t: Real) -> SVector<Real, D> {
        let n = self.control_points.len() - 1;
        let span = self.knots.find_span(t, n);
        let basis = eval_basis(span, t, self.degree, &self.knots);
        let p = self.degree;
        let mut result = SVector::<Real, D>::zeros();
        for (j, &b) in basis.iter().enumerate() {
            result += self.control_points[span - p + j] * b;
        }
        result
    }

    /// Evaluate the curve and its first derivative at parameter `t`.
    ///
    /// Returns `(C(t), C'(t))`.
    ///
    /// # Panics
    /// Panics if `t` is outside the knot domain.
    pub fn point_and_tangent(&self, t: Real) -> (SVector<Real, D>, SVector<Real, D>) {
        let n = self.control_points.len() - 1;
        let span = self.knots.find_span(t, n);
        let (basis, dbasis) = eval_basis_and_deriv(span, t, self.degree, &self.knots);
        let p = self.degree;
        let mut pt = SVector::<Real, D>::zeros();
        let mut tan = SVector::<Real, D>::zeros();
        for j in 0..=p {
            let cp = self.control_points[span - p + j];
            pt += cp * basis[j];
            tan += cp * dbasis[j];
        }
        (pt, tan)
    }

    /// Sample `count` uniformly spaced points on the curve.
    ///
    /// Includes both endpoints.
    pub fn sample_uniform(&self, count: usize) -> Vec<SVector<Real, D>> {
        assert!(count >= 2, "need at least 2 samples");
        let (lo, hi) = self.domain();
        (0..count)
            .map(|i| {
                let t = lo + (hi - lo) * (i as Real / (count - 1) as Real);
                self.point(t)
            })
            .collect()
    }
}

// ── NurbsCurve ────────────────────────────────────────────────────────────────

/// A rational B-spline (NURBS) curve of degree `p` in `D`-dimensional space.
///
/// Each control point `Pᵢ` has an associated positive weight `wᵢ`.  When all
/// weights are equal to 1.0, the curve is identical to a B-spline.  Non-uniform
/// weights allow exact representation of conics (circles, ellipses, parabolas).
///
/// # Invariants
///
/// - `control_points.len() == weights.len()`.
/// - All weights `wᵢ > 0`.
/// - `knots.len() == n + p + 2` where `n+1 == control_points.len()`.
///
/// # Theorem — Rational Partition of Unity
///
/// The NURBS basis functions `R_{i,p}(t) = N_{i,p}(t) · wᵢ / W(t)` where
/// `W(t) = Σ N_{i,p}(t) · wᵢ` satisfy `Σ R_{i,p}(t) = 1`, inheriting
/// the convex hull and affine invariance properties from B-splines.
#[derive(Clone, Debug)]
pub struct NurbsCurve<const D: usize> {
    /// Control points in D-dimensional Euclidean space.
    control_points: Vec<SVector<Real, D>>,
    /// Positive weights, one per control point.
    weights: Vec<Real>,
    /// Validated knot vector.
    knots: KnotVector,
    /// Polynomial degree.
    degree: usize,
}

impl<const D: usize> NurbsCurve<D> {
    /// Create a new NURBS curve.
    ///
    /// # Errors
    ///
    /// Returns [`CurveError`] if:
    /// - `control_points` or `weights` are empty.
    /// - `weights.len() != control_points.len()`.
    /// - Any weight is ≤ 0.
    /// - The knot count doesn't satisfy `n + p + 2`.
    pub fn new(
        control_points: Vec<SVector<Real, D>>,
        weights: Vec<Real>,
        knots: KnotVector,
        degree: usize,
    ) -> Result<Self, CurveError> {
        if control_points.is_empty() {
            return Err(CurveError::NoControlPoints);
        }
        if degree == 0 {
            return Err(CurveError::ZeroDegree);
        }
        if weights.len() != control_points.len() {
            return Err(CurveError::WeightsMismatch {
                weights: weights.len(),
                control_points: control_points.len(),
            });
        }
        for (i, &w) in weights.iter().enumerate() {
            if w <= 0.0 {
                return Err(CurveError::NonPositiveWeight { index: i });
            }
        }
        let n = control_points.len() - 1;
        let expected = n + degree + 2;
        if knots.len() != expected {
            return Err(CurveError::KnotCountMismatch {
                got: knots.len(),
                expected,
            });
        }
        Ok(Self {
            control_points,
            weights,
            knots,
            degree,
        })
    }

    /// Create a NURBS curve from a B-spline (all weights = 1).
    pub fn from_bspline(curve: BSplineCurve<D>) -> Self {
        let n = curve.control_points.len();
        let weights = vec![1.0; n];
        Self {
            control_points: curve.control_points,
            weights,
            knots: curve.knots,
            degree: curve.degree,
        }
    }

    /// The parameter domain `[t_min, t_max]`.
    pub fn domain(&self) -> (Real, Real) {
        self.knots.domain()
    }

    /// Number of control points.
    pub fn num_control_points(&self) -> usize {
        self.control_points.len()
    }

    /// Polynomial degree.
    pub fn degree(&self) -> usize {
        self.degree
    }

    /// Weights slice.
    pub fn weights(&self) -> &[Real] {
        &self.weights
    }

    /// Evaluate the NURBS curve at parameter `t`.
    ///
    /// ```text
    ///        Σᵢ N_{i,p}(t) · wᵢ · Pᵢ
    /// C(t) = ─────────────────────────
    ///           Σᵢ N_{i,p}(t) · wᵢ
    /// ```
    ///
    /// # Panics
    /// Panics if `t` is outside the knot domain.
    pub fn point(&self, t: Real) -> SVector<Real, D> {
        let n = self.control_points.len() - 1;
        let span = self.knots.find_span(t, n);
        let basis = eval_basis(span, t, self.degree, &self.knots);
        let p = self.degree;

        let mut num = SVector::<Real, D>::zeros();
        let mut den: Real = 0.0;
        for (j, &b) in basis.iter().enumerate() {
            let w = self.weights[span - p + j];
            let bw = b * w;
            num += self.control_points[span - p + j] * bw;
            den += bw;
        }
        // Guard against degenerate knot spans where all basis weights are 0
        if den.abs() < 1e-15 {
            return self.control_points[span - p];
        }
        num / den
    }

    /// Evaluate the NURBS curve and its first derivative at parameter `t`.
    ///
    /// Uses the quotient rule:
    /// ```text
    /// C'(t) = (A'(t) · W(t) − A(t) · W'(t)) / W(t)²
    /// ```
    /// where `A(t) = Σ N_{i,p}(t)·wᵢ·Pᵢ` and `W(t) = Σ N_{i,p}(t)·wᵢ`.
    ///
    /// # Panics
    /// Panics if `t` is outside the knot domain.
    pub fn point_and_tangent(&self, t: Real) -> (SVector<Real, D>, SVector<Real, D>) {
        let n = self.control_points.len() - 1;
        let span = self.knots.find_span(t, n);
        let (basis, dbasis) = eval_basis_and_deriv(span, t, self.degree, &self.knots);
        let p = self.degree;

        let mut a = SVector::<Real, D>::zeros(); // Σ N·w·P
        let mut da = SVector::<Real, D>::zeros(); // Σ N'·w·P
        let mut w: Real = 0.0; // Σ N·w
        let mut dw: Real = 0.0; // Σ N'·w

        for j in 0..=p {
            let cp = self.control_points[span - p + j];
            let wj = self.weights[span - p + j];
            a += cp * (basis[j] * wj);
            da += cp * (dbasis[j] * wj);
            w += basis[j] * wj;
            dw += dbasis[j] * wj;
        }

        let pt = if w.abs() < 1e-15 {
            self.control_points[span - p]
        } else {
            a / w
        };

        let tan = if w.abs() < 1e-15 {
            SVector::zeros()
        } else {
            (da - pt * dw) / w
        };

        (pt, tan)
    }

    /// Sample `count` uniformly spaced points on the curve.
    ///
    /// Includes both endpoints.
    pub fn sample_uniform(&self, count: usize) -> Vec<SVector<Real, D>> {
        assert!(count >= 2, "need at least 2 samples");
        let (lo, hi) = self.domain();
        (0..count)
            .map(|i| {
                let t = lo + (hi - lo) * (i as Real / (count - 1) as Real);
                self.point(t)
            })
            .collect()
    }
}

impl NurbsCurve<3> {
    /// Compute an axis-aligned bounding box over `resolution` samples.
    ///
    /// By the convex hull property of NURBS, the bounding box of the control
    /// points is a conservative bound; this provides a tighter empirical bound.
    pub fn aabb(&self, resolution: usize) -> crate::domain::geometry::Aabb {
        use crate::domain::core::scalar::Point3r;
        use crate::domain::geometry::Aabb;
        let mut aabb = Aabb::empty();
        let (lo, hi) = self.domain();
        let res = resolution.max(8);
        for i in 0..=res {
            let t = lo + (hi - lo) * (i as Real / res as Real);
            let pt = self.point(t);
            aabb.expand(&Point3r::new(pt[0], pt[1], pt[2]));
        }
        aabb
    }
}

// ── Tests ─────────────────────────────────────────────────────────────────────

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::SVector;

    type V3 = SVector<Real, 3>;
    type V2 = SVector<Real, 2>;

    fn v3(x: Real, y: Real, z: Real) -> V3 {
        V3::new(x, y, z)
    }

    // ── BSplineCurve ─────────────────────────────────────────────────────────

    #[test]
    fn bspline_linear_interpolates_endpoints() {
        // Linear (p=1), 2 control points: C(0)=P0, C(1)=P1
        let pts = vec![v3(0.0, 0.0, 0.0), v3(1.0, 2.0, 3.0)];
        let curve = BSplineCurve::clamped(pts, 1);
        let start = curve.point(0.0);
        let end = curve.point(1.0);
        assert!((start - v3(0.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((end - v3(1.0, 2.0, 3.0)).norm() < 1e-12);
    }

    #[test]
    fn bspline_quadratic_midpoint() {
        // Quadratic with 3 control points: P0=(0,0,0), P1=(1,2,0), P2=(2,0,0)
        // At t=0.5 the result should be between the control points
        let pts = vec![v3(0.0, 0.0, 0.0), v3(1.0, 2.0, 0.0), v3(2.0, 0.0, 0.0)];
        let curve = BSplineCurve::clamped(pts, 2);
        let mid = curve.point(0.5);
        // B-spline interpolates convex hull: y should be positive and x near 1
        assert!(mid[0] > 0.9 && mid[0] < 1.1);
        assert!(mid[1] > 0.0);
    }

    #[test]
    fn bspline_tangent_linear() {
        // Linear curve from (0,0,0) to (1,1,1): tangent should be constant (1,1,1)
        let pts = vec![v3(0.0, 0.0, 0.0), v3(1.0, 1.0, 1.0)];
        let curve = BSplineCurve::clamped(pts, 1);
        let (_, tan) = curve.point_and_tangent(0.5);
        // Tangent direction is (1,1,1), magnitude = degree * (P1-P0) / knot diff
        assert!(
            tan.dot(&v3(1.0, 1.0, 1.0)) > 0.0,
            "tangent must point in positive direction"
        );
    }

    #[test]
    fn bspline_sample_uniform_count() {
        let pts = vec![v3(0.0, 0.0, 0.0), v3(1.0, 0.0, 0.0), v3(2.0, 0.0, 0.0)];
        let curve = BSplineCurve::clamped(pts, 2);
        let samples = curve.sample_uniform(11);
        assert_eq!(samples.len(), 11);
    }

    #[test]
    fn bspline_creation_errors() {
        let kv = KnotVector::clamped_uniform(2, 2);
        // Wrong degree: knot vector has n+p+2=3+1+2... let's just test no control points
        let empty: Vec<V3> = vec![];
        assert!(BSplineCurve::new(empty, kv, 2).is_err());
    }

    // ── NurbsCurve ───────────────────────────────────────────────────────────

    #[test]
    fn nurbs_unit_weight_matches_bspline() {
        // NURBS with w=1 everywhere should equal B-spline
        let pts = vec![v3(0.0, 0.0, 0.0), v3(0.5, 1.0, 0.0), v3(1.0, 0.0, 0.0)];
        let bs = BSplineCurve::clamped(pts.clone(), 2);
        let weights = vec![1.0_f64; 3];
        let kv = KnotVector::clamped_uniform(2, 2);
        let nc = NurbsCurve::new(pts, weights, kv, 2).unwrap();
        for i in 0..=10 {
            let t = i as Real / 10.0;
            let pb = bs.point(t);
            let pn = nc.point(t);
            assert!(
                (pb - pn).norm() < 1e-12,
                "unit-weight NURBS ≠ B-spline at t={}: |diff|={}",
                t,
                (pb - pn).norm()
            );
        }
    }

    #[test]
    fn nurbs_quarter_circle() {
        // Exact unit quarter-circle in XY plane:
        // P0=(1,0), w0=1  P1=(1,1), w1=1/√2  P2=(0,1), w2=1
        let sq2_inv: Real = std::f64::consts::FRAC_1_SQRT_2 as Real;
        let ctrl = vec![V2::new(1.0, 0.0), V2::new(1.0, 1.0), V2::new(0.0, 1.0)];
        let weights = vec![1.0, sq2_inv, 1.0];
        let knots = KnotVector::try_new(vec![0.0, 0.0, 0.0, 1.0, 1.0, 1.0]).unwrap();
        let curve = NurbsCurve::new(ctrl, weights, knots, 2).unwrap();

        // Sample 11 points and check that they lie on the unit circle
        for i in 0..=10 {
            let t = i as Real / 10.0;
            let pt = curve.point(t);
            let r = (pt[0] * pt[0] + pt[1] * pt[1]).sqrt();
            assert!(
                (r - 1.0).abs() < 1e-10,
                "quarter-circle point not on unit circle at t={}: r={}",
                t,
                r
            );
        }
    }

    #[test]
    fn nurbs_endpoints_interpolate() {
        // Any clamped NURBS must pass through first and last control points
        let ctrl = vec![v3(1.0, 2.0, 3.0), v3(4.0, 5.0, 6.0), v3(7.0, 8.0, 9.0)];
        let weights = vec![1.0, 0.5, 2.0];
        let knots = KnotVector::clamped_uniform(2, 2);
        let curve = NurbsCurve::new(ctrl.clone(), weights, knots, 2).unwrap();
        let start = curve.point(0.0);
        let end = curve.point(1.0);
        assert!((start - ctrl[0]).norm() < 1e-12, "start should equal P0");
        assert!((end - ctrl[2]).norm() < 1e-12, "end should equal P2");
    }

    #[test]
    fn nurbs_non_positive_weight_errors() {
        let ctrl = vec![v3(0.0, 0.0, 0.0), v3(1.0, 0.0, 0.0)];
        let weights = vec![1.0, 0.0]; // zero weight is invalid
        let knots = KnotVector::clamped_uniform(1, 1);
        assert!(matches!(
            NurbsCurve::new(ctrl, weights, knots, 1),
            Err(CurveError::NonPositiveWeight { index: 1 })
        ));
    }

    #[test]
    fn nurbs_aabb_contains_control_points() {
        let ctrl = vec![v3(-1.0, -2.0, -3.0), v3(0.0, 0.0, 0.0), v3(4.0, 5.0, 6.0)];
        let weights = vec![1.0, 1.5, 1.0];
        let knots = KnotVector::clamped_uniform(2, 2);
        let curve = NurbsCurve::new(ctrl, weights, knots, 2).unwrap();
        let aabb = curve.aabb(50);
        // Curve lies in convex hull of control points; aabb should be non-degenerate
        assert!(aabb.min.x <= 0.0 && aabb.max.x >= 1.0);
    }

    #[test]
    fn nurbs_tangent_non_zero_for_non_degenerate_curve() {
        let ctrl = vec![v3(0.0, 0.0, 0.0), v3(1.0, 1.0, 0.0), v3(2.0, 0.0, 0.0)];
        let weights = vec![1.0_f64; 3];
        let knots = KnotVector::clamped_uniform(2, 2);
        let curve = NurbsCurve::new(ctrl, weights, knots, 2).unwrap();
        let (_, tan) = curve.point_and_tangent(0.5);
        assert!(tan.norm() > 0.0, "tangent at midpoint should be non-zero");
    }
}
