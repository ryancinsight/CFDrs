//! NURBS and B-Spline Surfaces -- tensor-product parameterisation.
//!
//! `BSplineSurface` is the non-rational case (all weights == 1).
//! `NurbsSurface` is the rational case with per-control-point positive weights.
//! Both are parameterised over a rectangular domain [u0,u1] x [v0,v1].

use crate::core::scalar::Real;
use super::basis::{eval_basis, eval_basis_and_deriv};
use super::knot::KnotVector;
use nalgebra::{Point3, UnitVector3, Vector3};

// ---------------------------------------------------------------------------
// Type aliases
// ---------------------------------------------------------------------------

/// 3-D point (matches `crate::core::scalar::Point3r`).
pub type Pt3 = Point3<Real>;

/// 3-D vector.
pub type Vec3 = Vector3<Real>;

// ---------------------------------------------------------------------------
// ControlGrid
// ---------------------------------------------------------------------------

/// A rectangular grid of 3-D control points stored row-major.
///
/// `grid.get(i, j)` returns the control point at row `i` (u-direction),
/// column `j` (v-direction).
#[derive(Clone, Debug)]
pub struct ControlGrid {
    data: Vec<Pt3>,
    n_rows: usize,
    n_cols: usize,
}

impl ControlGrid {
    /// Create from a flat row-major vector.
    ///
    /// # Panics
    /// Panics if `data.len() != n_rows * n_cols`.
    pub fn new(data: Vec<Pt3>, n_rows: usize, n_cols: usize) -> Self {
        assert_eq!(data.len(), n_rows * n_cols, "data.len() must equal n_rows * n_cols");
        Self { data, n_rows, n_cols }
    }

    /// Number of rows (u direction).
    pub fn n_rows(&self) -> usize { self.n_rows }

    /// Number of columns (v direction).
    pub fn n_cols(&self) -> usize { self.n_cols }

    /// Access control point at `(row, col)`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Pt3 {
        self.data[i * self.n_cols + j]
    }
}

// ---------------------------------------------------------------------------
// WeightGrid
// ---------------------------------------------------------------------------

/// A rectangular grid of positive NURBS weights, stored row-major.
#[derive(Clone, Debug)]
pub struct WeightGrid {
    data: Vec<Real>,
    n_rows: usize,
    n_cols: usize,
}

impl WeightGrid {
    /// All-ones weight grid (equivalent to B-spline).
    pub fn uniform(n_rows: usize, n_cols: usize) -> Self {
        Self { data: vec![1.0; n_rows * n_cols], n_rows, n_cols }
    }

    /// Create from a flat row-major vector.
    ///
    /// # Panics
    /// Panics if any weight <= 0 or length mismatch.
    pub fn new(data: Vec<Real>, n_rows: usize, n_cols: usize) -> Self {
        assert_eq!(data.len(), n_rows * n_cols);
        for (i, &w) in data.iter().enumerate() {
            assert!(w > 0.0, "weight[{}] = {} is not positive", i, w);
        }
        Self { data, n_rows, n_cols }
    }

    /// Access weight at `(row, col)`.
    #[inline]
    pub fn get(&self, i: usize, j: usize) -> Real {
        self.data[i * self.n_cols + j]
    }

    /// Dimensions.
    pub fn n_rows(&self) -> usize { self.n_rows }
    /// Dimensions.
    pub fn n_cols(&self) -> usize { self.n_cols }
}

// ---------------------------------------------------------------------------
// SurfaceError
// ---------------------------------------------------------------------------

/// Error returned when constructing a B-spline or NURBS surface.
#[derive(Clone, Debug, PartialEq)]
pub enum SurfaceError {
    /// Control grid is empty.
    EmptyControlGrid,
    /// Degree in the named direction is zero.
    ZeroDegree {
        /// Direction character: `'u'` or `'v'`.
        direction: char,
    },
    /// Knot count does not satisfy `n + p + 2`.
    KnotCountMismatch {
        /// Direction character: `'u'` or `'v'`.
        direction: char,
        /// Actual number of knots supplied.
        got: usize,
        /// Required number of knots (`n + degree + 2`).
        expected: usize,
    },
    /// Weight grid dimensions differ from the control grid.
    WeightGridMismatch,
}

impl std::fmt::Display for SurfaceError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SurfaceError::EmptyControlGrid => write!(f, "control grid is empty"),
            SurfaceError::ZeroDegree { direction } => {
                write!(f, "degree in {} direction must be >= 1", direction)
            }
            SurfaceError::KnotCountMismatch { direction, got, expected } => {
                write!(f, "knot-{}: expected {} knots, got {}", direction, expected, got)
            }
            SurfaceError::WeightGridMismatch => {
                write!(f, "weight grid dimensions differ from control grid")
            }
        }
    }
}

impl std::error::Error for SurfaceError {}

// ---------------------------------------------------------------------------
// BSplineSurface
// ---------------------------------------------------------------------------

/// A non-rational tensor-product B-spline surface.
///
/// Evaluation:
///   S(u, v) = sum_i sum_j  N_{i,p}(u) * N_{j,q}(v) * P_{ij}
///
/// where N_{i,p} and N_{j,q} are B-spline basis functions computed by
/// Cox-de Boor recursion.
#[derive(Clone, Debug)]
pub struct BSplineSurface {
    /// Control point grid, row-major in (u, v).
    pub control_grid: ControlGrid,
    /// Knot vector in the u direction.
    pub knots_u: KnotVector,
    /// Knot vector in the v direction.
    pub knots_v: KnotVector,
    /// Degree in the u direction.
    pub degree_u: usize,
    /// Degree in the v direction.
    pub degree_v: usize,
}

impl BSplineSurface {
    /// Create a B-spline surface, validating knot / control-point consistency.
    pub fn new(
        control_grid: ControlGrid,
        knots_u: KnotVector,
        knots_v: KnotVector,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self, SurfaceError> {
        if control_grid.n_rows() == 0 || control_grid.n_cols() == 0 {
            return Err(SurfaceError::EmptyControlGrid);
        }
        if degree_u == 0 { return Err(SurfaceError::ZeroDegree { direction: 'u' }); }
        if degree_v == 0 { return Err(SurfaceError::ZeroDegree { direction: 'v' }); }

        // n + p + 2  where  n_u = n_cols - 1,  n_v = n_rows - 1
        // (u varies along columns, v varies along rows)
        let exp_u = control_grid.n_cols() + degree_u + 1;
        if knots_u.len() != exp_u {
            return Err(SurfaceError::KnotCountMismatch {
                direction: 'u', got: knots_u.len(), expected: exp_u,
            });
        }
        let exp_v = control_grid.n_rows() + degree_v + 1;
        if knots_v.len() != exp_v {
            return Err(SurfaceError::KnotCountMismatch {
                direction: 'v', got: knots_v.len(), expected: exp_v,
            });
        }
        Ok(Self { control_grid, knots_u, knots_v, degree_u, degree_v })
    }

    /// Create with clamped uniform knot vectors constructed automatically.
    pub fn clamped(
        control_grid: ControlGrid,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self, SurfaceError> {
        if control_grid.n_rows() == 0 || control_grid.n_cols() == 0 {
            return Err(SurfaceError::EmptyControlGrid);
        }
        let n_u = control_grid.n_cols() - 1;  // u = cols
        let n_v = control_grid.n_rows() - 1;  // v = rows
        let ku = KnotVector::clamped_uniform(n_u, degree_u);
        let kv = KnotVector::clamped_uniform(n_v, degree_v);
        Self::new(control_grid, ku, kv, degree_u, degree_v)
    }

    /// Parameter domain `((u_min, u_max), (v_min, v_max))`.
    pub fn domain(&self) -> ((Real, Real), (Real, Real)) {
        (self.knots_u.domain(), self.knots_v.domain())
    }

    /// Evaluate the surface at `(u, v)`.
    pub fn point(&self, u: Real, v: Real) -> Pt3 {
        let n_u = self.control_grid.n_cols() - 1;  // u = cols
        let n_v = self.control_grid.n_rows() - 1;  // v = rows
        let su = self.knots_u.find_span(u, n_u);
        let sv = self.knots_v.find_span(v, n_v);
        let bu = eval_basis(su, u, self.degree_u, &self.knots_u);
        let bv = eval_basis(sv, v, self.degree_v, &self.knots_v);

        let pu = self.degree_u;
        let pv = self.degree_v;
        let mut res = Vec3::zeros();
        for (j, &nu) in bu.iter().enumerate() {
            for (k, &nv) in bv.iter().enumerate() {
                // get(row, col) = get(v_idx, u_idx)
                res += self.control_grid.get(sv - pv + k, su - pu + j).coords * (nu * nv);
            }
        }
        Pt3::from(res)
    }

    /// Evaluate surface point and partial derivatives `(S, dS/du, dS/dv)`.
    pub fn point_and_derivs(&self, u: Real, v: Real) -> (Pt3, Vec3, Vec3) {
        let n_u = self.control_grid.n_cols() - 1;  // u = cols
        let n_v = self.control_grid.n_rows() - 1;  // v = rows
        let su = self.knots_u.find_span(u, n_u);
        let sv = self.knots_v.find_span(v, n_v);
        let (bu, dbu) = eval_basis_and_deriv(su, u, self.degree_u, &self.knots_u);
        let (bv, dbv) = eval_basis_and_deriv(sv, v, self.degree_v, &self.knots_v);

        let pu = self.degree_u;
        let pv = self.degree_v;
        let mut s = Vec3::zeros();
        let mut ds_du = Vec3::zeros();
        let mut ds_dv = Vec3::zeros();

        for (j, (&nu, &dnu)) in bu.iter().zip(dbu.iter()).enumerate() {
            for (k, (&nv, &dnv)) in bv.iter().zip(dbv.iter()).enumerate() {
                let pt = self.control_grid.get(sv - pv + k, su - pu + j).coords;
                s      += pt * (nu * nv);
                ds_du  += pt * (dnu * nv);
                ds_dv  += pt * (nu * dnv);
            }
        }
        (Pt3::from(s), ds_du, ds_dv)
    }

    /// Unit surface normal at `(u, v)` = normalize(dS/du cross dS/dv).
    /// Returns `None` if the surface is degenerate at `(u, v)`.
    pub fn normal(&self, u: Real, v: Real) -> Option<UnitVector3<Real>> {
        let (_, du, dv) = self.point_and_derivs(u, v);
        UnitVector3::try_new(du.cross(&dv), 1e-15)
    }
}

// ---------------------------------------------------------------------------
// NurbsSurface
// ---------------------------------------------------------------------------

/// A rational tensor-product NURBS surface.
///
/// Evaluation (rational):
///   S(u,v) = (sum_ij N_i(u) * N_j(v) * w_ij * P_ij)
///           / (sum_ij N_i(u) * N_j(v) * w_ij)
///
/// Exact conics (spheres, cylinders) arise from specific weight configurations.
#[derive(Clone, Debug)]
pub struct NurbsSurface {
    /// Control point grid.
    pub control_grid: ControlGrid,
    /// Positive weight at each control point.
    pub weights: WeightGrid,
    /// Knot vector in the u direction.
    pub knots_u: KnotVector,
    /// Knot vector in the v direction.
    pub knots_v: KnotVector,
    /// Degree in the u direction.
    pub degree_u: usize,
    /// Degree in the v direction.
    pub degree_v: usize,
}

impl NurbsSurface {
    /// Create a NURBS surface, validating all dimensions.
    pub fn new(
        control_grid: ControlGrid,
        weights: WeightGrid,
        knots_u: KnotVector,
        knots_v: KnotVector,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self, SurfaceError> {
        // Delegate structural checks to BSplineSurface constructor
        let _test = BSplineSurface::new(
            control_grid.clone(), knots_u.clone(), knots_v.clone(), degree_u, degree_v,
        )?;
        if weights.n_rows() != control_grid.n_rows()
            || weights.n_cols() != control_grid.n_cols()
        {
            return Err(SurfaceError::WeightGridMismatch);
        }
        Ok(Self { control_grid, weights, knots_u, knots_v, degree_u, degree_v })
    }

    /// Create a NURBS surface from a B-spline (all weights = 1).
    pub fn from_bspline(s: BSplineSurface) -> Self {
        let w = WeightGrid::uniform(s.control_grid.n_rows(), s.control_grid.n_cols());
        Self {
            control_grid: s.control_grid,
            weights: w,
            knots_u: s.knots_u,
            knots_v: s.knots_v,
            degree_u: s.degree_u,
            degree_v: s.degree_v,
        }
    }

    /// Create with automatic clamped uniform knot vectors and uniform weights.
    pub fn clamped(
        control_grid: ControlGrid,
        degree_u: usize,
        degree_v: usize,
    ) -> Result<Self, SurfaceError> {
        let s = BSplineSurface::clamped(control_grid, degree_u, degree_v)?;
        Ok(Self::from_bspline(s))
    }

    /// Parameter domain `((u_min, u_max), (v_min, v_max))`.
    pub fn domain(&self) -> ((Real, Real), (Real, Real)) {
        (self.knots_u.domain(), self.knots_v.domain())
    }

    /// Evaluate the NURBS surface at `(u, v)`.
    pub fn point(&self, u: Real, v: Real) -> Pt3 {
        let (num, den) = self.rational_eval(u, v);
        if den.abs() < 1e-15 {
            return self.control_grid.get(0, 0);
        }
        Pt3::from(num / den)
    }

    /// Evaluate surface point and partial derivatives `(S, dS/du, dS/dv)`.
    ///
    /// Uses the quotient rule:
    ///   dS/du = (dA/du * W - A * dW/du) / W^2
    /// where A = sum N_i(u) N_j(v) w_ij P_ij and W = sum N_i N_j w_ij.
    pub fn point_and_derivs(&self, u: Real, v: Real) -> (Pt3, Vec3, Vec3) {
        let n_u = self.control_grid.n_cols() - 1;  // u = cols
        let n_v = self.control_grid.n_rows() - 1;  // v = rows
        let su = self.knots_u.find_span(u, n_u);
        let sv = self.knots_v.find_span(v, n_v);
        let (bu, dbu) = eval_basis_and_deriv(su, u, self.degree_u, &self.knots_u);
        let (bv, dbv) = eval_basis_and_deriv(sv, v, self.degree_v, &self.knots_v);

        let pu = self.degree_u;
        let pv = self.degree_v;

        let mut a = Vec3::zeros();
        let mut da_du = Vec3::zeros();
        let mut da_dv = Vec3::zeros();
        let mut w: Real = 0.0;
        let mut dw_du: Real = 0.0;
        let mut dw_dv: Real = 0.0;

        for (j, (&nu, &dnu)) in bu.iter().zip(dbu.iter()).enumerate() {
            for (k, (&nv, &dnv)) in bv.iter().zip(dbv.iter()).enumerate() {
                // get(row, col) = get(v_idx, u_idx)
                let wij = self.weights.get(sv - pv + k, su - pu + j);
                let pt = self.control_grid.get(sv - pv + k, su - pu + j).coords;
                a      += pt * (nu * nv * wij);
                da_du  += pt * (dnu * nv * wij);
                da_dv  += pt * (nu * dnv * wij);
                w      += nu * nv * wij;
                dw_du  += dnu * nv * wij;
                dw_dv  += nu * dnv * wij;
            }
        }

        if w.abs() < 1e-15 {
            return (self.control_grid.get(0, 0), Vec3::zeros(), Vec3::zeros());
        }
        let s = a / w;
        // Quotient rule: d(a/w)/du = (da/du - s * dw/du) / w
        let ds_du = (da_du - s * dw_du) / w;
        let ds_dv = (da_dv - s * dw_dv) / w;
        (Pt3::from(s), ds_du, ds_dv)
    }

    /// Unit surface normal at `(u, v)`.
    /// Returns `None` if degenerate (zero cross product).
    pub fn normal(&self, u: Real, v: Real) -> Option<UnitVector3<Real>> {
        let (_, du, dv) = self.point_and_derivs(u, v);
        UnitVector3::try_new(du.cross(&dv), 1e-15)
    }

    /// Axis-aligned bounding box from a resolution x resolution sample grid.
    pub fn aabb(&self, resolution: usize) -> crate::geometry::Aabb {
        use crate::geometry::Aabb;
        let mut aabb = Aabb::empty();
        let ((u0, u1), (v0, v1)) = self.domain();
        let res = resolution.max(4);
        for i in 0..=res {
            let u = u0 + (u1 - u0) * (i as Real / res as Real);
            for j in 0..=res {
                let v = v0 + (v1 - v0) * (j as Real / res as Real);
                let pt = self.point(u, v);
                aabb.expand(&pt);
            }
        }
        aabb
    }

    // -- internal --

    fn rational_eval(&self, u: Real, v: Real) -> (Vec3, Real) {
        let n_u = self.control_grid.n_cols() - 1;  // u = cols
        let n_v = self.control_grid.n_rows() - 1;  // v = rows
        let su = self.knots_u.find_span(u, n_u);
        let sv = self.knots_v.find_span(v, n_v);
        let bu = eval_basis(su, u, self.degree_u, &self.knots_u);
        let bv = eval_basis(sv, v, self.degree_v, &self.knots_v);

        let pu = self.degree_u;
        let pv = self.degree_v;
        let mut num = Vec3::zeros();
        let mut den: Real = 0.0;

        for (j, &nu) in bu.iter().enumerate() {
            for (k, &nv) in bv.iter().enumerate() {
                // get(row, col) = get(v_idx, u_idx)
                let wij = self.weights.get(sv - pv + k, su - pu + j);
                let bwij = nu * nv * wij;
                num += self.control_grid.get(sv - pv + k, su - pu + j).coords * bwij;
                den += bwij;
            }
        }
        (num, den)
    }
}

// ---------------------------------------------------------------------------
// Tests
// ---------------------------------------------------------------------------

#[cfg(test)]
mod tests {
    use super::*;

    fn pt(x: Real, y: Real, z: Real) -> Pt3 { Pt3::new(x, y, z) }

    fn grid_2x2(pts: &[Pt3]) -> ControlGrid {
        ControlGrid::new(pts.to_vec(), 2, 2)
    }

    // -- BSplineSurface --

    #[test]
    fn bspline_bilinear_corners() {
        let pts = vec![
            pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0), pt(1.0, 1.0, 0.0),
        ];
        let surf = BSplineSurface::clamped(grid_2x2(&pts), 1, 1).unwrap();
        assert!((surf.point(0.0, 0.0) - pt(0.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((surf.point(1.0, 0.0) - pt(1.0, 0.0, 0.0)).norm() < 1e-12);
        assert!((surf.point(0.0, 1.0) - pt(0.0, 1.0, 0.0)).norm() < 1e-12);
        assert!((surf.point(1.0, 1.0) - pt(1.0, 1.0, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn bspline_bilinear_midpoint() {
        let pts = vec![
            pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0), pt(1.0, 1.0, 0.0),
        ];
        let surf = BSplineSurface::clamped(grid_2x2(&pts), 1, 1).unwrap();
        let mid = surf.point(0.5, 0.5);
        assert!((mid - pt(0.5, 0.5, 0.0)).norm() < 1e-12);
    }

    #[test]
    fn bspline_normal_flat_patch() {
        let pts = vec![
            pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0), pt(1.0, 1.0, 0.0),
        ];
        let surf = BSplineSurface::clamped(grid_2x2(&pts), 1, 1).unwrap();
        let n = surf.normal(0.5, 0.5).expect("should have valid normal");
        assert!((n.dot(&Vec3::new(0.0, 0.0, 1.0)) - 1.0).abs() < 1e-10);
    }

    #[test]
    fn bspline_wrong_knot_count_errors() {
        let pts = vec![
            pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0), pt(1.0, 1.0, 0.0),
        ];
        let grid = grid_2x2(&pts);
        let ku = KnotVector::try_new(vec![0.0, 1.0]).unwrap(); // too short
        let kv = KnotVector::clamped_uniform(1, 1);
        assert!(BSplineSurface::new(grid, ku, kv, 1, 1).is_err());
    }

    // -- NurbsSurface --

    #[test]
    fn nurbs_unit_weight_matches_bspline() {
        let pts = vec![
            pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0), pt(1.0, 1.0, 0.0),
        ];
        let bs = BSplineSurface::clamped(ControlGrid::new(pts.clone(), 2, 2), 1, 1).unwrap();
        let ns = NurbsSurface::from_bspline(bs.clone());
        for i in 0..=5 {
            for j in 0..=5 {
                let u = i as Real / 5.0;
                let v = j as Real / 5.0;
                let pb = bs.point(u, v);
                let pn = ns.point(u, v);
                assert!(
                    (pb - pn).norm() < 1e-12,
                    "unit-weight NURBS != B-spline at ({}, {})", u, v
                );
            }
        }
    }

    #[test]
    fn nurbs_clamped_corners_interpolated() {
        let pts = vec![
            pt(0.0, 0.0, 1.0), pt(1.0, 0.0, 2.0),
            pt(0.0, 1.0, 3.0), pt(1.0, 1.0, 4.0),
        ];
        let surf = NurbsSurface::clamped(ControlGrid::new(pts.clone(), 2, 2), 1, 1).unwrap();
        assert!((surf.point(0.0, 0.0) - pts[0]).norm() < 1e-12);
        assert!((surf.point(1.0, 0.0) - pts[1]).norm() < 1e-12);
        assert!((surf.point(0.0, 1.0) - pts[2]).norm() < 1e-12);
        assert!((surf.point(1.0, 1.0) - pts[3]).norm() < 1e-12);
    }

    #[test]
    fn nurbs_aabb_non_degenerate() {
        let pts = vec![
            pt(-1.0, -1.0, 0.0), pt(1.0, -1.0, 0.0),
            pt(-1.0,  1.0, 0.0), pt(1.0,  1.0, 0.0),
        ];
        let surf = NurbsSurface::clamped(ControlGrid::new(pts, 2, 2), 1, 1).unwrap();
        let aabb = surf.aabb(8);
        assert!(aabb.max.x > 0.0 && aabb.min.x < 0.0);
    }

    #[test]
    fn nurbs_normal_not_zero() {
        let pts = vec![
            pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0), pt(1.0, 1.0, 0.0),
        ];
        let surf = NurbsSurface::clamped(ControlGrid::new(pts, 2, 2), 1, 1).unwrap();
        assert!(surf.normal(0.5, 0.5).is_some());
    }

    #[test]
    fn weight_grid_mismatch_errors() {
        let pts = vec![
            pt(0.0, 0.0, 0.0), pt(1.0, 0.0, 0.0),
            pt(0.0, 1.0, 0.0), pt(1.0, 1.0, 0.0),
        ];
        let grid = ControlGrid::new(pts, 2, 2);
        let wrong_weights = WeightGrid::uniform(3, 3); // 3x3 != 2x2
        let ku = KnotVector::clamped_uniform(1, 1);
        let kv = KnotVector::clamped_uniform(1, 1);
        assert!(NurbsSurface::new(grid, wrong_weights, ku, kv, 1, 1).is_err());
    }
}

