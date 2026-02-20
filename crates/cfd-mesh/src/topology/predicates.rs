//! Exact geometric predicates for topological validation.
//!
//! Exposes robust exact orientation tests to prevent floating-point
//! heuristics from causing degenerate topological failures like non-manifold
//! edge creation. These functions wrap adaptive multi-precision arithmetic.

use crate::core::scalar::Point3r;

/// Exact algebraic sign representing geometric orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, PartialOrd, Ord)]
pub enum Sign {
    Negative = -1,
    Zero = 0,
    Positive = 1,
}

impl Sign {
    /// Convert the exact expansion floating-point result into a strict sign.
    #[inline]
    pub fn from_exact_f64(v: f64) -> Self {
        if v > 0.0 {
            Sign::Positive
        } else if v < 0.0 {
            Sign::Negative
        } else {
            Sign::Zero
        }
    }
    
    pub fn is_positive(self) -> bool { self == Sign::Positive }
    pub fn is_negative(self) -> bool { self == Sign::Negative }
    pub fn is_zero(self) -> bool { self == Sign::Zero }
}

/// Exact 3D orientation predicate.
/// 
/// Returns whether the point `d` is strictly above, strictly below, or perfectly
/// coplanar with the oriented plane defined by `a`, `b`, and `c`.
/// 
/// This evaluation is mathematically exact and immune to floating-point epsilon noise.
#[inline]
pub fn orient3d(a: &Point3r, b: &Point3r, c: &Point3r, d: &Point3r) -> Sign {
    let pa = [a.x as f64, a.y as f64, a.z as f64];
    let pb = [b.x as f64, b.y as f64, b.z as f64];
    let pc = [c.x as f64, c.y as f64, c.z as f64];
    let pd = [d.x as f64, d.y as f64, d.z as f64];
    
    let det = geometry_predicates::orient3d(pa, pb, pc, pd);
    Sign::from_exact_f64(det)
}

/// Exact 2D orientation predicate (Sutherland-Hodgman / coplanar clipping).
/// 
/// Returns whether the point `c` lies strictly left, strictly right, or perfectly
/// collinear with the directed line from `a` to `b` in the 2D plane (X-Y).
#[inline]
pub fn orient2d(a: &Point3r, b: &Point3r, c: &Point3r) -> Sign {
    let pa = [a.x as f64, a.y as f64];
    let pb = [b.x as f64, b.y as f64];
    let pc = [c.x as f64, c.y as f64];
    
    let det = geometry_predicates::orient2d(pa, pb, pc);
    Sign::from_exact_f64(det)
}

/// Exact incircle predicate in 2D.
#[inline]
pub fn incircle2d(a: &Point3r, b: &Point3r, c: &Point3r, d: &Point3r) -> Sign {
    let pa = [a.x as f64, a.y as f64];
    let pb = [b.x as f64, b.y as f64];
    let pc = [c.x as f64, c.y as f64];
    let pd = [d.x as f64, d.y as f64];
    
    let det = geometry_predicates::incircle(pa, pb, pc, pd);
    Sign::from_exact_f64(det)
}

/// Exact insphere predicate in 3D.
#[inline]
pub fn insphere3d(a: &Point3r, b: &Point3r, c: &Point3r, d: &Point3r, e: &Point3r) -> Sign {
    let pa = [a.x as f64, a.y as f64, a.z as f64];
    let pb = [b.x as f64, b.y as f64, b.z as f64];
    let pc = [c.x as f64, c.y as f64, c.z as f64];
    let pd = [d.x as f64, d.y as f64, d.z as f64];
    let pe = [e.x as f64, e.y as f64, e.z as f64];
    
    let det = geometry_predicates::insphere(pa, pb, pc, pd, pe);
    Sign::from_exact_f64(det)
}
