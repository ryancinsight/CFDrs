//! Geometry abstractions for Method of Manufactured Solutions (MMS)
//!
//! This module provides geometry abstractions that allow MMS to work with
//! complex computational domains beyond simple rectangular boxes.
//!
//! References:
//! - Roache, P.J. (2002) "Code Verification by the Method of Manufactured Solutions"

use nalgebra::RealField;

// Submodules
pub mod annular;
pub mod circular;
pub mod rectangular;

// Public API exports
pub use self::annular::AnnularDomain;
pub use self::circular::CircularDomain;
pub use self::rectangular::RectangularDomain;

/// Point in 2D space
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point2D<T: RealField> {
    /// X-coordinate
    pub x: T,
    /// Y-coordinate
    pub y: T,
}

/// Point in 3D space
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct Point3D<T: RealField> {
    /// X-coordinate
    pub x: T,
    /// Y-coordinate
    pub y: T,
    /// Z-coordinate
    pub z: T,
}

/// Boundary condition specification for MMS
#[derive(Debug, Clone, PartialEq)]
pub enum BoundaryCondition<T: RealField> {
    /// Dirichlet boundary condition: φ = g
    Dirichlet(T),
    /// Neumann boundary condition: ∂φ/∂n = g
    Neumann(T),
    /// Robin boundary condition: αφ + β∂φ/∂n = γ
    Robin { alpha: T, beta: T, gamma: T },
}

/// Boundary face identifier
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryFace {
    /// Bottom boundary (y = y_min)
    Bottom,
    /// Top boundary (y = y_max)
    Top,
    /// Left boundary (x = x_min)
    Left,
    /// Right boundary (x = x_max)
    Right,
    /// Front boundary (z = z_min) - for 3D
    Front,
    /// Back boundary (z = z_max) - for 3D
    Back,
    /// Inner boundary (e.g., for annular domains)
    Inner,
    /// Outer boundary (e.g., for annular domains)
    Outer,
}

/// Geometry trait for computational domains used in MMS
pub trait Geometry<T: RealField + Copy> {
    /// Check if a point is inside the computational domain
    fn contains(&self, point: &Point2D<T>) -> bool;

    /// Get the distance to the nearest boundary
    fn distance_to_boundary(&self, point: &Point2D<T>) -> T;

    /// Get the normal vector at a boundary point
    fn boundary_normal(&self, point: &Point2D<T>) -> Option<Point2D<T>>;

    /// Get boundary condition for a specific boundary face
    fn boundary_condition(&self, face: BoundaryFace, s: T) -> BoundaryCondition<T>;

    /// Get domain bounds
    fn bounds(&self) -> (Point2D<T>, Point2D<T>);

    /// Get parametric coordinate along boundary (0 <= s <= 1)
    fn boundary_parameter(&self, face: BoundaryFace, point: &Point2D<T>) -> Option<T>;

    /// Check if a point is on a boundary face
    fn on_boundary(&self, point: &Point2D<T>, face: BoundaryFace, tolerance: T) -> bool;

    /// Get the area/volume of the domain
    fn measure(&self) -> T;
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_rectangular_domain_contains() {
        let domain: RectangularDomain<f64> = RectangularDomain::unit_square();

        assert!(domain.contains(&Point2D { x: 0.5, y: 0.5 }));
        assert!(domain.contains(&Point2D { x: 0.0, y: 0.0 }));
        assert!(domain.contains(&Point2D { x: 1.0, y: 1.0 }));

        assert!(!domain.contains(&Point2D { x: -0.1, y: 0.5 }));
        assert!(!domain.contains(&Point2D { x: 1.1, y: 0.5 }));
        assert!(!domain.contains(&Point2D { x: 0.5, y: -0.1 }));
        assert!(!domain.contains(&Point2D { x: 0.5, y: 1.1 }));
    }

    #[test]
    fn test_circular_domain_contains() {
        let domain: CircularDomain<f64> = CircularDomain::unit_disk();

        // Points inside the unit disk
        assert!(domain.contains(&Point2D { x: 0.0, y: 0.0 })); // center
        assert!(domain.contains(&Point2D { x: 0.5, y: 0.0 })); // on x-axis
        assert!(domain.contains(&Point2D { x: 0.0, y: 0.5 })); // on y-axis

        // Points outside the unit disk
        assert!(!domain.contains(&Point2D { x: 1.1, y: 0.0 })); // outside
        assert!(!domain.contains(&Point2D { x: 0.0, y: 1.1 })); // outside
        assert!(!domain.contains(&Point2D { x: 1.0, y: 1.0 })); // corner, distance = √2 ≈ 1.414
    }

    #[test]
    fn test_circular_domain_measure() {
        let domain: CircularDomain<f64> = CircularDomain::unit_disk();
        let expected_area = std::f64::consts::PI * 1.0 * 1.0; // πr²
        assert_relative_eq!(domain.measure(), expected_area, epsilon = 1e-10);
    }

    #[test]
    fn test_annular_domain_contains() {
        let domain: AnnularDomain<f64> = AnnularDomain::unit_annulus();

        // Points inside the annulus (between r=0.5 and r=1.0)
        assert!(domain.contains(&Point2D { x: 0.6, y: 0.0 })); // on x-axis
        assert!(domain.contains(&Point2D { x: 0.0, y: 0.8 })); // on y-axis
        assert!(domain.contains(&Point2D { x: 0.6, y: 0.8 })); // diagonal

        // Points inside inner circle (should not be contained)
        assert!(!domain.contains(&Point2D { x: 0.3, y: 0.0 })); // inside inner
        assert!(!domain.contains(&Point2D { x: 0.0, y: 0.4 })); // inside inner

        // Points outside outer circle (should not be contained)
        assert!(!domain.contains(&Point2D { x: 1.1, y: 0.0 })); // outside outer
        assert!(!domain.contains(&Point2D { x: 0.0, y: 1.2 })); // outside outer

        // Points on inner boundary (should be contained)
        assert!(domain.contains(&Point2D { x: 0.5, y: 0.0 })); // on inner boundary

        // Points on outer boundary (should be contained)
        assert!(domain.contains(&Point2D { x: 1.0, y: 0.0 })); // on outer boundary
    }

    #[test]
    fn test_annular_domain_distance_to_boundary() {
        let domain: AnnularDomain<f64> = AnnularDomain::unit_annulus();

        // Point inside annulus
        let point_inside = Point2D { x: 0.7, y: 0.0 };
        let dist_inside = domain.distance_to_boundary(&point_inside);
        // Distance to inner boundary: 0.7 - 0.5 = 0.2
        // Distance to outer boundary: 1.0 - 0.7 = 0.3
        // Should return minimum: 0.2
        assert_relative_eq!(dist_inside, 0.2, epsilon = 1e-10);

        // Point close to inner boundary
        let point_near_inner = Point2D { x: 0.45, y: 0.0 };
        let dist_near_inner = domain.distance_to_boundary(&point_near_inner);
        // Inside inner circle, distance to inner boundary: 0.5 - 0.45 = 0.05
        assert_relative_eq!(dist_near_inner, 0.05, epsilon = 1e-10);

        // Point close to outer boundary
        let point_near_outer = Point2D { x: 0.95, y: 0.0 };
        let dist_near_outer = domain.distance_to_boundary(&point_near_outer);
        // Distance to outer boundary: 1.0 - 0.95 = 0.05
        assert_relative_eq!(dist_near_outer, 0.05, epsilon = 1e-10);

        // Point outside outer boundary
        let point_outside = Point2D { x: 1.2, y: 0.0 };
        let dist_outside = domain.distance_to_boundary(&point_outside);
        // Outside outer circle, distance to outer boundary: 1.2 - 1.0 = 0.2
        assert_relative_eq!(dist_outside, 0.2, epsilon = 1e-10);
    }

    #[test]
    fn test_annular_domain_boundary_normal() {
        let domain: AnnularDomain<f64> = AnnularDomain::unit_annulus();

        // Point on inner boundary at (0.5, 0) should have inward normal (towards center)
        let point_inner = Point2D { x: 0.5, y: 0.0 };
        let normal_inner = domain.boundary_normal(&point_inner).unwrap();
        // Inward normal points towards center (0,0), so at (0.5,0) it should be (-1,0)
        assert_relative_eq!(normal_inner.x, -1.0, epsilon = 1e-10);
        assert_relative_eq!(normal_inner.y, 0.0, epsilon = 1e-10);

        // Point on outer boundary at (1,0) should have outward normal (away from center)
        let point_outer = Point2D { x: 1.0, y: 0.0 };
        let normal_outer = domain.boundary_normal(&point_outer).unwrap();
        // Outward normal points away from center, so at (1,0) it should be (1,0)
        assert_relative_eq!(normal_outer.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(normal_outer.y, 0.0, epsilon = 1e-10);

        // Point inside annulus should return None
        let point_inside = Point2D { x: 0.7, y: 0.0 };
        assert!(domain.boundary_normal(&point_inside).is_none());
    }

    #[test]
    fn test_annular_domain_measure() {
        let domain: AnnularDomain<f64> = AnnularDomain::unit_annulus();
        // Area = π(R² - r²) = π(1² - 0.5²) = π(1 - 0.25) = π(0.75)
        let expected_area = std::f64::consts::PI * (1.0 - 0.25);
        assert_relative_eq!(domain.measure(), expected_area, epsilon = 1e-10);

        // Test another annulus
        let domain2: AnnularDomain<f64> = AnnularDomain::new(0.0, 0.0, 1.0, 3.0);
        // Area = π(9 - 1) = π(8)
        let expected_area2 = std::f64::consts::PI * 8.0;
        assert_relative_eq!(domain2.measure(), expected_area2, epsilon = 1e-10);
    }

    #[test]
    fn test_annular_domain_on_boundary() {
        let domain: AnnularDomain<f64> = AnnularDomain::unit_annulus();

        // Points on inner boundary
        assert!(domain.on_boundary(&Point2D { x: 0.5, y: 0.0 }, BoundaryFace::Inner, 1e-10));
        assert!(domain.on_boundary(&Point2D { x: 0.0, y: 0.5 }, BoundaryFace::Inner, 1e-10));

        // Points on outer boundary
        assert!(domain.on_boundary(&Point2D { x: 1.0, y: 0.0 }, BoundaryFace::Outer, 1e-10));
        assert!(domain.on_boundary(&Point2D { x: 0.0, y: 1.0 }, BoundaryFace::Outer, 1e-10));

        // Points not on boundaries
        assert!(!domain.on_boundary(&Point2D { x: 0.7, y: 0.0 }, BoundaryFace::Inner, 1e-10));
        assert!(!domain.on_boundary(&Point2D { x: 0.7, y: 0.0 }, BoundaryFace::Outer, 1e-10));
        assert!(!domain.on_boundary(&Point2D { x: 1.0, y: 0.0 }, BoundaryFace::Bottom, 1e-10));
    }

    #[test]
    fn test_annular_domain_bounds() {
        let domain: AnnularDomain<f64> = AnnularDomain::unit_annulus();
        let (min_point, max_point) = domain.bounds();

        // Bounding box should be from (-1,-1) to (1,1) since outer radius is 1.0
        assert_relative_eq!(min_point.x, -1.0, epsilon = 1e-10);
        assert_relative_eq!(min_point.y, -1.0, epsilon = 1e-10);
        assert_relative_eq!(max_point.x, 1.0, epsilon = 1e-10);
        assert_relative_eq!(max_point.y, 1.0, epsilon = 1e-10);

        let domain2: AnnularDomain<f64> = AnnularDomain::new(2.0, 3.0, 1.0, 2.0);
        let (min_point2, max_point2) = domain2.bounds();

        // Center at (2,3), outer radius 2, so bounds should be (0,1) to (4,5)
        assert_relative_eq!(min_point2.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(min_point2.y, 1.0, epsilon = 1e-10);
        assert_relative_eq!(max_point2.x, 4.0, epsilon = 1e-10);
        assert_relative_eq!(max_point2.y, 5.0, epsilon = 1e-10);
    }
}
