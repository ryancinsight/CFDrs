//! Annular domain geometry implementation (ring-shaped region)

use super::{BoundaryCondition, BoundaryFace, Geometry2D, Point2D};
use crate::scalar;
use eunomia::FloatElement;
use eunomia::RealField;

/// Annular domain geometry (ring between two concentric circles)
#[derive(Debug, Clone)]
pub struct AnnularDomain<T: RealField> {
    /// Center x-coordinate
    center_x: T,
    /// Center y-coordinate
    center_y: T,
    /// Inner radius
    inner_radius: T,
    /// Outer radius
    outer_radius: T,
}

impl<T: RealField + Copy + FloatElement> AnnularDomain<T> {
    /// Create a new annular domain
    pub fn new(center_x: T, center_y: T, inner_radius: T, outer_radius: T) -> Self {
        assert!(
            inner_radius < outer_radius,
            "Inner radius must be less than outer radius"
        );
        Self {
            center_x,
            center_y,
            inner_radius,
            outer_radius,
        }
    }

    /// Create a unit annular domain (inner radius 0.5, outer radius 1.0)
    pub fn unit_annulus() -> Self {
        Self::new(
            scalar::zero(),
            scalar::zero(),
            scalar::from_f64::<T>(0.5),
            scalar::one(),
        )
    }

    /// Get the distance from center to a point
    fn distance_from_center(&self, point: &Point2D<T>) -> T {
        let dx = point.x - self.center_x;
        let dy = point.y - self.center_y;
        scalar::sqrt(dx * dx + dy * dy)
    }

    /// Get the angle (in radians) from center to a point
    fn angle_from_center(&self, point: &Point2D<T>) -> T {
        let dx = point.x - self.center_x;
        let dy = point.y - self.center_y;
        scalar::atan2(dy, dx)
    }
}

impl<T: RealField + Copy + FloatElement> Geometry2D<T> for AnnularDomain<T> {
    fn clone_box(&self) -> Box<dyn Geometry2D<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, point: &Point2D<T>) -> bool {
        let distance = self.distance_from_center(point);
        distance >= self.inner_radius && distance <= self.outer_radius
    }

    fn distance_to_boundary(&self, point: &Point2D<T>) -> T {
        let distance_from_center = self.distance_from_center(point);

        if distance_from_center <= self.inner_radius {
            // Inside inner circle, distance to inner boundary
            self.inner_radius - distance_from_center
        } else if distance_from_center >= self.outer_radius {
            // Outside outer circle, distance to outer boundary
            distance_from_center - self.outer_radius
        } else {
            // Between inner and outer, distance to nearest boundary
            let to_inner = distance_from_center - self.inner_radius;
            let to_outer = self.outer_radius - distance_from_center;
            if to_inner < to_outer {
                to_inner
            } else {
                to_outer
            }
        }
    }

    fn boundary_normal(&self, point: &Point2D<T>) -> Option<Point2D<T>> {
        let distance = self.distance_from_center(point);
        let tol = scalar::from_f64::<T>(1e-10);

        if scalar::abs(distance - self.inner_radius) < tol {
            // Point is on inner boundary, inward normal (towards center)
            let dx = self.center_x - point.x;
            let dy = self.center_y - point.y;
            let norm = scalar::sqrt(dx * dx + dy * dy);

            if norm > scalar::zero() {
                Some(Point2D::new(dx / norm, dy / norm))
            } else {
                // Point is at center, return arbitrary normal
                Some(Point2D::new(-scalar::one::<T>(), scalar::zero()))
            }
        } else if scalar::abs(distance - self.outer_radius) < tol {
            // Point is on outer boundary, outward normal (away from center)
            let dx = point.x - self.center_x;
            let dy = point.y - self.center_y;
            let norm = scalar::sqrt(dx * dx + dy * dy);

            if norm > scalar::zero() {
                Some(Point2D::new(dx / norm, dy / norm))
            } else {
                // Point is at center, return arbitrary normal
                Some(Point2D::new(scalar::one(), scalar::zero()))
            }
        } else {
            None
        }
    }

    fn boundary_condition(&self, face: BoundaryFace, _s: T) -> BoundaryCondition<T>
    where
        T: Copy,
    {
        // For annular domains, we have inner and outer boundaries
        match face {
            BoundaryFace::Inner => {
                // Default to Dirichlet with zero value for inner boundary
                BoundaryCondition::Dirichlet(scalar::zero())
            }
            BoundaryFace::Outer => {
                // Default to Dirichlet with zero value for outer boundary
                BoundaryCondition::Dirichlet(scalar::zero())
            }
            _ => BoundaryCondition::Dirichlet(scalar::zero()),
        }
    }

    fn bounds(&self) -> (Point2D<T>, Point2D<T>) {
        (
            Point2D::new(
                self.center_x - self.outer_radius,
                self.center_y - self.outer_radius,
            ),
            Point2D::new(
                self.center_x + self.outer_radius,
                self.center_y + self.outer_radius,
            ),
        )
    }

    fn boundary_parameter(&self, face: BoundaryFace, point: &Point2D<T>) -> Option<T> {
        match face {
            BoundaryFace::Inner => {
                let distance = self.distance_from_center(point);
                let tol = scalar::from_f64::<T>(1e-10);

                if scalar::abs(distance - self.inner_radius) < tol {
                    // Point is on inner boundary, return angular parameter [0, 2π)
                    let angle = self.angle_from_center(point);
                    // Normalize to [0, 2π)
                    let two_pi = scalar::from_f64::<T>(2.0 * std::f64::consts::PI);
                    Some(if angle >= scalar::zero() {
                        angle
                    } else {
                        angle + two_pi
                    })
                } else {
                    None
                }
            }
            BoundaryFace::Outer => {
                let distance = self.distance_from_center(point);
                let tol = scalar::from_f64::<T>(1e-10);

                if scalar::abs(distance - self.outer_radius) < tol {
                    // Point is on outer boundary, return angular parameter [0, 2π)
                    let angle = self.angle_from_center(point);
                    // Normalize to [0, 2π)
                    let two_pi = scalar::from_f64::<T>(2.0 * std::f64::consts::PI);
                    Some(if angle >= scalar::zero() {
                        angle
                    } else {
                        angle + two_pi
                    })
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    fn on_boundary(&self, point: &Point2D<T>, face: BoundaryFace, tolerance: T) -> bool {
        let distance = self.distance_from_center(point);

        match face {
            BoundaryFace::Inner => scalar::abs(distance - self.inner_radius) < tolerance,
            BoundaryFace::Outer => scalar::abs(distance - self.outer_radius) < tolerance,
            _ => false,
        }
    }

    fn measure(&self) -> T {
        // Area of annulus: π(R² - r²) where R is outer radius, r is inner radius
        let pi = scalar::from_f64::<T>(std::f64::consts::PI);
        let outer_area = pi * self.outer_radius * self.outer_radius;
        let inner_area = pi * self.inner_radius * self.inner_radius;
        outer_area - inner_area
    }
}
