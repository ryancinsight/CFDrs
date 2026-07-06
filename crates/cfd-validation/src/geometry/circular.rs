//! Circular domain geometry implementation

use super::{BoundaryCondition, BoundaryFace, Geometry2D, Point2D};
use crate::scalar;
use eunomia::FloatElement;
use eunomia::RealField;

/// Circular domain geometry (disk)
#[derive(Debug, Clone)]
pub struct CircularDomain<T: RealField> {
    /// Center x-coordinate
    center_x: T,
    /// Center y-coordinate
    center_y: T,
    /// Radius of the circular domain
    radius: T,
}

impl<T: RealField + Copy + FloatElement> CircularDomain<T> {
    /// Create a new circular domain
    pub fn new(center_x: T, center_y: T, radius: T) -> Self {
        Self {
            center_x,
            center_y,
            radius,
        }
    }

    /// Create a unit disk centered at origin
    pub fn unit_disk() -> Self
    where
        T: From<f64>,
    {
        Self::new(scalar::zero(), scalar::zero(), scalar::one())
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

impl<T: RealField + Copy + FloatElement> Geometry2D<T> for CircularDomain<T> {
    fn clone_box(&self) -> Box<dyn Geometry2D<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, point: &Point2D<T>) -> bool {
        self.distance_from_center(point) <= self.radius
    }

    fn distance_to_boundary(&self, point: &Point2D<T>) -> T {
        let distance_from_center = self.distance_from_center(point);
        if distance_from_center <= self.radius {
            // Inside the circle, distance to boundary is radius - distance_from_center
            self.radius - distance_from_center
        } else {
            // Outside the circle, distance to boundary is distance_from_center - radius
            distance_from_center - self.radius
        }
    }

    fn boundary_normal(&self, point: &Point2D<T>) -> Option<Point2D<T>> {
        let distance = self.distance_from_center(point);
        let tol = scalar::from_f64::<T>(1e-10);

        if scalar::abs(distance - self.radius) < tol {
            // Point is on the boundary, compute outward normal
            let dx = point.x - self.center_x;
            let dy = point.y - self.center_y;
            let norm = scalar::sqrt(dx * dx + dy * dy);

            if norm > scalar::zero() {
                Some(Point2D::new(dx / norm, dy / norm))
            } else {
                // Point is at the center, return arbitrary normal (e.g., along x-axis)
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
        // For circular domains, we only have one boundary face: the circumference
        match face {
            BoundaryFace::Outer => {
                // Default to Dirichlet with zero value
                // Specific MMS implementations can override this
                BoundaryCondition::Dirichlet(scalar::zero())
            }
            _ => BoundaryCondition::Dirichlet(scalar::zero()),
        }
    }

    fn bounds(&self) -> (Point2D<T>, Point2D<T>) {
        (
            Point2D::new(self.center_x - self.radius, self.center_y - self.radius),
            Point2D::new(self.center_x + self.radius, self.center_y + self.radius),
        )
    }

    fn boundary_parameter(&self, face: BoundaryFace, point: &Point2D<T>) -> Option<T> {
        match face {
            BoundaryFace::Outer => {
                let distance = self.distance_from_center(point);
                let tol = scalar::from_f64::<T>(1e-10);

                if scalar::abs(distance - self.radius) < tol {
                    // Point is on boundary, return angular parameter [0, 2π)
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
        match face {
            BoundaryFace::Outer => {
                let distance = self.distance_from_center(point);
                scalar::abs(distance - self.radius) < tolerance
            }
            _ => false,
        }
    }

    fn measure(&self) -> T {
        // Area of circle: πr²
        let pi = scalar::from_f64::<T>(std::f64::consts::PI);
        pi * self.radius * self.radius
    }
}
