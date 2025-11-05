//! Circular domain geometry implementation

use super::{BoundaryCondition, BoundaryFace, Geometry, Point2D};
use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;
use num_traits::Float;

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

impl<T: RealField + Float> CircularDomain<T> {
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
        Self::new(T::zero(), T::zero(), T::one())
    }

    /// Get the distance from center to a point
    fn distance_from_center(&self, point: &Point2D<T>) -> T {
        let dx = point.x.clone() - self.center_x.clone();
        let dy = point.y.clone() - self.center_y.clone();
        Float::sqrt(dx.clone() * dx + dy.clone() * dy)
    }

    /// Get the angle (in radians) from center to a point
    fn angle_from_center(&self, point: &Point2D<T>) -> T {
        let dx = point.x.clone() - self.center_x.clone();
        let dy = point.y.clone() - self.center_y.clone();
        RealField::atan2(dy, dx)
    }
}

impl<T: RealField + num_traits::Float> Geometry<T> for CircularDomain<T> {
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
        let tol = T::from_f64(1e-10).unwrap_or(T::from_f64(1e-10).unwrap_or_else(|| T::zero()));

        if Float::abs(distance - self.radius) < tol {
            // Point is on the boundary, compute outward normal
            let dx = point.x.clone() - self.center_x.clone();
            let dy = point.y.clone() - self.center_y.clone();
            let norm = Float::sqrt(dx.clone() * dx.clone() + dy.clone() * dy.clone());

            if norm > T::zero() {
                Some(Point2D {
                    x: dx / norm.clone(),
                    y: dy / norm,
                })
            } else {
                // Point is at the center, return arbitrary normal (e.g., along x-axis)
                Some(Point2D { x: T::one(), y: T::zero() })
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
                BoundaryCondition::Dirichlet(T::zero())
            }
            _ => BoundaryCondition::Dirichlet(T::zero()),
        }
    }

    fn bounds(&self) -> (Point2D<T>, Point2D<T>) {
        (
            Point2D {
                x: self.center_x - self.radius,
                y: self.center_y - self.radius,
            },
            Point2D {
                x: self.center_x + self.radius,
                y: self.center_y + self.radius,
            },
        )
    }

    fn boundary_parameter(&self, face: BoundaryFace, point: &Point2D<T>) -> Option<T> {
        match face {
            BoundaryFace::Outer => {
                let distance = self.distance_from_center(point);
                let tol = T::from_f64(1e-10).unwrap_or(T::from_f64(1e-10).unwrap_or_else(|| T::zero()));

                if Float::abs(distance - self.radius) < tol {
                    // Point is on boundary, return angular parameter [0, 2π)
                    let angle = self.angle_from_center(point);
                    // Normalize to [0, 2π)
                    let two_pi = T::from_f64(2.0 * std::f64::consts::PI)
                        .unwrap_or(T::from_f64(2.0 * std::f64::consts::PI).unwrap_or_else(|| T::from_f64(6.28).unwrap_or(T::zero())));
                    Some(if angle >= T::zero() { angle } else { angle + two_pi })
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
                Float::abs(distance - self.radius) < tolerance
            }
            _ => false,
        }
    }

    fn measure(&self) -> T {
        // Area of circle: πr²
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or(T::from_f64(3.14159).unwrap_or(T::zero()));
        pi * self.radius.clone() * self.radius
    }
}
