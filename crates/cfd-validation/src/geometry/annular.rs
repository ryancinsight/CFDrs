//! Annular domain geometry implementation (ring-shaped region)

use super::{BoundaryCondition, BoundaryFace, Geometry, Point2D};
use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;
use num_traits::Float;

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

impl<T: RealField + Float> AnnularDomain<T> {
    /// Create a new annular domain
    pub fn new(center_x: T, center_y: T, inner_radius: T, outer_radius: T) -> Self {
        assert!(inner_radius < outer_radius, "Inner radius must be less than outer radius");
        Self {
            center_x,
            center_y,
            inner_radius,
            outer_radius,
        }
    }

    /// Create a unit annular domain (inner radius 0.5, outer radius 1.0)
    pub fn unit_annulus() -> Self
    where
        T: SafeFromF64,
    {
        Self::new(
            T::zero(),
            T::zero(),
            SafeFromF64::from_f64_or_zero(0.5),
            T::one()
        )
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

impl<T: RealField + Float + SafeFromF64> Geometry<T> for AnnularDomain<T> {
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
            RealField::min(to_inner, to_outer)
        }
    }

    fn boundary_normal(&self, point: &Point2D<T>) -> Option<Point2D<T>> {
        let distance = self.distance_from_center(point);
        let tol = SafeFromF64::from_f64_or_zero(1e-10);

        if Float::abs(distance - self.inner_radius) < tol {
            // Point is on inner boundary, inward normal (towards center)
            let dx = self.center_x.clone() - point.x.clone();
            let dy = self.center_y.clone() - point.y.clone();
            let norm = Float::sqrt(dx.clone() * dx + dy.clone() * dy);

            if norm > T::zero() {
                Some(Point2D {
                    x: dx / norm.clone(),
                    y: dy / norm,
                })
            } else {
                // Point is at center, return arbitrary normal
                Some(Point2D { x: -T::one(), y: T::zero() })
            }
        } else if Float::abs(distance - self.outer_radius) < tol {
            // Point is on outer boundary, outward normal (away from center)
            let dx = point.x.clone() - self.center_x.clone();
            let dy = point.y.clone() - self.center_y.clone();
            let norm = Float::sqrt(dx.clone() * dx + dy.clone() * dy);

            if norm > T::zero() {
                Some(Point2D {
                    x: dx / norm.clone(),
                    y: dy / norm,
                })
            } else {
                // Point is at center, return arbitrary normal
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
        // For annular domains, we have inner and outer boundaries
        match face {
            BoundaryFace::Inner => {
                // Default to Dirichlet with zero value for inner boundary
                BoundaryCondition::Dirichlet(T::zero())
            }
            BoundaryFace::Outer => {
                // Default to Dirichlet with zero value for outer boundary
                BoundaryCondition::Dirichlet(T::zero())
            }
            _ => BoundaryCondition::Dirichlet(T::zero()),
        }
    }

    fn bounds(&self) -> (Point2D<T>, Point2D<T>) {
        (
            Point2D {
                x: self.center_x.clone() - self.outer_radius.clone(),
                y: self.center_y.clone() - self.outer_radius.clone(),
            },
            Point2D {
                x: self.center_x.clone() + self.outer_radius.clone(),
                y: self.center_y.clone() + self.outer_radius.clone(),
            },
        )
    }

    fn boundary_parameter(&self, face: BoundaryFace, point: &Point2D<T>) -> Option<T> {
        match face {
            BoundaryFace::Inner => {
                let distance = self.distance_from_center(point);
                let tol = SafeFromF64::from_f64_or_zero(1e-10);

                if Float::abs(distance - self.inner_radius) < tol {
                    // Point is on inner boundary, return angular parameter [0, 2π)
                    let angle = self.angle_from_center(point);
                    // Normalize to [0, 2π)
                    let two_pi = T::from_f64(2.0 * std::f64::consts::PI).unwrap_or(T::from_f64(6.28).unwrap_or(T::zero()));
                    Some(if angle >= T::zero() { angle } else { angle + two_pi })
                } else {
                    None
                }
            }
            BoundaryFace::Outer => {
                let distance = self.distance_from_center(point);
                let tol = SafeFromF64::from_f64_or_zero(1e-10);

                if Float::abs(distance - self.outer_radius) < tol {
                    // Point is on outer boundary, return angular parameter [0, 2π)
                    let angle = self.angle_from_center(point);
                    // Normalize to [0, 2π)
                    let two_pi = T::from_f64(2.0 * std::f64::consts::PI).unwrap_or(T::from_f64(6.28).unwrap_or(T::zero()));
                    Some(if angle >= T::zero() { angle } else { angle + two_pi })
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
            BoundaryFace::Inner => {
                Float::abs(distance - self.inner_radius) < tolerance
            }
            BoundaryFace::Outer => {
                Float::abs(distance - self.outer_radius) < tolerance
            }
            _ => false,
        }
    }

    fn measure(&self) -> T {
        // Area of annulus: π(R² - r²) where R is outer radius, r is inner radius
        let pi = T::from_f64(std::f64::consts::PI).unwrap_or(T::from_f64(3.14159).unwrap_or(T::zero()));
        let outer_area = pi.clone() * self.outer_radius.clone() * self.outer_radius.clone();
        let inner_area = pi * self.inner_radius.clone() * self.inner_radius;
        outer_area - inner_area
    }
}
