//! Rectangular domain geometry implementation

use super::{BoundaryCondition, BoundaryFace, Geometry, Point2D};
use cfd_core::conversion::SafeFromF64;
use nalgebra::RealField;

/// Rectangular domain geometry
#[derive(Debug, Clone)]
pub struct RectangularDomain<T: RealField> {
    /// Minimum x-coordinate
    x_min: T,
    /// Maximum x-coordinate
    x_max: T,
    /// Minimum y-coordinate
    y_min: T,
    /// Maximum y-coordinate
    y_max: T,
}

impl<T: RealField> RectangularDomain<T> {
    /// Create a new rectangular domain
    pub fn new(x_min: T, x_max: T, y_min: T, y_max: T) -> Self {
        Self {
            x_min,
            x_max,
            y_min,
            y_max,
        }
    }

    /// Create a unit square centered at origin
    pub fn unit_square() -> Self
    where
        T: From<f64>,
    {
        Self::new(T::zero(), T::one(), T::zero(), T::one())
    }
}

impl<T: RealField + Copy> Geometry<T> for RectangularDomain<T> {
    fn clone_box(&self) -> Box<dyn Geometry<T>> {
        Box::new(self.clone())
    }

    fn contains(&self, point: &Point2D<T>) -> bool {
        point.x >= self.x_min
            && point.x <= self.x_max
            && point.y >= self.y_min
            && point.y <= self.y_max
    }

    fn distance_to_boundary(&self, point: &Point2D<T>) -> T {
        let dx_left = point.x - self.x_min;
        let dx_right = self.x_max - point.x;
        let dy_bottom = point.y - self.y_min;
        let dy_top = self.y_max - point.y;

        let min_x = if dx_left < dx_right {
            dx_left
        } else {
            dx_right
        };
        let min_y = if dy_bottom < dy_top {
            dy_bottom
        } else {
            dy_top
        };
        if min_x < min_y {
            min_x
        } else {
            min_y
        }
    }

    fn boundary_normal(&self, point: &Point2D<T>) -> Option<Point2D<T>> {
        let tol = T::from_f64_or_zero(1e-10);

        if (point.x - self.x_min).abs() < tol {
            Some(Point2D {
                x: -T::one(),
                y: T::zero(),
            }) // Left boundary
        } else if (point.x - self.x_max).abs() < tol {
            Some(Point2D {
                x: T::one(),
                y: T::zero(),
            }) // Right boundary
        } else if (point.y - self.y_min).abs() < tol {
            Some(Point2D {
                x: T::zero(),
                y: -T::one(),
            }) // Bottom boundary
        } else if (point.y - self.y_max).abs() < tol {
            Some(Point2D {
                x: T::zero(),
                y: T::one(),
            }) // Top boundary
        } else {
            None
        }
    }

    fn boundary_condition(&self, face: BoundaryFace, _s: T) -> BoundaryCondition<T>
    where
        T: Copy,
    {
        // Default to Dirichlet with zero value for all boundaries
        // Specific MMS implementations can override this
        match face {
            BoundaryFace::Bottom | BoundaryFace::Top | BoundaryFace::Left | BoundaryFace::Right => {
                BoundaryCondition::Dirichlet(T::zero())
            }
            _ => BoundaryCondition::Dirichlet(T::zero()),
        }
    }

    fn bounds(&self) -> (Point2D<T>, Point2D<T>) {
        (
            Point2D {
                x: self.x_min,
                y: self.y_min,
            },
            Point2D {
                x: self.x_max,
                y: self.y_max,
            },
        )
    }

    fn boundary_parameter(&self, face: BoundaryFace, point: &Point2D<T>) -> Option<T> {
        let tol = T::from_f64_or_zero(1e-10);
        match face {
            BoundaryFace::Bottom => {
                if (point.y - self.y_min).abs() < tol {
                    Some((point.x - self.x_min) / (self.x_max - self.x_min))
                } else {
                    None
                }
            }
            BoundaryFace::Top => {
                if (point.y - self.y_max).abs() < tol {
                    Some((point.x - self.x_min) / (self.x_max - self.x_min))
                } else {
                    None
                }
            }
            BoundaryFace::Left => {
                if (point.x - self.x_min).abs() < tol {
                    Some((point.y - self.y_min) / (self.y_max - self.y_min))
                } else {
                    None
                }
            }
            BoundaryFace::Right => {
                if (point.x - self.x_max).abs() < tol {
                    Some((point.y - self.y_min) / (self.y_max - self.y_min))
                } else {
                    None
                }
            }
            _ => None,
        }
    }

    fn on_boundary(&self, point: &Point2D<T>, face: BoundaryFace, tolerance: T) -> bool {
        match face {
            BoundaryFace::Bottom => (point.y - self.y_min).abs() < tolerance,
            BoundaryFace::Top => (point.y - self.y_max).abs() < tolerance,
            BoundaryFace::Left => (point.x - self.x_min).abs() < tolerance,
            BoundaryFace::Right => (point.x - self.x_max).abs() < tolerance,
            _ => false,
        }
    }

    fn measure(&self) -> T {
        (self.x_max - self.x_min) * (self.y_max - self.y_min)
    }
}
