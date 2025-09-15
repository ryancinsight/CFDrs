//! Boundary geometry definitions

use super::specification::BoundaryConditionSpec;
use nalgebra::{Point3, RealField, Vector3};
use serde::{Deserialize, Serialize};

/// Boundary region specification
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BoundaryRegion<T: RealField + Copy> {
    /// Region identifier
    pub id: String,
    /// Region geometry
    pub geometry: BoundaryGeometry<T>,
    /// Associated boundary condition
    pub condition: Option<BoundaryConditionSpec<T>>,
}

/// Boundary geometry types
///
/// Defines the geometric representation of boundary regions for different dimensionalities.
/// Used to specify where boundary conditions should be applied in the computational domain.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub enum BoundaryGeometry<T: RealField + Copy> {
    /// Point boundary (0D) - single point in space
    Point(Point3<T>),
    /// Line boundary (1D) - line segment defined by start and end points
    Line {
        /// Starting point of the line segment
        start: Point3<T>,
        /// Ending point of the line segment
        end: Point3<T>,
    },
    /// Surface boundary (2D) - polygonal surface defined by vertices
    Surface {
        /// Vertices defining the boundary surface (ordered)
        vertices: Vec<Point3<T>>,
    },
    /// Volume boundary (3D) - volumetric region defined by vertices
    Volume {
        /// Vertices defining the boundary volume
        vertices: Vec<Point3<T>>,
    },
    /// Planar boundary - infinite plane defined by point and normal
    Plane {
        /// Point on the plane
        point: Point3<T>,
        /// Normal vector to the plane
        normal: Vector3<T>,
    },
    /// Cylindrical boundary
    Cylinder {
        /// Center of cylinder base
        center: Point3<T>,
        /// Cylinder axis direction
        axis: Vector3<T>,
        /// Cylinder radius
        radius: T,
        /// Cylinder height
        height: T,
    },
    /// Spherical boundary
    Sphere {
        /// Sphere center
        center: Point3<T>,
        /// Sphere radius
        radius: T,
    },
}

impl<T: RealField + Copy> BoundaryGeometry<T> {
    /// Check if a point is on or near the boundary
    pub fn contains_point(&self, point: &Point3<T>, tolerance: T) -> bool {
        match self {
            Self::Point(p) => (point - p).norm() <= tolerance,

            Self::Line { start, end } => {
                let line_vec = end - start;
                let point_vec = point - start;
                let t = point_vec.dot(&line_vec) / line_vec.norm_squared();

                if t >= T::zero() && t <= T::one() {
                    let closest = start + line_vec.scale(t);
                    (point - closest).norm() <= tolerance
                } else {
                    false
                }
            }

            Self::Plane { point: p, normal } => {
                let distance = (point - p).dot(normal).abs();
                distance <= tolerance
            }

            Self::Sphere { center, radius } => {
                let distance = (point - center).norm();
                (distance - *radius).abs() <= tolerance
            }

            Self::Cylinder {
                center,
                axis,
                radius,
                height,
            } => {
                let to_point = point - center;
                let axis_projection = to_point.dot(axis);

                if axis_projection >= T::zero() && axis_projection <= *height {
                    let radial_vec = to_point - axis.scale(axis_projection);
                    (radial_vec.norm() - *radius).abs() <= tolerance
                } else {
                    false
                }
            }

            _ => false, // Complex geometries need more sophisticated checks
        }
    }

    /// Get the dimension of the boundary geometry
    pub fn dimension(&self) -> usize {
        match self {
            Self::Point(_) => 0,
            Self::Line { .. } => 1,
            Self::Surface { .. } | Self::Plane { .. } => 2,
            Self::Volume { .. } | Self::Cylinder { .. } | Self::Sphere { .. } => 3,
        }
    }

    /// Calculate the measure (length/area/volume) of the geometry
    pub fn measure(&self) -> T {
        match self {
            Self::Point(_) => T::zero(),

            Self::Line { start, end } => (end - start).norm(),

            Self::Sphere { radius, .. } => {
                let four = T::from_f64(4.0).unwrap_or_else(T::one);
                let three = T::from_f64(3.0).unwrap_or_else(T::one);
                let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::one);
                four / three * pi * radius.powi(3)
            }

            Self::Cylinder { radius, height, .. } => {
                let pi = T::from_f64(std::f64::consts::PI).unwrap_or_else(T::one);
                pi * radius.powi(2) * *height
            }

            _ => T::zero(), // Complex geometries need specialized calculations
        }
    }
}

impl<T: RealField + Copy> BoundaryRegion<T> {
    /// Create a new boundary region
    pub fn new(id: String, geometry: BoundaryGeometry<T>) -> Self {
        Self {
            id,
            geometry,
            condition: None,
        }
    }

    /// Attach a boundary condition to this region
    #[must_use]
    pub fn with_condition(mut self, condition: BoundaryConditionSpec<T>) -> Self {
        self.condition = Some(condition);
        self
    }

    /// Check if this region has an assigned boundary condition
    pub fn has_condition(&self) -> bool {
        self.condition.is_some()
    }
}
