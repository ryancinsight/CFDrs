//! Dimension specifications for engineering drawings.

use serde::{Deserialize, Serialize};

/// Tolerance specification for a dimension.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct Tolerance {
    /// Upper deviation (positive value).
    pub upper: f64,
    /// Lower deviation (negative value stored as positive magnitude).
    pub lower: f64,
}

/// A linear dimension between two points.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct LinearDimension {
    /// Start point in 3D model space.
    pub start: [f64; 3],
    /// End point in 3D model space.
    pub end: [f64; 3],
    /// Perpendicular offset of the dimension line from geometry (mm on sheet).
    pub offset_mm: f64,
    /// Optional text override (replaces computed value).
    pub text_override: Option<String>,
    /// Optional tolerance.
    pub tolerance: Option<Tolerance>,
}

impl LinearDimension {
    /// Compute the measured distance between start and end points.
    #[must_use]
    pub fn measured_value(&self) -> f64 {
        let dx = self.end[0] - self.start[0];
        let dy = self.end[1] - self.start[1];
        let dz = self.end[2] - self.start[2];
        (dx * dx + dy * dy + dz * dz).sqrt()
    }
}

/// An angular dimension between two lines.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct AngularDimension {
    /// Vertex of the angle (where the two lines meet) in 3D.
    pub vertex: [f64; 3],
    /// Direction of the first line from the vertex.
    pub direction_a: [f64; 3],
    /// Direction of the second line from the vertex.
    pub direction_b: [f64; 3],
    /// Radius of the arc used to display the angle on the sheet (mm).
    pub arc_radius_mm: f64,
    /// Optional text override.
    pub text_override: Option<String>,
}

impl AngularDimension {
    /// Compute the angle in degrees between the two directions.
    #[must_use]
    pub fn measured_angle_deg(&self) -> f64 {
        let a = nalgebra::Vector3::new(
            self.direction_a[0],
            self.direction_a[1],
            self.direction_a[2],
        );
        let b = nalgebra::Vector3::new(
            self.direction_b[0],
            self.direction_b[1],
            self.direction_b[2],
        );
        let na = a.norm();
        let nb = b.norm();
        if na < 1e-12 || nb < 1e-12 {
            return 0.0;
        }
        let cos_angle = (a.dot(&b) / (na * nb)).clamp(-1.0, 1.0);
        cos_angle.acos().to_degrees()
    }
}

/// A diameter dimension on a circular feature.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct DiameterDimension {
    /// Center of the circle in 3D model space.
    pub center: [f64; 3],
    /// Diameter value in model units.
    pub diameter: f64,
    /// Optional text override.
    pub text_override: Option<String>,
    /// Optional tolerance.
    pub tolerance: Option<Tolerance>,
}

/// A radius dimension on a circular feature.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub struct RadiusDimension {
    /// Center of the arc in 3D model space.
    pub center: [f64; 3],
    /// Radius value in model units.
    pub radius: f64,
    /// Optional text override.
    pub text_override: Option<String>,
}

/// Union of all dimension types.
#[derive(Clone, Debug, Serialize, Deserialize)]
pub enum DimensionSpec {
    /// Linear distance between two points.
    Linear(LinearDimension),
    /// Angle between two lines.
    Angular(AngularDimension),
    /// Diameter of a circle.
    Diameter(DiameterDimension),
    /// Radius of an arc.
    Radius(RadiusDimension),
}
