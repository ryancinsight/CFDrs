//! Primitive creation dialog — collects shape parameters from the user.

use crate::application::mesh_ops::primitives::PrimitiveSpec;

/// Available primitive shape types for the dialog.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PrimitiveType {
    Cube,
    Cylinder,
    Sphere,
    Cone,
    Torus,
    Pipe,
    Ellipsoid,
    Capsule,
    Frustum,
    GeodesicSphere,
    RoundedCube,
    Elbow,
}

impl PrimitiveType {
    /// All available primitive types.
    #[must_use]
    pub fn all() -> &'static [Self] {
        &[
            Self::Cube,
            Self::Cylinder,
            Self::Sphere,
            Self::Cone,
            Self::Torus,
            Self::Pipe,
            Self::Ellipsoid,
            Self::Capsule,
            Self::Frustum,
            Self::GeodesicSphere,
            Self::RoundedCube,
            Self::Elbow,
        ]
    }

    /// Display name for this type.
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Cube => "Cube",
            Self::Cylinder => "Cylinder",
            Self::Sphere => "Sphere",
            Self::Cone => "Cone",
            Self::Torus => "Torus",
            Self::Pipe => "Pipe",
            Self::Ellipsoid => "Ellipsoid",
            Self::Capsule => "Capsule",
            Self::Frustum => "Frustum",
            Self::GeodesicSphere => "Geodesic Sphere",
            Self::RoundedCube => "Rounded Cube",
            Self::Elbow => "Elbow",
        }
    }
}

/// Parameters collected from the primitive dialog.
pub struct PrimitiveDialogResult {
    /// The name for the new object.
    pub name: String,
    /// The primitive specification.
    pub spec: PrimitiveSpec,
}

/// Create a default `PrimitiveSpec` for the given type.
#[must_use]
pub fn default_spec(ptype: PrimitiveType) -> PrimitiveSpec {
    match ptype {
        PrimitiveType::Cube => PrimitiveSpec::Cube {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
        },
        PrimitiveType::Cylinder => PrimitiveSpec::Cylinder {
            radius: 0.5,
            height: 1.0,
            segments: 32,
        },
        PrimitiveType::Sphere => PrimitiveSpec::Sphere {
            radius: 0.5,
            segments: 32,
            stacks: 16,
        },
        PrimitiveType::Cone => PrimitiveSpec::Cone {
            radius: 0.5,
            height: 1.0,
            segments: 32,
        },
        PrimitiveType::Torus => PrimitiveSpec::Torus {
            major_radius: 1.0,
            minor_radius: 0.25,
            major_segments: 32,
            minor_segments: 16,
        },
        PrimitiveType::Pipe => PrimitiveSpec::Pipe {
            outer_radius: 0.5,
            inner_radius: 0.4,
            height: 1.0,
            segments: 32,
        },
        PrimitiveType::Ellipsoid => PrimitiveSpec::Ellipsoid {
            semi_x: 0.5,
            semi_y: 0.35,
            semi_z: 0.25,
            segments: 32,
            stacks: 16,
        },
        PrimitiveType::Capsule => PrimitiveSpec::Capsule {
            radius: 0.3,
            cylinder_height: 0.6,
            segments: 32,
            hemisphere_stacks: 8,
        },
        PrimitiveType::Frustum => PrimitiveSpec::Frustum {
            bottom_radius: 0.5,
            top_radius: 0.25,
            height: 1.0,
            segments: 32,
        },
        PrimitiveType::GeodesicSphere => PrimitiveSpec::GeodesicSphere {
            radius: 0.5,
            frequency: 3,
        },
        PrimitiveType::RoundedCube => PrimitiveSpec::RoundedCube {
            width: 1.0,
            height: 1.0,
            depth: 1.0,
            corner_radius: 0.15,
            corner_segments: 4,
        },
        PrimitiveType::Elbow => PrimitiveSpec::Elbow {
            tube_radius: 0.15,
            bend_radius: 0.5,
            bend_angle: std::f64::consts::FRAC_PI_2,
            tube_segments: 16,
            arc_segments: 16,
        },
    }
}
