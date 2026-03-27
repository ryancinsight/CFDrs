//! Named engineering views — standard orthographic and isometric viewpoints.

use std::f64::consts::{FRAC_PI_2, FRAC_PI_4, PI};

/// Standard engineering views matching SolidWorks/Fusion 360 conventions.
///
/// Orthographic views align one axis with the screen normal; isometric views
/// combine two rotations for a 3D perspective.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
pub enum NamedView {
    Front,
    Back,
    Top,
    Bottom,
    Left,
    Right,
    IsoFrontTopRight,
    IsoFrontTopLeft,
    IsoBackTopRight,
    IsoBackTopLeft,
}

/// Camera angles for a named view (Y-up, right-handed).
pub struct ViewAngles {
    /// Horizontal rotation in radians (0 = looking along +Z).
    pub azimuth: f64,
    /// Vertical rotation in radians from the XZ plane toward +Y.
    pub elevation: f64,
}

impl NamedView {
    /// Azimuth and elevation angles for this view.
    #[must_use]
    pub fn angles(self) -> ViewAngles {
        match self {
            Self::Front => ViewAngles { azimuth: 0.0, elevation: 0.0 },
            Self::Back => ViewAngles { azimuth: PI, elevation: 0.0 },
            Self::Top => ViewAngles { azimuth: 0.0, elevation: FRAC_PI_2 - 0.001 },
            Self::Bottom => ViewAngles { azimuth: 0.0, elevation: -(FRAC_PI_2 - 0.001) },
            Self::Left => ViewAngles { azimuth: -FRAC_PI_2, elevation: 0.0 },
            Self::Right => ViewAngles { azimuth: FRAC_PI_2, elevation: 0.0 },
            Self::IsoFrontTopRight => ViewAngles {
                azimuth: FRAC_PI_4,
                elevation: ISO_ELEVATION,
            },
            Self::IsoFrontTopLeft => ViewAngles {
                azimuth: -FRAC_PI_4,
                elevation: ISO_ELEVATION,
            },
            Self::IsoBackTopRight => ViewAngles {
                azimuth: PI - FRAC_PI_4,
                elevation: ISO_ELEVATION,
            },
            Self::IsoBackTopLeft => ViewAngles {
                azimuth: -(PI - FRAC_PI_4),
                elevation: ISO_ELEVATION,
            },
        }
    }

    /// Whether this view is an orthographic (non-perspective) view.
    #[must_use]
    pub fn is_orthographic(self) -> bool {
        matches!(
            self,
            Self::Front | Self::Back | Self::Top | Self::Bottom | Self::Left | Self::Right
        )
    }

    /// Display label for the view (shown on the view cube).
    #[must_use]
    pub fn label(self) -> &'static str {
        match self {
            Self::Front => "FRONT",
            Self::Back => "BACK",
            Self::Top => "TOP",
            Self::Bottom => "BOTTOM",
            Self::Left => "LEFT",
            Self::Right => "RIGHT",
            Self::IsoFrontTopRight => "ISO",
            Self::IsoFrontTopLeft => "ISO",
            Self::IsoBackTopRight => "ISO",
            Self::IsoBackTopLeft => "ISO",
        }
    }

    /// All six orthographic views.
    #[must_use]
    pub fn orthographic_views() -> [Self; 6] {
        [Self::Front, Self::Back, Self::Top, Self::Bottom, Self::Left, Self::Right]
    }
}

/// Standard isometric elevation: arctan(1/√2) ≈ 35.264°.
const ISO_ELEVATION: f64 = 0.6154797086703874;

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn front_view_looks_along_positive_z() {
        let angles = NamedView::Front.angles();
        assert_relative_eq!(angles.azimuth, 0.0, epsilon = 1e-12);
        assert_relative_eq!(angles.elevation, 0.0, epsilon = 1e-12);
    }

    #[test]
    fn orthographic_views_are_orthographic() {
        for view in NamedView::orthographic_views() {
            assert!(view.is_orthographic());
        }
        assert!(!NamedView::IsoFrontTopRight.is_orthographic());
    }

    #[test]
    fn iso_elevation_is_correct() {
        let expected = (1.0_f64 / 2.0_f64.sqrt()).atan();
        assert_relative_eq!(ISO_ELEVATION, expected, epsilon = 1e-10);
    }
}
