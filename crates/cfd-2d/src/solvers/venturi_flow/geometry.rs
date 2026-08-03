use crate::scalar::Cfd2dScalar;
use crate::scalar::{self, from_f64};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Venturi throat geometry configuration
///
/// Defines the shape of a Venturi in 2D with:
/// - Inlet section (constant width)
/// - Converging section (linear taper)
/// - Throat section (constant width)
/// - Diverging section (linear expansion)
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiGeometry<T: Cfd2dScalar + Copy> {
    /// Inlet width \[m]
    pub w_inlet: T,
    /// Throat width \[m]
    pub w_throat: T,
    /// Inlet section length \[m]
    pub l_inlet: T,
    /// Converging section length \[m]
    pub l_converge: T,
    /// Throat section length \[m]
    pub l_throat: T,
    /// Diverging (recovery) section length \[m]
    pub l_diverge: T,
    /// Channel height (constant) \[m]
    pub height: T,
}

impl<T: Cfd2dScalar + Copy + FloatElement> VenturiGeometry<T> {
    /// Create standard ISO 5167 Venturi with area ratio 0.5
    pub fn iso_5167_standard() -> Self {
        Self {
            w_inlet: from_f64::<T>(10e-3),
            w_throat: from_f64::<T>(7.07e-3),
            l_inlet: from_f64::<T>(10e-3),
            l_converge: from_f64::<T>(1e-3),
            l_throat: from_f64::<T>(2e-3),
            l_diverge: from_f64::<T>(3e-3),
            height: from_f64::<T>(1.0e-3),
        }
    }

    /// Create custom Venturi geometry
    pub fn new(
        w_inlet: T,
        w_throat: T,
        l_inlet: T,
        l_converge: T,
        l_throat: T,
        l_diverge: T,
        height: T,
    ) -> Self {
        Self {
            w_inlet,
            w_throat,
            l_inlet,
            l_converge,
            l_throat,
            l_diverge,
            height,
        }
    }

    /// Explicitly opt into the built-in center-clustering recommendation for
    /// `VenturiSolver2D::new_stretched`.
    #[must_use]
    pub fn recommended_center_clustering_beta(&self) -> T {
        let four = from_f64::<T>(4.0);
        let max_beta = from_f64::<T>(0.9);
        let candidate = scalar::one::<T>() - four * self.w_throat / self.w_inlet;
        <T as NumericElement>::min_scalar(
            <T as NumericElement>::max_scalar(candidate, scalar::zero::<T>()),
            max_beta,
        )
    }

    /// Calculate area ratio (A_throat / A_inlet)
    pub fn area_ratio(&self) -> T {
        self.w_throat / self.w_inlet
    }

    /// Calculate total length
    pub fn total_length(&self) -> T {
        self.l_inlet + self.l_converge + self.l_throat + self.l_diverge
    }

    /// Get inlet cross-sectional area \[m²]
    pub fn area_inlet(&self) -> T {
        self.w_inlet * self.height
    }

    /// Get throat cross-sectional area \[m²]
    pub fn area_throat(&self) -> T {
        self.w_throat * self.height
    }

    /// Check if a point (x, y) is within the fluid domain
    pub fn contains(&self, x: T, y: T) -> bool {
        let x_inlet_end = self.l_inlet;
        let x_converge_end = x_inlet_end + self.l_converge;
        let x_throat_end = x_converge_end + self.l_throat;
        let x_diverge_end = x_throat_end + self.l_diverge;

        if x < scalar::zero::<T>() || x > x_diverge_end {
            return false;
        }

        let w_local = if x < x_inlet_end {
            self.w_inlet
        } else if x < x_converge_end {
            let frac = (x - x_inlet_end) / self.l_converge;
            self.w_inlet + frac * (self.w_throat - self.w_inlet)
        } else if x < x_throat_end {
            self.w_throat
        } else if x < x_diverge_end {
            let frac = (x - x_throat_end) / self.l_diverge;
            self.w_throat + frac * (self.w_inlet - self.w_throat)
        } else {
            self.w_inlet
        };

        let half_w = w_local / from_f64::<T>(2.0);
        eunomia::NumericElement::abs(y) <= half_w
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn test_venturi_geometry_iso() {
        let geom = VenturiGeometry::<f64>::iso_5167_standard();

        assert!(geom.w_inlet > geom.w_throat);
        assert!(geom.area_ratio() < 1.0);
        assert_relative_eq!(geom.area_ratio(), 0.707, epsilon = 0.01);
    }
}
