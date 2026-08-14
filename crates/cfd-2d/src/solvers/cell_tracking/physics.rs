use aequitas::systems::si::quantities::{
    Dimensionless, DynamicViscosity, Length, MassDensity, Velocity,
};

use super::population::OutletZone;

/// Velocity field interface for the cell tracker.
pub trait VelocityFieldInterpolator {
    /// Interpolate velocity at a physical position.
    fn velocity_at(&self, x: Length, y: Length) -> (Velocity, Velocity);

    /// Check if a physical position is inside the fluid domain.
    fn is_fluid(&self, x: Length, y: Length) -> bool;

    /// Return the physical domain bounding box.
    fn bounds(&self) -> (Length, Length, Length, Length);
}

/// Configuration for the Lagrangian tracker.
#[derive(Debug, Clone)]
pub struct CellTrackerConfig {
    /// Fluid dynamic viscosity.
    pub viscosity: DynamicViscosity,
    /// Fluid mass density.
    pub fluid_density: MassDensity,
    /// Hydraulic diameter of the parent channel used for lift scaling.
    pub hydraulic_diameter: Length,
    /// Maximum streamwise velocity used for lift scaling. If zero, it is
    /// estimated from the velocity field at the inlet centerline.
    pub u_max: Velocity,
    /// Outlet zones for classifying exits. If empty, a default center/
    /// peripheral split at the midpoint is used.
    pub outlet_zones: Vec<OutletZone>,
    /// X-coordinate of the bifurcation split plane. Zero disables the
    /// junction routing correction.
    pub split_x: Length,
    /// Y-coordinate of the dividing streamline at the split plane.
    pub dividing_streamline_y: Length,
    /// Pries Phase Separation Model parameters for the bifurcation.
    pub psm_params: Option<PsmBifurcationParams>,
}

impl Default for CellTrackerConfig {
    fn default() -> Self {
        Self {
            viscosity: DynamicViscosity::from_base(3.5e-3),
            fluid_density: MassDensity::from_base(1025.0),
            hydraulic_diameter: Length::from_base(1.0e-3),
            u_max: Velocity::from_base(0.0),
            outlet_zones: Vec::new(),
            split_x: Length::from_base(0.0),
            dividing_streamline_y: Length::from_base(0.0),
            psm_params: None,
        }
    }
}

/// Parameters for the Pries PSM at a specific bifurcation.
#[derive(Debug, Clone)]
pub struct PsmBifurcationParams {
    /// Fractional blood flow into the wide (center) daughter.
    pub flow_fraction_wide: Dimensionless,
    /// Hydraulic diameter of the wide daughter.
    pub wide_daughter_diameter: Length,
    /// Hydraulic diameter of the narrow daughter.
    pub narrow_daughter_diameter: Length,
    /// Feed hematocrit.
    pub feed_hematocrit: Dimensionless,
}

/// Poiseuille flow in a 2D channel used by analytical cell-tracking tests.
pub struct PoiseuilleFlow2D {
    /// Maximum center-line velocity.
    pub u_max: Velocity,
    /// Channel extent along the X axis.
    pub width: Length,
    /// Channel extent along the Y axis.
    pub height: Length,
}

impl VelocityFieldInterpolator for PoiseuilleFlow2D {
    fn velocity_at(&self, _x: Length, y: Length) -> (Velocity, Velocity) {
        let y = y.into_base();
        let height = self.height.into_base();
        let s = (2.0 * y / height - 1.0).clamp(-1.0, 1.0);
        (
            Velocity::from_base(self.u_max.into_base() * (1.0 - s * s)),
            Velocity::from_base(0.0),
        )
    }

    fn is_fluid(&self, x: Length, y: Length) -> bool {
        let x = x.into_base();
        let y = y.into_base();
        x >= 0.0 && x <= self.width.into_base() && y >= 0.0 && y <= self.height.into_base()
    }

    fn bounds(&self) -> (Length, Length, Length, Length) {
        (
            Length::from_base(0.0),
            self.width,
            Length::from_base(0.0),
            self.height,
        )
    }
}

/// Asymmetric bifurcation flow with a wide and a narrow daughter.
pub struct AsymmetricBifurcationFlow {
    /// Width of the parent channel.
    pub parent_width: Length,
    /// Height of the parent channel.
    pub parent_height: Length,
    /// Width of the wide daughter channel.
    pub wide_daughter_width: Length,
    /// Width of the narrow daughter channel.
    pub narrow_daughter_width: Length,
    /// Length of the flow domain.
    pub length: Length,
    /// Input velocity profile scale.
    pub u_inlet: Velocity,
    /// X-coordinate of the split point.
    pub x_split: Length,
}

impl AsymmetricBifurcationFlow {
    /// Return the flow fraction entering the wide daughter.
    #[must_use]
    pub fn flow_fraction_wide(&self) -> Dimensionless {
        let q_wide = self.wide_daughter_width.into_base().powi(3);
        let q_narrow = self.narrow_daughter_width.into_base().powi(3);
        Dimensionless::from_base(q_wide / (q_wide + q_narrow))
    }

    /// Return the dividing streamline position at the split plane.
    #[must_use]
    pub fn dividing_streamline_y(&self) -> Length {
        Length::from_base(
            self.parent_height.into_base()
                * inverse_poiseuille_flow_fraction_above(self.flow_fraction_wide().into_base()),
        )
    }
}

#[inline]
fn inverse_poiseuille_flow_fraction_above(flow_fraction_above: f64) -> f64 {
    let clamped_fraction = flow_fraction_above.clamp(0.0, 1.0);

    if clamped_fraction <= 0.0 {
        return 1.0;
    }
    if clamped_fraction >= 1.0 {
        return 0.0;
    }

    let inverse_cdf = 0.5 - ((2.0 * clamped_fraction - 1.0).asin() / 3.0).sin();
    inverse_cdf.clamp(0.0, 1.0)
}

impl VelocityFieldInterpolator for AsymmetricBifurcationFlow {
    fn velocity_at(&self, x: Length, y: Length) -> (Velocity, Velocity) {
        let x = x.into_base();
        let y = y.into_base();
        let h = self.parent_height.into_base();
        let f_wide = self.flow_fraction_wide().into_base();
        let y_div = self.dividing_streamline_y().into_base();
        let u_inlet = self.u_inlet.into_base();
        let x_split = self.x_split.into_base();

        let (u, v) = if x < x_split {
            let s = (2.0 * y / h - 1.0).clamp(-1.0, 1.0);
            let u = u_inlet * (1.0 - s * s);

            let approach_factor = ((x - x_split + h * 0.5) / (h * 0.5)).clamp(0.0, 1.0);
            let v = if approach_factor > 0.0 {
                let target_y = if y > y_div {
                    (y_div + h) * 0.5
                } else {
                    y_div * 0.5
                };
                let steering = (target_y - y) * approach_factor * 2.0 * u_inlet / h;
                steering.clamp(-u_inlet * 0.3, u_inlet * 0.3)
            } else {
                0.0
            };

            (u, v)
        } else if y >= y_div {
            let dh = h - y_div;
            if dh < 1e-12 {
                return (Velocity::from_base(0.0), Velocity::from_base(0.0));
            }
            let s = (2.0 * (y - y_div) / dh - 1.0).clamp(-1.0, 1.0);
            let u_d = f_wide * u_inlet * h / dh;
            (u_d * (1.0 - s * s), 0.0)
        } else {
            if y_div < 1e-12 {
                return (Velocity::from_base(0.0), Velocity::from_base(0.0));
            }
            let s = (2.0 * y / y_div - 1.0).clamp(-1.0, 1.0);
            let u_d = (1.0 - f_wide) * u_inlet * h / y_div;
            (u_d * (1.0 - s * s), 0.0)
        };

        (Velocity::from_base(u), Velocity::from_base(v))
    }

    fn is_fluid(&self, x: Length, y: Length) -> bool {
        let x = x.into_base();
        let y = y.into_base();
        x >= 0.0 && x <= self.length.into_base() && y >= 0.0 && y <= self.parent_height.into_base()
    }

    fn bounds(&self) -> (Length, Length, Length, Length) {
        (
            Length::from_base(0.0),
            self.length,
            Length::from_base(0.0),
            self.parent_height,
        )
    }
}

#[cfg(test)]
mod tests {
    use aequitas::systems::si::quantities::{Length, Velocity};

    use super::{inverse_poiseuille_flow_fraction_above, AsymmetricBifurcationFlow};

    #[test]
    fn dividing_streamline_inverts_poiseuille_cdf() {
        let flow = AsymmetricBifurcationFlow {
            parent_width: Length::from_base(2.0e-3),
            parent_height: Length::from_base(2.0e-3),
            wide_daughter_width: Length::from_base(1.1e-3),
            narrow_daughter_width: Length::from_base(0.9e-3),
            length: Length::from_base(0.015),
            u_inlet: Velocity::from_base(0.05),
            x_split: Length::from_base(0.005),
        };

        let y_div = flow.dividing_streamline_y().into_base();
        let x = y_div / flow.parent_height.into_base();
        let flow_above = 1.0 - 3.0 * x * x + 2.0 * x * x * x;

        assert!((flow_above - flow.flow_fraction_wide().into_base()).abs() < 1.0e-12);
        assert!(x > 0.0 && x < 1.0);
    }

    #[test]
    fn inverse_poiseuille_flow_fraction_above_is_monotone() {
        let low = inverse_poiseuille_flow_fraction_above(0.25);
        let mid = inverse_poiseuille_flow_fraction_above(0.50);
        let high = inverse_poiseuille_flow_fraction_above(0.75);

        assert!(low > mid && mid > high);
        assert!((mid - 0.5).abs() < 1.0e-12);
    }
}
