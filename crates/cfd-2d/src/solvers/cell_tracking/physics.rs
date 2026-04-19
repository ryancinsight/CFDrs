use super::population::OutletZone;

/// Velocity field interface for the cell tracker.
pub trait VelocityFieldInterpolator {
    /// Interpolate velocity at (x, y).  Returns (u, v) in m/s.
    fn velocity_at(&self, x: f64, y: f64) -> (f64, f64);

    /// Check if (x, y) is inside the fluid domain.
    fn is_fluid(&self, x: f64, y: f64) -> bool;

    /// Domain bounding box: (x_min, x_max, y_min, y_max) in meters.
    fn bounds(&self) -> (f64, f64, f64, f64);
}

/// Configuration for the Lagrangian tracker.
#[derive(Debug, Clone)]
pub struct CellTrackerConfig {
    /// Fluid viscosity [Pa.s].  Default: 3.5e-3 (blood).
    pub viscosity: f64,
    /// Fluid density [kg/m3].  Default: 1025.0 (plasma).
    pub fluid_density: f64,
    /// Hydraulic diameter of the parent channel [m], used for lift scaling.
    pub hydraulic_diameter_m: f64,
    /// Max streamwise velocity [m/s] for lift scaling.  If zero, estimated
    /// from the velocity field at the inlet centerline.
    pub u_max: f64,
    /// Outlet zones for classifying exits.  If empty, a default center/peripheral
    /// split at y_mid +/- 25% is used.
    pub outlet_zones: Vec<OutletZone>,
    /// X-coordinate of the bifurcation split plane [m].  Set to 0 to disable
    /// the junction routing correction.
    pub split_x: f64,
    /// Y-coordinate of the dividing streamline at the split plane [m].
    pub dividing_streamline_y: f64,
    /// Pries Phase Separation Model parameters for the bifurcation.
    /// If None, the junction routing correction is disabled.
    pub psm_params: Option<PsmBifurcationParams>,
}

impl Default for CellTrackerConfig {
    fn default() -> Self {
        Self {
            viscosity: 3.5e-3,
            fluid_density: 1025.0,
            hydraulic_diameter_m: 1.0e-3,
            u_max: 0.0,
            outlet_zones: Vec::new(),
            split_x: 0.0,
            dividing_streamline_y: 0.0,
            psm_params: None,
        }
    }
}

/// Parameters for the Pries PSM at a specific bifurcation.
#[derive(Debug, Clone)]
pub struct PsmBifurcationParams {
    /// Fractional blood flow into the wide (center) daughter.
    pub flow_fraction_wide: f64,
    /// Hydraulic diameter of the wide daughter [m].
    pub wide_daughter_dh: f64,
    /// Hydraulic diameter of the narrow daughter [m].
    pub narrow_daughter_dh: f64,
    /// Feed hematocrit.
    pub feed_hematocrit: f64,
}

/// Poiseuille flow in a 2D channel: u(y) = U_max * (1 - (2y/H - 1)^2).
/// Used for unit testing the cell tracker against known analytical solutions.
pub struct PoiseuilleFlow2D {
    /// Maximum center-line velocity [m/s].
    pub u_max: f64,
    /// Channel extent along X axis [m].
    pub width: f64,
    /// Channel extent along Y axis (height) [m].
    pub height: f64,
}

impl VelocityFieldInterpolator for PoiseuilleFlow2D {
    fn velocity_at(&self, _x: f64, y: f64) -> (f64, f64) {
        let s = (2.0 * y / self.height - 1.0).clamp(-1.0, 1.0);
        (self.u_max * (1.0 - s * s), 0.0)
    }
    fn is_fluid(&self, x: f64, y: f64) -> bool {
        x >= 0.0 && x <= self.width && y >= 0.0 && y <= self.height
    }
    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, self.width, 0.0, self.height)
    }
}

/// Asymmetric bifurcation flow: parent channel splits into a wide and narrow
/// daughter at x = x_split.  The wide daughter gets more flow (proportional
/// to width^3 per Hagen-Poiseuille).
pub struct AsymmetricBifurcationFlow {
    /// The width of the parent channel [m].
    pub parent_width_m: f64,
    /// The height of the parent channel [m].
    pub parent_height_m: f64,
    /// The width of the wide daughter channel [m].
    pub wide_daughter_width_m: f64,
    /// The width of the narrow daughter channel [m].
    pub narrow_daughter_width_m: f64,
    /// Length of the domain [m].
    pub length_m: f64,
    /// Input velocity profile scale [m/s].
    pub u_inlet: f64,
    /// The x-coordinate of the split point [m].
    pub x_split: f64,
}

impl AsymmetricBifurcationFlow {
    /// Flow fraction to the wide daughter (Q_wide / Q_total).
    /// HP scaling: Q proportional to w^3 for equal height/length.
    pub fn flow_fraction_wide(&self) -> f64 {
        let q_wide = self.wide_daughter_width_m.powi(3);
        let q_narrow = self.narrow_daughter_width_m.powi(3);
        q_wide / (q_wide + q_narrow)
    }

    /// The dividing streamline y-position at the split plane.
    /// In the parent, the flow above this line goes to the wide daughter
    /// and below goes to the narrow daughter. For a Poiseuille profile,
    /// the exact cumulative flow fraction above x = y/H is
    /// `F_above(x) = 1 - 3x^2 + 2x^3`; this method inverts that relation
    /// so the split line matches the requested wide-daughter flow fraction.
    pub fn dividing_streamline_y(&self) -> f64 {
        self.parent_height_m * inverse_poiseuille_flow_fraction_above(self.flow_fraction_wide())
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
    fn velocity_at(&self, x: f64, y: f64) -> (f64, f64) {
        let h = self.parent_height_m;
        let f_wide = self.flow_fraction_wide();
        let y_div = self.dividing_streamline_y();

        if x < self.x_split {
            // Parent: Poiseuille profile + cross-stream steering near split.
            let s = (2.0 * y / h - 1.0).clamp(-1.0, 1.0);
            let u = self.u_inlet * (1.0 - s * s);

            // Near the split plane, add cross-stream velocity that redirects
            // streamlines toward their destination daughter channel.
            // Cells above the dividing streamline get pushed toward the wide
            // daughter (upper region); cells below get pushed to narrow (lower).
            let approach_factor = ((x - self.x_split + h * 0.5) / (h * 0.5)).clamp(0.0, 1.0);
            let v = if approach_factor > 0.0 {
                let target_y = if y > y_div {
                    // Going to wide daughter: steer toward center of wide region
                    (y_div + h) * 0.5
                } else {
                    // Going to narrow daughter: steer toward center of narrow
                    y_div * 0.5
                };
                let steering = (target_y - y) * approach_factor * 2.0 * self.u_inlet / h;
                steering.clamp(-self.u_inlet * 0.3, self.u_inlet * 0.3)
            } else {
                0.0
            };

            (u, v)
        } else {
            // Daughters: Poiseuille in each sub-channel.
            // Wide daughter: y in [y_div, H].
            // Narrow daughter: y in [0, y_div].
            if y >= y_div {
                let dh = h - y_div;
                if dh < 1e-12 {
                    return (0.0, 0.0);
                }
                let s = (2.0 * (y - y_div) / dh - 1.0).clamp(-1.0, 1.0);
                // U_max = (3/2) * Q / A.  Q = f_wide * Q_parent.
                // Q_parent = (2/3) * u_inlet * H (Poiseuille integral).
                // U_max_daughter = (3/2) * f_wide * (2/3) * u_inlet * H / dh
                //                = f_wide * u_inlet * H / dh
                let u_d = f_wide * self.u_inlet * h / dh;
                (u_d * (1.0 - s * s), 0.0)
            } else {
                if y_div < 1e-12 {
                    return (0.0, 0.0);
                }
                let s = (2.0 * y / y_div - 1.0).clamp(-1.0, 1.0);
                let u_d = (1.0 - f_wide) * self.u_inlet * h / y_div;
                (u_d * (1.0 - s * s), 0.0)
            }
        }
    }

    fn is_fluid(&self, x: f64, y: f64) -> bool {
        x >= 0.0 && x <= self.length_m && y >= 0.0 && y <= self.parent_height_m
    }

    fn bounds(&self) -> (f64, f64, f64, f64) {
        (0.0, self.length_m, 0.0, self.parent_height_m)
    }
}

#[cfg(test)]
mod tests {
    use super::{inverse_poiseuille_flow_fraction_above, AsymmetricBifurcationFlow};

    #[test]
    fn dividing_streamline_inverts_poiseuille_cdf() {
        let flow = AsymmetricBifurcationFlow {
            parent_width_m: 2.0e-3,
            parent_height_m: 2.0e-3,
            wide_daughter_width_m: 1.1e-3,
            narrow_daughter_width_m: 0.9e-3,
            length_m: 0.015,
            u_inlet: 0.05,
            x_split: 0.005,
        };

        let y_div = flow.dividing_streamline_y();
        let x = y_div / flow.parent_height_m;
        let flow_above = 1.0 - 3.0 * x * x + 2.0 * x * x * x;

        assert!((flow_above - flow.flow_fraction_wide()).abs() < 1.0e-12);
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
