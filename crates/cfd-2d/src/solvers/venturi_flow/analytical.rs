use super::geometry::VenturiGeometry;
use crate::scalar::Cfd2dScalar;
use crate::scalar::{self, from_f64};
use eunomia::FloatElement;

/// Analytical Venturi solution based on Bernoulli equation
pub struct BernoulliVenturi<T: Cfd2dScalar + Copy> {
    geometry: VenturiGeometry<T>,
    /// Inlet velocity \[m/s]
    pub u_inlet: T,
    /// Inlet pressure \[Pa]
    pub p_inlet: T,
    /// Fluid density [kg/m³]
    pub rho: T,
}

impl<T: Cfd2dScalar + Copy + FloatElement> BernoulliVenturi<T> {
    /// Create new Bernoulli solution
    pub fn new(geometry: VenturiGeometry<T>, u_inlet: T, p_inlet: T, rho: T) -> Self {
        Self {
            geometry,
            u_inlet,
            p_inlet,
            rho,
        }
    }

    /// Calculate velocity at throat (from mass conservation)
    pub fn velocity_throat(&self) -> T {
        self.u_inlet / self.geometry.area_ratio()
    }

    /// Calculate pressure at throat (from Bernoulli equation)
    pub fn pressure_throat(&self) -> T {
        let one_half = from_f64::<T>(0.5);
        let u_throat = self.velocity_throat();
        let dynamic_pressure_inlet = one_half * self.rho * self.u_inlet * self.u_inlet;
        let dynamic_pressure_throat = one_half * self.rho * u_throat * u_throat;

        self.p_inlet + (dynamic_pressure_inlet - dynamic_pressure_throat)
    }

    /// Calculate pressure coefficient at throat
    pub fn pressure_coefficient_throat(&self) -> T {
        let ar = self.geometry.area_ratio();
        let one = scalar::one::<T>();
        one - (one / ar) * (one / ar)
    }

    /// Calculate pressure recovery in diffuser (outlet)
    pub fn pressure_recovery_ideal(&self) -> T {
        scalar::zero::<T>()
    }
}

/// Venturi with viscous friction loss correction
pub struct ViscousVenturi<T: Cfd2dScalar + Copy> {
    bernoulli: BernoulliVenturi<T>,
    /// Loss coefficient (typical: 0.1-0.3)
    pub loss_coefficient: T,
}

impl<T: Cfd2dScalar + Copy + FloatElement> ViscousVenturi<T> {
    /// Create Venturi with friction losses
    pub fn new(
        geometry: VenturiGeometry<T>,
        u_inlet: T,
        p_inlet: T,
        rho: T,
        loss_coefficient: T,
    ) -> Self {
        let bernoulli = BernoulliVenturi::new(geometry, u_inlet, p_inlet, rho);
        Self {
            bernoulli,
            loss_coefficient,
        }
    }

    /// Calculate outlet pressure with friction loss
    pub fn pressure_outlet_with_loss(&self) -> T {
        let one_half = from_f64::<T>(0.5);
        let dynamic_pressure_inlet =
            one_half * self.bernoulli.rho * self.bernoulli.u_inlet * self.bernoulli.u_inlet;

        self.bernoulli.p_inlet - self.loss_coefficient * dynamic_pressure_inlet
    }

    /// Calculate pressure recovery coefficient (real)
    pub fn pressure_recovery_coefficient(&self) -> T {
        -self.loss_coefficient
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::assert_relative_eq;

    #[test]
    fn test_bernoulli_venturi_mass_conservation() {
        let geom = VenturiGeometry::<f64>::iso_5167_standard();
        let bernoulli = BernoulliVenturi::new(geom.clone(), 1.0, 101325.0, 1000.0);

        let q_inlet = geom.area_inlet() * bernoulli.u_inlet;
        let q_throat = geom.area_throat() * bernoulli.velocity_throat();

        assert_relative_eq!(q_inlet, q_throat, epsilon = 1e-10);
    }

    #[test]
    fn test_bernoulli_pressure_drop() {
        let geom = VenturiGeometry::<f64>::iso_5167_standard();
        let bernoulli = BernoulliVenturi::new(geom, 2.0, 101325.0, 1000.0);

        let p_throat = bernoulli.pressure_throat();
        assert!(p_throat < bernoulli.p_inlet);
    }

    #[test]
    fn test_viscous_venturi_recovery_loss() {
        let geom = VenturiGeometry::<f64>::iso_5167_standard();
        let viscous = ViscousVenturi::new(geom, 1.0, 101325.0, 1000.0, 0.15);

        let p_outlet = viscous.pressure_outlet_with_loss();
        assert!(p_outlet < viscous.bernoulli.p_inlet);
    }
}
