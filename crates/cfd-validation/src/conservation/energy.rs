//! Energy conservation checker for CFD simulations
//!
//! Validates that total energy is conserved in the system

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use crate::scalar;
use cfd_core::error::Result;
use eunomia::FloatElement;
use eunomia::RealField;
use leto::Array2;

/// Energy conservation checker
pub struct EnergyConservationChecker<T: RealField + Copy> {
    tolerance: T,
    nx: usize,
    ny: usize,
    density: T,
    specific_heat: T,
}

impl<T: RealField + Copy + FloatElement> EnergyConservationChecker<T> {
    /// Create new energy conservation checker
    pub fn new(tolerance: T, nx: usize, ny: usize, density: T, specific_heat: T) -> Self {
        Self {
            tolerance,
            nx,
            ny,
            density,
            specific_heat,
        }
    }

    /// Check energy conservation for a 2D flow field with temperature
    ///
    /// Validates the energy equation:
    /// ∂(ρcₚT)/∂t + ∇·(ρcₚuT) = ∇·(k∇T) + Φ + Q
    pub fn check_energy_2d(
        &self,
        temperature: &Array2<T>,
        temperature_prev: &Array2<T>,
        u: &Array2<T>,
        v: &Array2<T>,
        thermal_conductivity: T,
        dt: T,
        dx: T,
        dy: T,
        source: Option<&Array2<T>>,
    ) -> Result<ConservationReport<T>> {
        assert_eq!(temperature.shape()[0], self.nx);
        assert_eq!(temperature.shape()[1], self.ny);
        assert_eq!(u.shape()[0], self.nx);
        assert_eq!(v.shape()[0], self.nx);

        let mut max_residual = scalar::zero::<T>();
        let mut total_residual = scalar::zero::<T>();
        let mut total_energy = scalar::zero::<T>();
        let mut total_energy_prev = scalar::zero::<T>();
        let mut count = 0;

        // Check energy conservation at interior points
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let t_center = temperature[[i, j]];
                let t_prev = temperature_prev[[i, j]];

                // Accumulate total energy
                total_energy += self.density * self.specific_heat * t_center * dx * dy;
                total_energy_prev += self.density * self.specific_heat * t_prev * dx * dy;

                // Time derivative: ∂(ρcₚT)/∂t
                let dtdt = self.density * self.specific_heat * (t_center - t_prev) / dt;

                // Convective term: ∇·(ρcₚuT) using central differences
                let conv = self.density
                    * self.specific_heat
                    * (u[[i, j]] * (temperature[[i + 1, j]] - temperature[[i - 1, j]])
                        / (scalar::from_f64::<T>(2.0) * dx)
                        + v[[i, j]] * (temperature[[i, j + 1]] - temperature[[i, j - 1]])
                            / (scalar::from_f64::<T>(2.0) * dy));

                // Diffusive term: ∇·(k∇T) using central differences
                let diff = thermal_conductivity
                    * ((temperature[[i + 1, j]] - scalar::from_f64::<T>(2.0) * t_center
                        + temperature[[i - 1, j]])
                        / (dx * dx)
                        + (temperature[[i, j + 1]] - scalar::from_f64::<T>(2.0) * t_center
                            + temperature[[i, j - 1]])
                            / (dy * dy));

                // Source term
                let q = source.map_or(scalar::zero::<T>(), |s| s[[i, j]]);

                // Energy residual
                let residual = dtdt + conv - diff - q;

                let abs_residual = scalar::abs(residual);
                max_residual = scalar::max(max_residual, abs_residual);
                total_residual += abs_residual;
                count += 1;
            }
        }

        let avg_residual = if count > 0 {
            total_residual / scalar::from_usize::<T>(count)
        } else {
            scalar::zero::<T>()
        };

        // Global energy change
        let energy_change = scalar::abs(total_energy - total_energy_prev);
        let relative_change = if total_energy_prev == scalar::zero::<T>() {
            scalar::zero::<T>()
        } else {
            energy_change / scalar::abs(total_energy_prev)
        };

        let mut report = ConservationReport::new(
            "Energy Conservation (2D)".to_string(),
            max_residual,
            self.tolerance,
        );

        report.add_detail("max_residual".to_string(), max_residual);
        report.add_detail("avg_residual".to_string(), avg_residual);
        report.add_detail("total_energy".to_string(), total_energy);
        report.add_detail("energy_change".to_string(), energy_change);
        report.add_detail("relative_change".to_string(), relative_change);
        report.add_detail(
            "grid_points_checked".to_string(),
            scalar::from_usize::<T>(count),
        );

        Ok(report)
    }

    /// Check total kinetic energy conservation
    pub fn check_kinetic_energy(
        &self,
        u: &Array2<T>,
        v: &Array2<T>,
        u_prev: &Array2<T>,
        v_prev: &Array2<T>,
        dx: T,
        dy: T,
    ) -> Result<ConservationReport<T>> {
        let mut ke_current = scalar::zero::<T>();
        let mut ke_prev = scalar::zero::<T>();

        let half = scalar::from_f64::<T>(0.5);

        for i in 0..self.nx {
            for j in 0..self.ny {
                let u_current = u[[i, j]];
                let v_current = v[[i, j]];
                let u_previous = u_prev[[i, j]];
                let v_previous = v_prev[[i, j]];
                let u2 = u_current * u_current + v_current * v_current;
                let u2_prev = u_previous * u_previous + v_previous * v_previous;

                ke_current += half * self.density * u2 * dx * dy;
                ke_prev += half * self.density * u2_prev * dx * dy;
            }
        }

        let ke_change = scalar::abs(ke_current - ke_prev);
        let relative_change = if ke_prev == scalar::zero::<T>() {
            scalar::zero::<T>()
        } else {
            ke_change / scalar::abs(ke_prev)
        };

        let mut report = ConservationReport::new(
            "Kinetic Energy Conservation".to_string(),
            relative_change,
            self.tolerance,
        );

        report.add_detail("kinetic_energy".to_string(), ke_current);
        report.add_detail("ke_change".to_string(), ke_change);
        report.add_detail("relative_change".to_string(), relative_change);

        Ok(report)
    }
}

impl<T: RealField + Copy + FloatElement> ConservationChecker<T> for EnergyConservationChecker<T> {
    type FlowField = Array2<T>;

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        // For generic check, assume steady state temperature field
        let temperature = field;
        let u = Array2::zeros([self.nx, self.ny]);
        let v = Array2::zeros([self.nx, self.ny]);
        let thermal_conductivity = scalar::from_f64::<T>(0.025);
        let dt = scalar::from_f64::<T>(1e-3);
        let dx = scalar::one::<T>();
        let dy = scalar::one::<T>();

        self.check_energy_2d(
            temperature,
            temperature, // Use same field for steady state
            &u,
            &v,
            thermal_conductivity,
            dt,
            dx,
            dy,
            None,
        )
    }

    fn name(&self) -> &'static str {
        "Energy Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use leto::Array2;

    #[test]
    fn test_uniform_temperature_conservation() {
        let nx = 10;
        let ny = 10;
        let checker = EnergyConservationChecker::<f64>::new(1e-10, nx, ny, 1.0, 1000.0);

        // Uniform temperature field (steady state)
        let temperature = Array2::from_elem([nx, ny], 300.0);
        let u = Array2::zeros([nx, ny]);
        let v = Array2::zeros([nx, ny]);

        let report = checker
            .check_energy_2d(
                &temperature,
                &temperature,
                &u,
                &v,
                0.025,
                1e-3,
                0.1,
                0.1,
                None,
            )
            .expect("Energy conservation check should succeed with valid inputs");

        // For uniform steady temperature, residual should be zero
        assert!(report.error < 1e-10);
        assert!(report.is_conserved);
    }

    #[test]
    fn test_linear_temperature_profile() {
        let nx = 10;
        let ny = 10;
        let checker = EnergyConservationChecker::<f64>::new(1e-6, nx, ny, 1.0, 1000.0);

        // Linear temperature gradient
        let mut temperature = Array2::zeros([nx, ny]);
        for i in 0..nx {
            for j in 0..ny {
                temperature[[i, j]] = 300.0 + 10.0 * (i as f64 / (nx - 1) as f64);
            }
        }

        let u = Array2::zeros([nx, ny]);
        let v = Array2::zeros([nx, ny]);

        let report = checker
            .check_energy_2d(
                &temperature,
                &temperature,
                &u,
                &v,
                0.025,
                1e-3,
                0.1,
                0.1,
                None,
            )
            .expect("Linear temperature profile energy check should succeed");

        // For steady linear profile with no flow, should be conserved
        assert!(report.error < 1e-10);
        assert!(report.is_conserved);
    }

    #[test]
    fn test_kinetic_energy_conservation() {
        let nx = 10;
        let ny = 10;
        let checker = EnergyConservationChecker::<f64>::new(1e-10, nx, ny, 1.0, 1000.0);

        // Uniform velocity field
        let u = Array2::from_elem([nx, ny], 1.0);
        let v = Array2::from_elem([nx, ny], 0.5);

        // Check kinetic energy conservation (steady state)
        let report = checker
            .check_kinetic_energy(&u, &v, &u, &v, 0.1, 0.1)
            .expect("Kinetic energy conservation check should succeed");

        // For steady flow, kinetic energy should be conserved
        assert!(report.error < 1e-10);
        assert!(report.is_conserved);
    }
}
