//! Energy conservation checker for CFD simulations
//!
//! Validates that total energy is conserved in the system

use super::report::ConservationReport;
use super::traits::ConservationChecker;
use cfd_core::error::Result;
use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;

/// Energy conservation checker
pub struct EnergyConservationChecker<T: RealField + Copy> {
    tolerance: T,
    nx: usize,
    ny: usize,
    density: T,
    specific_heat: T,
}

impl<T: RealField + Copy + FromPrimitive> EnergyConservationChecker<T> {
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
        temperature: &DMatrix<T>,
        temperature_prev: &DMatrix<T>,
        u: &DMatrix<T>,
        v: &DMatrix<T>,
        thermal_conductivity: T,
        dt: T,
        dx: T,
        dy: T,
        source: Option<&DMatrix<T>>,
    ) -> Result<ConservationReport<T>> {
        assert_eq!(temperature.nrows(), self.nx);
        assert_eq!(temperature.ncols(), self.ny);
        assert_eq!(u.nrows(), self.nx);
        assert_eq!(v.nrows(), self.nx);

        let mut max_residual = T::zero();
        let mut total_residual = T::zero();
        let mut total_energy = T::zero();
        let mut total_energy_prev = T::zero();
        let mut count = 0;

        // Check energy conservation at interior points
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let t_center = temperature[(i, j)];
                let t_prev = temperature_prev[(i, j)];

                // Accumulate total energy
                total_energy =
                    total_energy + self.density * self.specific_heat * t_center * dx * dy;
                total_energy_prev =
                    total_energy_prev + self.density * self.specific_heat * t_prev * dx * dy;

                // Time derivative: ∂(ρcₚT)/∂t
                let dtdt = self.density * self.specific_heat * (t_center - t_prev) / dt;

                // Convective term: ∇·(ρcₚuT) using central differences
                let conv = self.density
                    * self.specific_heat
                    * (u[(i, j)] * (temperature[(i + 1, j)] - temperature[(i - 1, j)])
                        / (T::from_f64(2.0).unwrap_or(T::one()) * dx)
                        + v[(i, j)] * (temperature[(i, j + 1)] - temperature[(i, j - 1)])
                            / (T::from_f64(2.0).unwrap_or(T::one()) * dy));

                // Diffusive term: ∇·(k∇T) using central differences
                let diff = thermal_conductivity
                    * ((temperature[(i + 1, j)] - T::from_f64(2.0).unwrap_or(T::one()) * t_center
                        + temperature[(i - 1, j)])
                        / (dx * dx)
                        + (temperature[(i, j + 1)]
                            - T::from_f64(2.0).unwrap_or(T::one()) * t_center
                            + temperature[(i, j - 1)])
                            / (dy * dy));

                // Source term
                let q = source.map(|s| s[(i, j)]).unwrap_or(T::zero());

                // Energy residual
                let residual = dtdt + conv - diff - q;

                max_residual = max_residual.max(residual.abs());
                total_residual = total_residual + residual.abs();
                count += 1;
            }
        }

        let avg_residual = if count > 0 {
            total_residual / T::from_usize(count).unwrap_or(T::one())
        } else {
            T::zero()
        };

        // Global energy change
        let energy_change = (total_energy - total_energy_prev).abs();
        let relative_change = if total_energy_prev != T::zero() {
            energy_change / total_energy_prev.abs()
        } else {
            T::zero()
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
            T::from_usize(count).unwrap_or(T::zero()),
        );

        Ok(report)
    }

    /// Check total kinetic energy conservation
    pub fn check_kinetic_energy(
        &self,
        u: &DMatrix<T>,
        v: &DMatrix<T>,
        u_prev: &DMatrix<T>,
        v_prev: &DMatrix<T>,
        dx: T,
        dy: T,
    ) -> Result<ConservationReport<T>> {
        let mut ke_current = T::zero();
        let mut ke_prev = T::zero();

        let half = T::from_f64(0.5).unwrap_or(T::one() / (T::one() + T::one()));

        for i in 0..self.nx {
            for j in 0..self.ny {
                let u2 = u[(i, j)].powi(2) + v[(i, j)].powi(2);
                let u2_prev = u_prev[(i, j)].powi(2) + v_prev[(i, j)].powi(2);

                ke_current = ke_current + half * self.density * u2 * dx * dy;
                ke_prev = ke_prev + half * self.density * u2_prev * dx * dy;
            }
        }

        let ke_change = (ke_current - ke_prev).abs();
        let relative_change = if ke_prev != T::zero() {
            ke_change / ke_prev.abs()
        } else {
            T::zero()
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

impl<T: RealField + Copy + FromPrimitive> ConservationChecker<T> for EnergyConservationChecker<T> {
    type FlowField = DMatrix<T>;

    fn check_conservation(&self, field: &Self::FlowField) -> Result<ConservationReport<T>> {
        // For generic check, assume steady state temperature field
        let temperature = field;
        let temperature_prev = field.clone();
        let u = DMatrix::zeros(self.nx, self.ny);
        let v = DMatrix::zeros(self.nx, self.ny);
        let thermal_conductivity = T::from_f64(0.025).unwrap_or(T::one());
        let dt = T::from_f64(1e-3).unwrap_or(T::one());
        let dx = T::one();
        let dy = T::one();

        self.check_energy_2d(
            temperature,
            &temperature_prev,
            &u,
            &v,
            thermal_conductivity,
            dt,
            dx,
            dy,
            None,
        )
    }

    fn name(&self) -> &str {
        "Energy Conservation"
    }

    fn tolerance(&self) -> T {
        self.tolerance
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DMatrix;

    #[test]
    fn test_uniform_temperature_conservation() {
        let nx = 10;
        let ny = 10;
        let checker = EnergyConservationChecker::<f64>::new(1e-10, nx, ny, 1.0, 1000.0);

        // Uniform temperature field (steady state)
        let temperature = DMatrix::from_element(nx, ny, 300.0);
        let u = DMatrix::zeros(nx, ny);
        let v = DMatrix::zeros(nx, ny);

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
        let mut temperature = DMatrix::zeros(nx, ny);
        for i in 0..nx {
            for j in 0..ny {
                temperature[(i, j)] = 300.0 + 10.0 * (i as f64 / (nx - 1) as f64);
            }
        }

        let u = DMatrix::zeros(nx, ny);
        let v = DMatrix::zeros(nx, ny);

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
        let u = DMatrix::from_element(nx, ny, 1.0);
        let v = DMatrix::from_element(nx, ny, 0.5);

        // Check kinetic energy conservation (steady state)
        let report = checker
            .check_kinetic_energy(&u, &v, &u, &v, 0.1, 0.1)
            .expect("Kinetic energy conservation check should succeed");

        // For steady flow, kinetic energy should be conserved
        assert!(report.error < 1e-10);
        assert!(report.is_conserved);
    }
}
