//! Update and trait-impl methods for `SpalartAllmaras`.

use super::helpers::wall_distance_field_2d;
use super::SpalartAllmaras;
use cfd_core::{error::Result, physics::constants::mathematical::numeric::TWO};
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use tracing::instrument;

impl<T: RealField + FromPrimitive + Copy> SpalartAllmaras<T> {
    /// Update modified turbulent viscosity field with explicit wall distances
    ///
    /// This method allows providing custom wall distances (e.g., for DES where d is replaced by min(d, Cdes*Delta))
    #[instrument(skip(
        self,
        nu_tilde,
        velocity,
        molecular_viscosity,
        wall_distances,
        dx,
        dy,
        dt
    ))]
    pub fn update_with_distance(
        &self,
        nu_tilde: &mut [T],
        velocity: &[Vector2<T>],
        molecular_viscosity: T,
        wall_distances: &[T],
        dx: T,
        dy: T,
        dt: T,
    ) -> Result<()> {
        // Temporary storage for new values
        let mut nu_new = vec![T::zero(); nu_tilde.len()];

        for j in 1..self.ny - 1 {
            for i in 1..self.nx - 1 {
                let idx = j * self.nx + i;

                // Compute velocity gradients
                let du_dx = (velocity[idx + 1].x - velocity[idx - 1].x)
                    / (T::from_f64(TWO).expect("analytical constant conversion") * dx);
                let du_dy = (velocity[idx + self.nx].x - velocity[idx - self.nx].x)
                    / (T::from_f64(TWO).expect("analytical constant conversion") * dy);
                let dv_dx = (velocity[idx + 1].y - velocity[idx - 1].y)
                    / (T::from_f64(TWO).expect("analytical constant conversion") * dx);
                let dv_dy = (velocity[idx + self.nx].y - velocity[idx - self.nx].y)
                    / (T::from_f64(TWO).expect("analytical constant conversion") * dy);

                let velocity_gradient = [[du_dx, du_dy], [dv_dx, dv_dy]];

                // Calculate vorticity
                let vorticity = self.vorticity_magnitude(&velocity_gradient);

                // Calculate modified vorticity
                let s_tilde = self.modified_vorticity(
                    vorticity,
                    nu_tilde[idx],
                    molecular_viscosity,
                    wall_distances[idx],
                );

                // Production term
                let production = self.production(nu_tilde[idx], s_tilde);

                // Wall destruction function
                let fw =
                    self.wall_destruction_function(nu_tilde[idx], s_tilde, wall_distances[idx]);

                // Destruction term
                let destruction = self.destruction(nu_tilde[idx], wall_distances[idx], fw);

                // Diffusion term: (1/σ)∇·[(ν+ν̃)∇ν̃]
                let nu_total_e = molecular_viscosity + nu_tilde[idx + 1];
                let nu_total_w = molecular_viscosity + nu_tilde[idx - 1];
                let nu_total_n = molecular_viscosity + nu_tilde[idx + self.nx];
                let nu_total_s = molecular_viscosity + nu_tilde[idx - self.nx];

                let dnu_dx_e = (nu_tilde[idx + 1] - nu_tilde[idx]) / dx;
                let dnu_dx_w = (nu_tilde[idx] - nu_tilde[idx - 1]) / dx;
                let dnu_dy_n = (nu_tilde[idx + self.nx] - nu_tilde[idx]) / dy;
                let dnu_dy_s = (nu_tilde[idx] - nu_tilde[idx - self.nx]) / dy;

                let diffusion_x = (nu_total_e * dnu_dx_e - nu_total_w * dnu_dx_w) / dx;
                let diffusion_y = (nu_total_n * dnu_dy_n - nu_total_s * dnu_dy_s) / dy;
                let diffusion = (diffusion_x + diffusion_y) / self.sigma;

                // Cross-diffusion term: (Cb2/σ)|∇ν̃|²
                let dnu_dx = (nu_tilde[idx + 1] - nu_tilde[idx - 1])
                    / (T::from_f64(TWO).expect("analytical constant conversion") * dx);
                let dnu_dy = (nu_tilde[idx + self.nx] - nu_tilde[idx - self.nx])
                    / (T::from_f64(TWO).expect("analytical constant conversion") * dy);
                let grad_nu_sq = dnu_dx * dnu_dx + dnu_dy * dnu_dy;
                let cross_diffusion = (self.cb2 / self.sigma) * grad_nu_sq;

                // Trip term (zero for fully turbulent)
                let trip = self.trip_term(nu_tilde[idx], wall_distances[idx]);

                // Explicit time integration
                nu_new[idx] = nu_tilde[idx]
                    + dt * (production - destruction + diffusion + cross_diffusion + trip);

                // Ensure positive values
                nu_new[idx] = nu_new[idx].max(T::zero());
            }
        }

        // Update field with efficient slice copy
        nu_tilde.copy_from_slice(&nu_new);

        // Apply boundary conditions
        self.apply_boundary_conditions(nu_tilde);

        Ok(())
    }

    /// Update modified turbulent viscosity field
    ///
    /// Solves: ∂ν̃/∂t = P - D + (1/σ)∇·[(ν+ν̃)∇ν̃] + (Cb2/σ)|∇ν̃|²
    #[instrument(skip(self, nu_tilde, velocity, molecular_viscosity, dx, dy, dt))]
    pub fn update(
        &self,
        nu_tilde: &mut [T],
        velocity: &[Vector2<T>],
        molecular_viscosity: T,
        dx: T,
        dy: T,
        dt: T,
    ) -> Result<()> {
        // Calculate wall distances
        let wall_distances = wall_distance_field_2d(self.nx, self.ny, dx, dy);

        // Update with calculated wall distances
        self.update_with_distance(
            nu_tilde,
            velocity,
            molecular_viscosity,
            &wall_distances,
            dx,
            dy,
            dt,
        )
    }
}

// Implement TurbulenceModel trait for Spalart-Allmaras
impl<T: RealField + FromPrimitive + Copy + num_traits::ToPrimitive>
    crate::physics::turbulence::TurbulenceModel<T> for SpalartAllmaras<T>
{
    fn turbulent_viscosity(&self, _k: T, epsilon_or_omega: T, density: T) -> T {
        // For SA model, k is not used, epsilon_or_omega represents ν̃ (modified viscosity)
        let nu_tilde = epsilon_or_omega;
        let molecular_viscosity = T::from_f64(1e-5).expect("analytical constant conversion"); // Typical air viscosity
        density * self.eddy_viscosity(nu_tilde, molecular_viscosity)
    }

    fn production_term(
        &self,
        velocity_gradient: &[[T; 2]; 2],
        _turbulent_viscosity: T,
        turbulence_variable: T,
        wall_distance: T,
        molecular_viscosity: T,
    ) -> T {
        // For SA model, production is calculated as P = Cb1 * S̃ * ν̃
        // turbulence_variable corresponds to ν̃ (nu_tilde)
        let nu_tilde = turbulence_variable;
        let vorticity = self.vorticity_magnitude(velocity_gradient);

        let s_tilde =
            self.modified_vorticity(vorticity, nu_tilde, molecular_viscosity, wall_distance);

        self.cb1 * s_tilde * nu_tilde
    }

    fn dissipation_term(&self, nu_tilde: T, wall_distance: T) -> T {
        use crate::physics::turbulence::constants::EPSILON_MIN;
        let wall_distance =
            wall_distance.max(T::from_f64(EPSILON_MIN).expect("analytical constant conversion"));
        let ratio = nu_tilde / wall_distance;
        self.cw1 * ratio * ratio
    }

    fn update(
        &mut self,
        k: &mut [T],
        epsilon_or_omega: &mut [T],
        velocity: &[Vector2<T>],
        _density: T,
        molecular_viscosity: T,
        dt: T,
        dx: T,
        dy: T,
    ) -> Result<()> {
        // For SA model, k is not used, epsilon_or_omega represents ν̃ (modified viscosity)
        // Update the SA transport equation
        SpalartAllmaras::update(
            self,
            epsilon_or_omega,
            velocity,
            molecular_viscosity,
            dx,
            dy,
            dt,
        )?;

        // k is not updated in SA model (single equation model)
        // Set k to zero or some reference value if needed
        for k_val in k.iter_mut() {
            *k_val = T::zero();
        }

        Ok(())
    }

    fn name(&self) -> &'static str {
        "Spalart-Allmaras"
    }

    fn is_valid_for_reynolds(&self, reynolds: T) -> bool {
        // SA model is valid for moderate to high Reynolds numbers
        reynolds > T::from_f64(1e4).expect("analytical constant conversion")
    }
}
