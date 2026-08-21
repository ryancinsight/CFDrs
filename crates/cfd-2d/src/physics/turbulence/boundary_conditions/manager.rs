use crate::physics::turbulence::constants::{C_MU, EPSILON_MIN, OMEGA_MIN};
use crate::physics::turbulence::wall_functions::WallTreatment;
use cfd_core::error::Error;
use eunomia::{FloatElement, NumericElement, RealField};

/// Turbulence boundary condition types
#[derive(Debug, Clone)]
pub enum TurbulenceBoundaryCondition<T: RealField> {
    /// Wall boundary with specified wall function treatment
    Wall {
        /// Wall function treatment for near-wall turbulence modeling
        wall_treatment: WallTreatment<T>,
    },
    /// Inlet with specified turbulence intensity and length scale
    Inlet {
        /// Turbulence intensity at inlet (dimensionless)
        turbulence_intensity: T,
        /// Turbulent length scale at inlet
        turbulence_length_scale: T,
        /// Reference velocity for turbulence scaling
        reference_velocity: T,
    },
    /// Outlet with zero gradient (natural boundary condition)
    Outlet,
    /// Periodic boundary condition
    Periodic,
}

/// Turbulence boundary condition manager
pub struct TurbulenceBoundaryManager<T: RealField> {
    pub(crate) nx: usize,
    pub(crate) ny: usize,
    dx: T,
    dy: T,
    pub(crate) wall_distances: Vec<T>,
}

impl<T: RealField> TurbulenceBoundaryManager<T> {
    /// Create a new boundary condition manager.
    ///
    /// # Panics
    /// Panics if any invariant is violated (see [`Self::try_new`]).
    #[must_use]
    pub fn new(nx: usize, ny: usize, dx: T, dy: T) -> Self {
        Self::try_new(nx, ny, dx, dy).unwrap_or_else(|error| {
            panic!("TurbulenceBoundaryManager::new called with invalid arguments: {error}");
        })
    }

    /// Create a new boundary condition manager with invariant validation.
    ///
    /// # Errors
    /// Returns `Error::InvalidConfiguration` if `nx == 0`, `ny == 0`,
    /// `dx` non-finite or non-positive, or `dy` non-finite or non-positive.
    pub fn try_new(nx: usize, ny: usize, dx: T, dy: T) -> cfd_core::error::Result<Self> {
        if nx == 0 {
            return Err(Error::InvalidConfiguration(
                "TurbulenceBoundaryManager::try_new: nx must be at least 1".to_string(),
            ));
        }
        if ny == 0 {
            return Err(Error::InvalidConfiguration(
                "TurbulenceBoundaryManager::try_new: ny must be at least 1".to_string(),
            ));
        }
        if !<T as NumericElement>::is_finite(dx) || dx <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidConfiguration(format!(
                "TurbulenceBoundaryManager::try_new: dx must be finite and positive, got {dx:?}"
            )));
        }
        if !<T as NumericElement>::is_finite(dy) || dy <= <T as NumericElement>::ZERO {
            return Err(Error::InvalidConfiguration(format!(
                "TurbulenceBoundaryManager::try_new: dy must be finite and positive, got {dy:?}"
            )));
        }
        let mut manager = Self {
            nx,
            ny,
            dx,
            dy,
            wall_distances: vec![T::ZERO; nx * ny],
        };
        manager.calculate_wall_distances();
        Ok(manager)
    }

    /// Calculate wall distances for all grid points
    fn calculate_wall_distances(&mut self) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let idx = j * self.nx + i;

                let half = T::from_f64(0.5);
                let dist_left = (T::from_f64(i as f64) + half) * self.dx;
                let dist_right = (T::from_f64((self.nx - 1 - i) as f64) + half) * self.dx;
                let dist_bottom = (T::from_f64(j as f64) + half) * self.dy;
                let dist_top = (T::from_f64((self.ny - 1 - j) as f64) + half) * self.dy;

                let min_dist = dist_left
                    .min_scalar(dist_right)
                    .min_scalar(dist_bottom)
                    .min_scalar(dist_top);
                self.wall_distances[idx] = min_dist;
            }
        }
    }

    /// Apply k-ε boundary conditions
    pub fn apply_k_epsilon_boundaries(
        &self,
        k: &mut [T],
        epsilon: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        self.apply_wall_boundaries_k_epsilon(k, epsilon, boundaries);
        self.apply_inlet_boundaries_k_epsilon(k, epsilon, boundaries);
        self.apply_outlet_boundaries(k, epsilon, boundaries);

        let eps_min = T::from_f64(EPSILON_MIN);
        for i in 0..k.len() {
            k[i] = k[i].max_scalar(T::ZERO);
            epsilon[i] = epsilon[i].max_scalar(eps_min);
        }
    }

    /// Apply k-ω SST boundary conditions
    pub fn apply_k_omega_sst_boundaries(
        &self,
        k: &mut [T],
        omega: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        self.apply_wall_boundaries_k_omega(k, omega, boundaries);
        self.apply_inlet_boundaries_k_omega(k, omega, boundaries);
        self.apply_outlet_boundaries(k, omega, boundaries);

        let omega_min = T::from_f64(OMEGA_MIN);
        for i in 0..k.len() {
            k[i] = k[i].max_scalar(T::ZERO);
            omega[i] = omega[i].max_scalar(omega_min);
        }
    }

    /// Apply Spalart-Allmaras boundary conditions
    pub fn apply_spalart_allmaras_boundaries(
        &self,
        nu_tilde: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        self.apply_wall_boundaries_sa(nu_tilde, boundaries);
        self.apply_inlet_boundaries_sa(nu_tilde, boundaries);
        self.apply_outlet_boundaries_sa(nu_tilde, boundaries);

        for val in nu_tilde.iter_mut() {
            *val = (*val).max_scalar(T::ZERO);
        }
    }

    /// Apply inlet boundaries for k-ε model
    fn apply_inlet_boundaries_k_epsilon(
        &self,
        k: &mut [T],
        epsilon: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Inlet {
                turbulence_intensity,
                turbulence_length_scale,
                reference_velocity,
            } = bc
            {
                // k = (3/2) * (I * U_ref)²
                let k_inlet = T::from_f64(1.5)
                    * *turbulence_intensity
                    * *turbulence_intensity
                    * *reference_velocity
                    * *reference_velocity;

                // ε = C_μ^{3/4} * k^{3/2} / l
                let c_mu_34 = <T as FloatElement>::powf(T::from_f64(C_MU), T::from_f64(0.75));
                let k_32 = <T as FloatElement>::powf(k_inlet, T::from_f64(1.5));
                let eps_inlet = c_mu_34 * k_32 / *turbulence_length_scale;

                match name.as_str() {
                    "west" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx;
                            k[idx] = k_inlet;
                            epsilon[idx] = eps_inlet;
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx + (self.nx - 1);
                            k[idx] = k_inlet;
                            epsilon[idx] = eps_inlet;
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            k[i] = k_inlet;
                            epsilon[i] = eps_inlet;
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            k[base_idx + i] = k_inlet;
                            epsilon[base_idx + i] = eps_inlet;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply inlet boundaries for k-ω SST model
    fn apply_inlet_boundaries_k_omega(
        &self,
        k: &mut [T],
        omega: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Inlet {
                turbulence_intensity,
                turbulence_length_scale,
                reference_velocity,
            } = bc
            {
                // k = (3/2) * (I * U_ref)²
                let k_inlet = T::from_f64(1.5)
                    * *turbulence_intensity
                    * *turbulence_intensity
                    * *reference_velocity
                    * *reference_velocity;

                // ω = √k / (C_μ^{1/4} * l)
                let c_mu_14 = <T as FloatElement>::powf(T::from_f64(C_MU), T::from_f64(0.25));
                let omega_inlet =
                    <T as NumericElement>::sqrt(k_inlet) / (c_mu_14 * *turbulence_length_scale);

                match name.as_str() {
                    "west" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx;
                            k[idx] = k_inlet;
                            omega[idx] = omega_inlet;
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx + (self.nx - 1);
                            k[idx] = k_inlet;
                            omega[idx] = omega_inlet;
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            k[i] = k_inlet;
                            omega[i] = omega_inlet;
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            k[base_idx + i] = k_inlet;
                            omega[base_idx + i] = omega_inlet;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply inlet boundaries for Spalart-Allmaras model
    fn apply_inlet_boundaries_sa(
        &self,
        nu_tilde: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Inlet {
                turbulence_intensity,
                turbulence_length_scale,
                reference_velocity,
            } = bc
            {
                // For SA model, ν̃_inlet ≈ (3/2) * I² * U_ref * l / C_μ^{3/4}
                let factor = T::from_f64(1.5)
                    * *turbulence_intensity
                    * *turbulence_intensity
                    * *reference_velocity
                    * *turbulence_length_scale;
                let c_mu_inv = T::from_f64(1.0 / C_MU);
                let nu_tilde_inlet =
                    factor * <T as FloatElement>::powf(c_mu_inv, T::from_f64(0.75));

                match name.as_str() {
                    "west" => {
                        for j in 0..self.ny {
                            nu_tilde[j * self.nx] = nu_tilde_inlet;
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            nu_tilde[j * self.nx + (self.nx - 1)] = nu_tilde_inlet;
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            nu_tilde[i] = nu_tilde_inlet;
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            nu_tilde[base_idx + i] = nu_tilde_inlet;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply outlet boundaries (zero gradient)
    fn apply_outlet_boundaries(
        &self,
        field1: &mut [T],
        field2: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Outlet = bc {
                match name.as_str() {
                    "west" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx;
                            let idx_next = idx + 1;
                            if idx_next < field1.len() {
                                field1[idx] = field1[idx_next];
                                field2[idx] = field2[idx_next];
                            }
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx + (self.nx - 1);
                            let idx_prev = idx - 1;
                            field1[idx] = field1[idx_prev];
                            field2[idx] = field2[idx_prev];
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            let idx_next = i + self.nx;
                            if idx_next < field1.len() {
                                field1[i] = field1[idx_next];
                                field2[i] = field2[idx_next];
                            }
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            let idx = base_idx + i;
                            let idx_prev = idx - self.nx;
                            field1[idx] = field1[idx_prev];
                            field2[idx] = field2[idx_prev];
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply outlet boundaries for single field (SA model)
    fn apply_outlet_boundaries_sa(
        &self,
        field: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Outlet = bc {
                match name.as_str() {
                    "west" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx;
                            let idx_next = idx + 1;
                            if idx_next < field.len() {
                                field[idx] = field[idx_next];
                            }
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx + (self.nx - 1);
                            let idx_prev = idx - 1;
                            field[idx] = field[idx_prev];
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            let idx_next = i + self.nx;
                            if idx_next < field.len() {
                                field[i] = field[idx_next];
                            }
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            let idx = base_idx + i;
                            let idx_prev = idx - self.nx;
                            field[idx] = field[idx_prev];
                        }
                    }
                    _ => {}
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    /// **Positive**: `try_new` accepts valid arguments.
    #[test]
    fn turbulence_boundary_manager_try_new_accepts_valid_arguments() {
        let mgr = TurbulenceBoundaryManager::<f64>::try_new(8, 12, 1e-3, 1e-3)
            .expect("valid must succeed");
        assert_eq!(mgr.nx, 8);
        assert_eq!(mgr.ny, 12);
        assert_eq!(mgr.wall_distances.len(), 96);
    }

    /// **Adversarial**: zero `nx` is rejected.
    #[test]
    fn turbulence_boundary_manager_try_new_rejects_zero_nx() {
        match TurbulenceBoundaryManager::<f64>::try_new(0, 10, 1.0, 1.0) {
            Err(e) => assert!(e.to_string().contains("nx"), "error must mention nx: {e}"),
            Ok(_) => panic!("zero nx must be rejected"),
        }
    }

    /// **Adversarial**: zero `dy` is rejected.
    #[test]
    fn turbulence_boundary_manager_try_new_rejects_zero_dy() {
        match TurbulenceBoundaryManager::<f64>::try_new(8, 12, 1.0, 0.0) {
            Err(e) => assert!(e.to_string().contains("dy"), "error must mention dy: {e}"),
            Ok(_) => panic!("zero dy must be rejected"),
        }
    }

    /// **Boundary**: `new` panics on invalid `dx` (thin wrapper contract).
    #[test]
    #[should_panic(expected = "dx")]
    fn turbulence_boundary_manager_new_panics_on_invalid_dx() {
        let _ = TurbulenceBoundaryManager::<f64>::new(8, 12, 0.0, 1.0);
    }
}
