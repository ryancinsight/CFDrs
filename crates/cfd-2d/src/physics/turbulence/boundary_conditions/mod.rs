//! Turbulence boundary condition implementations
//!
//! This module provides physically accurate boundary conditions for turbulence models:
//! - Wall boundary conditions for k, ε, ω, ν̃
//! - Inlet/outlet boundary conditions for turbulence quantities
//! - Wall distance calculation for wall functions
//!
//! # Theorem
//! The turbulence model must satisfy the realizability conditions for the Reynolds stress tensor.
//!
//! **Proof sketch**:
//! For any turbulent flow, the Reynolds stress tensor $\tau_{ij} = -\rho \overline{u_i^\prime u_j^\prime}$
//! must be positive semi-definite. This requires that the turbulent kinetic energy $k \ge 0$
//! and the normal stresses $\overline{u_i^\prime u_i^\prime} \ge 0$. The implemented model
//! enforces these constraints either through exact transport equations or bounded eddy-viscosity
//! formulations, ensuring physical realizability and numerical stability.

mod wall;

use super::constants::{C_MU, EPSILON_MIN, OMEGA_MIN};
use super::wall_functions::WallTreatment;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Turbulence boundary condition types
#[derive(Debug, Clone)]
pub enum TurbulenceBoundaryCondition<T: RealField + Copy> {
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
pub struct TurbulenceBoundaryManager<T: RealField + Copy> {
    pub(crate) nx: usize,
    pub(crate) ny: usize,
    dx: T,
    dy: T,
    pub(crate) wall_distances: Vec<T>,
}

impl<T: RealField + FromPrimitive + Copy> TurbulenceBoundaryManager<T> {
    /// Create a new boundary condition manager
    pub fn new(nx: usize, ny: usize, dx: T, dy: T) -> Self {
        let mut manager = Self {
            nx,
            ny,
            dx,
            dy,
            wall_distances: vec![T::zero(); nx * ny],
        };
        manager.calculate_wall_distances();
        manager
    }

    /// Calculate wall distances for all grid points
    fn calculate_wall_distances(&mut self) {
        for j in 0..self.ny {
            for i in 0..self.nx {
                let idx = j * self.nx + i;

                let dist_left = (T::from_usize(i).expect("analytical constant conversion")
                    + T::from_f64(0.5).expect("analytical constant conversion"))
                    * self.dx;
                let dist_right = (T::from_usize(self.nx - 1 - i)
                    .expect("analytical constant conversion")
                    + T::from_f64(0.5).expect("analytical constant conversion"))
                    * self.dx;
                let dist_bottom = (T::from_usize(j).expect("analytical constant conversion")
                    + T::from_f64(0.5).expect("analytical constant conversion"))
                    * self.dy;
                let dist_top = (T::from_usize(self.ny - 1 - j)
                    .expect("analytical constant conversion")
                    + T::from_f64(0.5).expect("analytical constant conversion"))
                    * self.dy;

                let min_dist = dist_left.min(dist_right).min(dist_bottom).min(dist_top);
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

        let eps_min = T::from_f64(EPSILON_MIN).expect("analytical constant conversion");
        for i in 0..k.len() {
            k[i] = k[i].max(T::zero());
            epsilon[i] = epsilon[i].max(eps_min);
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

        let omega_min = T::from_f64(OMEGA_MIN).expect("analytical constant conversion");
        for i in 0..k.len() {
            k[i] = k[i].max(T::zero());
            omega[i] = omega[i].max(omega_min);
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
            *val = (*val).max(T::zero());
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
                let k_inlet = T::from_f64(1.5).expect("analytical constant conversion")
                    * *turbulence_intensity
                    * *turbulence_intensity
                    * *reference_velocity
                    * *reference_velocity;

                // ε = C_μ^{3/4} * k^{3/2} / l
                let c_mu_34 = T::from_f64(C_MU)
                    .expect("analytical constant conversion")
                    .powf(T::from_f64(0.75).expect("analytical constant conversion"));
                let k_32 = k_inlet.powf(T::from_f64(1.5).expect("analytical constant conversion"));
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
                let k_inlet = T::from_f64(1.5).expect("analytical constant conversion")
                    * *turbulence_intensity
                    * *turbulence_intensity
                    * *reference_velocity
                    * *reference_velocity;

                // ω = √k / (C_μ^{1/4} * l)
                let c_mu_14 = T::from_f64(C_MU)
                    .expect("analytical constant conversion")
                    .powf(T::from_f64(0.25).expect("analytical constant conversion"));
                let omega_inlet = k_inlet.sqrt() / (c_mu_14 * *turbulence_length_scale);

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
                let factor = T::from_f64(1.5).expect("analytical constant conversion")
                    * *turbulence_intensity
                    * *turbulence_intensity
                    * *reference_velocity
                    * *turbulence_length_scale;
                let c_mu_inv = T::from_f64(1.0 / C_MU).expect("analytical constant conversion");
                let nu_tilde_inlet = factor
                    * c_mu_inv.powf(T::from_f64(0.75).expect("analytical constant conversion"));

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
    use crate::physics::turbulence::wall_functions::WallFunction;

    #[test]
    fn test_wall_distance_calculation() {
        let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);

        assert!(manager.wall_distances[0] > 0.0);
        assert!(manager.wall_distances[0] < 0.1);

        let center_idx = 2 * 5 + 2;
        assert!(manager.wall_distances[center_idx] > manager.wall_distances[0]);
    }

    #[test]
    fn test_k_epsilon_wall_boundaries() {
        let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);
        let boundaries = vec![(
            "south".to_string(),
            TurbulenceBoundaryCondition::Wall {
                wall_treatment: WallTreatment::new(WallFunction::Standard),
            },
        )];

        let mut k = vec![1.0; 25];
        let mut epsilon = vec![1.0; 25];

        manager.apply_k_epsilon_boundaries(&mut k, &mut epsilon, &boundaries);

        for i in 0..5 {
            assert_eq!(k[i], 0.0);
            assert!(epsilon[i] > 0.0 && epsilon[i] < 1e-6);
        }
    }

    #[test]
    fn test_k_omega_wall_boundaries() {
        let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);
        let boundaries = vec![(
            "south".to_string(),
            TurbulenceBoundaryCondition::Wall {
                wall_treatment: WallTreatment::new(WallFunction::Standard),
            },
        )];

        let mut k = vec![1.0; 25];
        let mut omega = vec![1.0; 25];

        manager.apply_k_omega_sst_boundaries(&mut k, &mut omega, &boundaries);

        for i in 0..5 {
            assert_eq!(k[i], 0.0);
            assert!(omega[i] > 1.0);
        }
    }

    #[test]
    fn test_sa_wall_boundaries() {
        let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);
        let boundaries = vec![(
            "south".to_string(),
            TurbulenceBoundaryCondition::Wall {
                wall_treatment: WallTreatment::new(WallFunction::Standard),
            },
        )];

        let mut nu_tilde = vec![1.0; 25];

        manager.apply_spalart_allmaras_boundaries(&mut nu_tilde, &boundaries);

        for i in 0..5 {
            assert_eq!(nu_tilde[i], 0.0);
        }
    }

    #[test]
    fn test_inlet_boundaries() {
        let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);
        let boundaries = vec![(
            "west".to_string(),
            TurbulenceBoundaryCondition::Inlet {
                turbulence_intensity: 0.05,
                turbulence_length_scale: 0.01,
                reference_velocity: 1.0,
            },
        )];

        let mut k = vec![0.0; 25];
        let mut epsilon = vec![0.0; 25];

        manager.apply_k_epsilon_boundaries(&mut k, &mut epsilon, &boundaries);

        for j in 0..5 {
            let idx = j * 5;
            assert!(k[idx] > 0.0);
            assert!(epsilon[idx] > 0.0);
        }
    }
}
