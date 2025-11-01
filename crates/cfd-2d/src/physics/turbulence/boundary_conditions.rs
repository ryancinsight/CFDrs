//! Turbulence boundary condition implementations
//!
//! This module provides physically accurate boundary conditions for turbulence models:
//! - Wall boundary conditions for k, ε, ω, ν̃
//! - Inlet/outlet boundary conditions for turbulence quantities
//! - Wall distance calculation for wall functions

use super::constants::*;
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
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    wall_distances: Vec<T>,
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
        // Simplified wall distance calculation
        // In production code, this should use the Eikonal equation solver
        // For now, use minimum distance to domain boundaries

        for j in 0..self.ny {
            for i in 0..self.nx {
                let idx = j * self.nx + i;

                // Distance to left/right walls (from cell center)
                let dist_left = (T::from_usize(i).unwrap_or_else(T::one) + T::from_f64(0.5).unwrap_or_else(T::one)) * self.dx;
                let dist_right = (T::from_usize(self.nx - 1 - i).unwrap_or_else(T::one) + T::from_f64(0.5).unwrap_or_else(T::one)) * self.dx;

                // Distance to bottom/top walls (from cell center)
                let dist_bottom = (T::from_usize(j).unwrap_or_else(T::one) + T::from_f64(0.5).unwrap_or_else(T::one)) * self.dy;
                let dist_top = (T::from_usize(self.ny - 1 - j).unwrap_or_else(T::one) + T::from_f64(0.5).unwrap_or_else(T::one)) * self.dy;

                // Minimum distance to any wall
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
        // Apply wall boundaries (k=0, ε=minimum or wall function)
        self.apply_wall_boundaries_k_epsilon(k, epsilon, boundaries);

        // Apply inlet boundaries
        self.apply_inlet_boundaries_k_epsilon(k, epsilon, boundaries);

        // Apply outlet boundaries (zero gradient)
        self.apply_outlet_boundaries(k, epsilon, boundaries);

        // Ensure positivity
        let eps_min = T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero);
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
        // Apply wall boundaries
        self.apply_wall_boundaries_k_omega(k, omega, boundaries);

        // Apply inlet boundaries
        self.apply_inlet_boundaries_k_omega(k, omega, boundaries);

        // Apply outlet boundaries
        self.apply_outlet_boundaries(k, omega, boundaries);

        // Ensure positivity
        let omega_min = T::from_f64(OMEGA_MIN).unwrap_or_else(T::zero);
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
        // Apply wall boundaries (ν̃ = 0)
        self.apply_wall_boundaries_sa(nu_tilde, boundaries);

        // Apply inlet boundaries
        self.apply_inlet_boundaries_sa(nu_tilde, boundaries);

        // Apply outlet boundaries
        self.apply_outlet_boundaries_sa(nu_tilde, boundaries);

        // Ensure non-negative values
        for val in nu_tilde.iter_mut() {
            *val = (*val).max(T::zero());
        }
    }

    /// Apply wall boundaries for k-ε model
    fn apply_wall_boundaries_k_epsilon(
        &self,
        k: &mut [T],
        epsilon: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        let eps_min = T::from_f64(EPSILON_MIN).unwrap_or_else(T::zero);

        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Wall { .. } = bc {
                match name.as_str() {
                    "west" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx;
                            k[idx] = T::zero();
                            epsilon[idx] = eps_min;
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx + (self.nx - 1);
                            k[idx] = T::zero();
                            epsilon[idx] = eps_min;
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            k[i] = T::zero();
                            epsilon[i] = eps_min;
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            k[base_idx + i] = T::zero();
                            epsilon[base_idx + i] = eps_min;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply wall boundaries for k-ω SST model
    fn apply_wall_boundaries_k_omega(
        &self,
        k: &mut [T],
        omega: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Wall { wall_treatment: _ } = bc {
                match name.as_str() {
                    "west" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx;
                            k[idx] = T::zero();
                            // ω_wall = 6ν/(β₁ y_wall²) for SST model
                            let y_wall = self.wall_distances[idx].max(T::from_f64(1e-6).unwrap_or_else(T::one));
                            let omega_wall = T::from_f64(6.0).unwrap_or_else(T::one)
                                / (T::from_f64(SST_BETA_1).unwrap_or_else(T::one) * y_wall * y_wall);
                            omega[idx] = omega_wall;
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx + (self.nx - 1);
                            k[idx] = T::zero();
                            let y_wall = self.wall_distances[idx].max(T::from_f64(1e-6).unwrap_or_else(T::one));
                            let omega_wall = T::from_f64(6.0).unwrap_or_else(T::one)
                                / (T::from_f64(SST_BETA_1).unwrap_or_else(T::one) * y_wall * y_wall);
                            omega[idx] = omega_wall;
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            k[i] = T::zero();
                            let y_wall = self.wall_distances[i].max(T::from_f64(1e-6).unwrap_or_else(T::one));
                            let omega_wall = T::from_f64(6.0).unwrap_or_else(T::one)
                                / (T::from_f64(SST_BETA_1).unwrap_or_else(T::one) * y_wall * y_wall);
                            omega[i] = omega_wall;
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            let idx = base_idx + i;
                            k[idx] = T::zero();
                            let y_wall = self.wall_distances[idx].max(T::from_f64(1e-6).unwrap_or_else(T::one));
                            let omega_wall = T::from_f64(6.0).unwrap_or_else(T::one)
                                / (T::from_f64(SST_BETA_1).unwrap_or_else(T::one) * y_wall * y_wall);
                            omega[idx] = omega_wall;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply wall boundaries for Spalart-Allmaras model
    fn apply_wall_boundaries_sa(
        &self,
        nu_tilde: &mut [T],
        boundaries: &[(String, TurbulenceBoundaryCondition<T>)],
    ) {
        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Wall { .. } = bc {
                match name.as_str() {
                    "west" => {
                        for j in 0..self.ny {
                            nu_tilde[j * self.nx] = T::zero();
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            nu_tilde[j * self.nx + (self.nx - 1)] = T::zero();
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            nu_tilde[i] = T::zero();
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            nu_tilde[base_idx + i] = T::zero();
                        }
                    }
                    _ => {}
                }
            }
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
                let k_inlet = T::from_f64(1.5).unwrap_or_else(T::one)
                    * *turbulence_intensity * *turbulence_intensity
                    * *reference_velocity * *reference_velocity;

                // ε = C_μ^{3/4} * k^{3/2} / l
                let c_mu_34 = T::from_f64(C_MU).unwrap_or_else(T::one).powf(T::from_f64(0.75).unwrap_or_else(T::one));
                let k_32 = k_inlet.powf(T::from_f64(1.5).unwrap_or_else(T::one));
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
                let k_inlet = T::from_f64(1.5).unwrap_or_else(T::one)
                    * *turbulence_intensity * *turbulence_intensity
                    * *reference_velocity * *reference_velocity;

                // ω = √k / (C_μ^{1/4} * l)
                let c_mu_14 = T::from_f64(C_MU).unwrap_or_else(T::one).powf(T::from_f64(0.25).unwrap_or_else(T::one));
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
                let factor = T::from_f64(1.5).unwrap_or_else(T::one)
                    * *turbulence_intensity * *turbulence_intensity
                    * *reference_velocity * *turbulence_length_scale;
                let c_mu_inv = T::from_f64(1.0 / C_MU).unwrap_or_else(T::one);
                let nu_tilde_inlet = factor * c_mu_inv.powf(T::from_f64(0.75).unwrap_or_else(T::one));

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
    fn apply_outlet_boundaries(&self, field1: &mut [T], field2: &mut [T], boundaries: &[(String, TurbulenceBoundaryCondition<T>)]) {
        for (name, bc) in boundaries {
            if let TurbulenceBoundaryCondition::Outlet = bc {
                match name.as_str() {
                    "west" => {
                        // Zero gradient: field[i,j] = field[i+1,j]
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
                        // Zero gradient: field[i,j] = field[i-1,j]
                        for j in 0..self.ny {
                            let idx = j * self.nx + (self.nx - 1);
                            let idx_prev = idx - 1;
                            field1[idx] = field1[idx_prev];
                            field2[idx] = field2[idx_prev];
                        }
                    }
                    "south" => {
                        // Zero gradient: field[i,j] = field[i,j+1]
                        for i in 0..self.nx {
                            let idx_next = i + self.nx;
                            if idx_next < field1.len() {
                                field1[i] = field1[idx_next];
                                field2[i] = field2[idx_next];
                            }
                        }
                    }
                    "north" => {
                        // Zero gradient: field[i,j] = field[i,j-1]
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
    fn apply_outlet_boundaries_sa(&self, field: &mut [T], boundaries: &[(String, TurbulenceBoundaryCondition<T>)]) {
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
    use approx::assert_relative_eq;

    #[test]
    fn test_wall_distance_calculation() {
        let manager = TurbulenceBoundaryManager::<f64>::new(5, 5, 0.1, 0.1);

        // Corner should have minimum distance
        assert!(manager.wall_distances[0] > 0.0);
        assert!(manager.wall_distances[0] < 0.1); // Less than one cell

        // Center should have maximum distance
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

        // South wall (bottom) should have k=0, epsilon=minimum
        for i in 0..5 {
            assert_eq!(k[i], 0.0);
            assert!(epsilon[i] > 0.0 && epsilon[i] < 1e-6); // epsilon_min
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

        // South wall should have k=0, omega=large value (wall function)
        for i in 0..5 {
            assert_eq!(k[i], 0.0);
            assert!(omega[i] > 1.0); // Large omega near wall
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

        // South wall should have nu_tilde=0
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
                turbulence_intensity: 0.05, // 5%
                turbulence_length_scale: 0.01,
                reference_velocity: 1.0,
            },
        )];

        let mut k = vec![0.0; 25];
        let mut epsilon = vec![0.0; 25];

        manager.apply_k_epsilon_boundaries(&mut k, &mut epsilon, &boundaries);

        // West boundary should have non-zero turbulence
        for j in 0..5 {
            let idx = j * 5;
            assert!(k[idx] > 0.0);
            assert!(epsilon[idx] > 0.0);
        }
    }
}
