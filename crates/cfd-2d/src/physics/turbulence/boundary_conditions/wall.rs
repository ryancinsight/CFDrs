//! Wall-specific boundary condition implementations for turbulence models.
//!
//! Provides wall boundary conditions for k-ε, k-ω SST, and Spalart-Allmaras models.

use super::{TurbulenceBoundaryCondition, TurbulenceBoundaryManager};
use crate::physics::turbulence::constants::{EPSILON_MIN, SST_BETA_1};
use nalgebra::RealField;
use num_traits::FromPrimitive;

impl<T: RealField + FromPrimitive + Copy> TurbulenceBoundaryManager<T> {
    /// Apply wall boundaries for k-ε model
    pub(super) fn apply_wall_boundaries_k_epsilon(
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
    pub(super) fn apply_wall_boundaries_k_omega(
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
                            let y_wall = self.wall_distances[idx]
                                .max(T::from_f64(1e-6).unwrap_or_else(T::one));
                            let omega_wall = T::from_f64(6.0).unwrap_or_else(T::one)
                                / (T::from_f64(SST_BETA_1).unwrap_or_else(T::one)
                                    * y_wall
                                    * y_wall);
                            omega[idx] = omega_wall;
                        }
                    }
                    "east" => {
                        for j in 0..self.ny {
                            let idx = j * self.nx + (self.nx - 1);
                            k[idx] = T::zero();
                            let y_wall = self.wall_distances[idx]
                                .max(T::from_f64(1e-6).unwrap_or_else(T::one));
                            let omega_wall = T::from_f64(6.0).unwrap_or_else(T::one)
                                / (T::from_f64(SST_BETA_1).unwrap_or_else(T::one)
                                    * y_wall
                                    * y_wall);
                            omega[idx] = omega_wall;
                        }
                    }
                    "south" => {
                        for i in 0..self.nx {
                            k[i] = T::zero();
                            let y_wall = self.wall_distances[i]
                                .max(T::from_f64(1e-6).unwrap_or_else(T::one));
                            let omega_wall = T::from_f64(6.0).unwrap_or_else(T::one)
                                / (T::from_f64(SST_BETA_1).unwrap_or_else(T::one)
                                    * y_wall
                                    * y_wall);
                            omega[i] = omega_wall;
                        }
                    }
                    "north" => {
                        let base_idx = (self.ny - 1) * self.nx;
                        for i in 0..self.nx {
                            let idx = base_idx + i;
                            k[idx] = T::zero();
                            let y_wall = self.wall_distances[idx]
                                .max(T::from_f64(1e-6).unwrap_or_else(T::one));
                            let omega_wall = T::from_f64(6.0).unwrap_or_else(T::one)
                                / (T::from_f64(SST_BETA_1).unwrap_or_else(T::one)
                                    * y_wall
                                    * y_wall);
                            omega[idx] = omega_wall;
                        }
                    }
                    _ => {}
                }
            }
        }
    }

    /// Apply wall boundaries for Spalart-Allmaras model
    pub(super) fn apply_wall_boundaries_sa(
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
}
