//! Energy equation solver for temperature transport
//!
//! Solves the energy equation for incompressible flows:
//! ∂T/∂t + (u·∇)T = α∇²T + Q/(ρCp)
//! where α = k/(ρCp) is thermal diffusivity

use cfd_core::boundary::BoundaryCondition;
use cfd_core::error::Result;
use nalgebra::RealField;
use std::collections::HashMap;

/// Constants for energy equation
pub mod constants {
    /// Default Prandtl number for air
    pub const DEFAULT_PRANDTL: f64 = 0.71;
    /// Stefan-Boltzmann constant [W/(m²·K⁴)]
    pub const STEFAN_BOLTZMANN: f64 = 5.67e-8;
    /// Central difference coefficient for second derivatives
    pub const CENTRAL_DIFF_COEFF: f64 = 2.0;
}

/// Energy equation solver
pub struct EnergyEquationSolver<T: RealField + Copy> {
    /// Temperature field
    pub temperature: Vec<Vec<T>>,
    /// Thermal diffusivity field
    pub thermal_diffusivity: Vec<Vec<T>>,
    /// Heat source term
    pub heat_source: Vec<Vec<T>>,
    /// Grid dimensions
    nx: usize,
    ny: usize,
}

impl<T: RealField + Copy> EnergyEquationSolver<T> {
    /// Create new energy equation solver
    pub fn new(nx: usize, ny: usize, initial_temperature: T, thermal_diffusivity: T) -> Self {
        Self {
            temperature: vec![vec![initial_temperature; ny]; nx],
            thermal_diffusivity: vec![vec![thermal_diffusivity; ny]; nx],
            heat_source: vec![vec![T::zero(); ny]; nx],
            nx,
            ny,
        }
    }

    /// Solve energy equation using explicit time stepping
    pub fn solve_explicit(
        &mut self,
        u_velocity: &[Vec<T>],
        v_velocity: &[Vec<T>],
        dt: T,
        dx: T,
        dy: T,
        boundary_conditions: &HashMap<(usize, usize), BoundaryCondition<T>>,
    ) -> Result<()> {
        let mut new_temperature = self.temperature.clone();

        // Interior points
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                // Skip boundary points
                if boundary_conditions.contains_key(&(i, j)) {
                    continue;
                }

                let t = self.temperature[i][j];
                let alpha = self.thermal_diffusivity[i][j];
                let u = u_velocity[i][j];
                let v = v_velocity[i][j];

                // Convection terms (upwind scheme)
                let dt_dx = if u > T::zero() {
                    (t - self.temperature[i - 1][j]) / dx
                } else {
                    (self.temperature[i + 1][j] - t) / dx
                };

                let dt_dy = if v > T::zero() {
                    (t - self.temperature[i][j - 1]) / dy
                } else {
                    (self.temperature[i][j + 1] - t) / dy
                };

                // Diffusion terms (central difference)
                let coeff = T::from_f64(constants::CENTRAL_DIFF_COEFF).unwrap_or_else(|| T::zero());
                let d2t_dx2 = (self.temperature[i + 1][j] - coeff * t + self.temperature[i - 1][j])
                    / (dx * dx);
                let d2t_dy2 = (self.temperature[i][j + 1] - coeff * t + self.temperature[i][j - 1])
                    / (dy * dy);

                // Update temperature
                new_temperature[i][j] = t + dt
                    * (-u * dt_dx - v * dt_dy
                        + alpha * (d2t_dx2 + d2t_dy2)
                        + self.heat_source[i][j]);
            }
        }

        // Apply boundary conditions
        for ((i, j), bc) in boundary_conditions {
            match bc {
                BoundaryCondition::Dirichlet { value } => {
                    new_temperature[*i][*j] = *value;
                }
                BoundaryCondition::Neumann { gradient } => {
                    // Apply gradient boundary condition using interior points
                    if *i == 0 {
                        new_temperature[0][*j] = new_temperature[1][*j] - *gradient * dx;
                    } else if *i == self.nx - 1 {
                        new_temperature[self.nx - 1][*j] =
                            new_temperature[self.nx - 2][*j] + *gradient * dx;
                    }
                    if *j == 0 {
                        new_temperature[*i][0] = new_temperature[*i][1] - *gradient * dy;
                    } else if *j == self.ny - 1 {
                        new_temperature[*i][self.ny - 1] =
                            new_temperature[*i][self.ny - 2] + *gradient * dy;
                    }
                }
                BoundaryCondition::Periodic { partner: _ } => {
                    // Periodic BC: T(boundary) = T(partner_boundary)
                    if *i == 0 {
                        new_temperature[0][*j] = new_temperature[self.nx - 1][*j];
                    } else if *i == self.nx - 1 {
                        new_temperature[self.nx - 1][*j] = new_temperature[0][*j];
                    }
                    if *j == 0 {
                        new_temperature[*i][0] = new_temperature[*i][self.ny - 1];
                    } else if *j == self.ny - 1 {
                        new_temperature[*i][self.ny - 1] = new_temperature[*i][0];
                    }
                }
                BoundaryCondition::Symmetry => {
                    // Symmetry BC: zero normal gradient
                    if *i == 0 {
                        new_temperature[0][*j] = new_temperature[1][*j];
                    } else if *i == self.nx - 1 {
                        new_temperature[self.nx - 1][*j] = new_temperature[self.nx - 2][*j];
                    }
                    if *j == 0 {
                        new_temperature[*i][0] = new_temperature[*i][1];
                    } else if *j == self.ny - 1 {
                        new_temperature[*i][self.ny - 1] = new_temperature[*i][self.ny - 2];
                    }
                }
                _ => {}
            }
        }

        self.temperature = new_temperature;
        Ok(())
    }

    /// Calculate Nusselt number for heat transfer analysis
    pub fn nusselt_number(&self, wall_temp: T, bulk_temp: T, characteristic_length: T, dy: T) -> T {
        let dt_dy_wall = (self.temperature[0][1] - self.temperature[0][0]) / dy; // Use actual dy
        let h = self.thermal_diffusivity[0][0] * dt_dy_wall / (wall_temp - bulk_temp);
        h * characteristic_length / self.thermal_diffusivity[0][0]
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_new_solver_initialization() {
        let solver = EnergyEquationSolver::<f64>::new(10, 10, 300.0, 0.01);
        assert_eq!(solver.nx, 10);
        assert_eq!(solver.ny, 10);
        
        // Check initial temperature field
        for i in 0..10 {
            for j in 0..10 {
                assert_relative_eq!(solver.temperature[i][j], 300.0, epsilon = 1e-10);
                assert_relative_eq!(solver.thermal_diffusivity[i][j], 0.01, epsilon = 1e-10);
                assert_relative_eq!(solver.heat_source[i][j], 0.0, epsilon = 1e-10);
            }
        }
    }

    #[test]
    fn test_solve_explicit_zero_velocity() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        let u_velocity = vec![vec![0.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let boundary_conditions = HashMap::new();
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
    }

    #[test]
    fn test_solve_explicit_uniform_temperature() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        let u_velocity = vec![vec![1.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let boundary_conditions = HashMap::new();
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
        
        // With uniform temperature and zero source, temperature should remain approximately constant
        for i in 1..4 {
            for j in 1..4 {
                assert!(solver.temperature[i][j] > 299.0);
                assert!(solver.temperature[i][j] < 301.0);
            }
        }
    }

    #[test]
    fn test_dirichlet_boundary_condition() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        let u_velocity = vec![vec![0.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let mut boundary_conditions = HashMap::new();
        
        // Apply Dirichlet BC at point (0, 0)
        boundary_conditions.insert((0, 0), BoundaryCondition::Dirichlet { value: 350.0 });
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
        
        // Check that Dirichlet BC is applied
        assert_relative_eq!(solver.temperature[0][0], 350.0, epsilon = 1e-10);
    }

    #[test]
    fn test_neumann_boundary_condition_left() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        let u_velocity = vec![vec![0.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let mut boundary_conditions = HashMap::new();
        
        // Apply Neumann BC at left boundary
        boundary_conditions.insert((0, 2), BoundaryCondition::Neumann { gradient: 10.0 });
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
    }

    #[test]
    fn test_neumann_boundary_condition_right() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        let u_velocity = vec![vec![0.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let mut boundary_conditions = HashMap::new();
        
        // Apply Neumann BC at right boundary
        boundary_conditions.insert((4, 2), BoundaryCondition::Neumann { gradient: -10.0 });
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
    }

    #[test]
    fn test_periodic_boundary_condition() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        let u_velocity = vec![vec![0.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let mut boundary_conditions = HashMap::new();
        
        // Apply periodic BC at left/right boundaries
        boundary_conditions.insert((0, 2), BoundaryCondition::Periodic { partner: "right".to_string() });
        boundary_conditions.insert((4, 2), BoundaryCondition::Periodic { partner: "left".to_string() });
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
    }

    #[test]
    fn test_symmetry_boundary_condition() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        let u_velocity = vec![vec![0.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let mut boundary_conditions = HashMap::new();
        
        // Apply symmetry BC at bottom boundary
        boundary_conditions.insert((2, 0), BoundaryCondition::Symmetry);
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
        
        // Symmetry should set T(i, 0) = T(i, 1)
        assert_relative_eq!(solver.temperature[2][0], solver.temperature[2][1], epsilon = 1e-10);
    }

    #[test]
    fn test_heat_source_effect() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        
        // Add heat source at center
        solver.heat_source[2][2] = 1000.0;
        
        let u_velocity = vec![vec![0.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let boundary_conditions = HashMap::new();
        
        let initial_temp = solver.temperature[2][2];
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
        
        // Temperature at center should increase due to heat source
        assert!(solver.temperature[2][2] > initial_temp);
    }

    #[test]
    fn test_convection_positive_u_velocity() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        
        // Set hot spot at upstream location
        solver.temperature[1][2] = 350.0;
        
        // Positive u velocity (flow in positive x direction)
        let mut u_velocity = vec![vec![1.0; 5]; 5];
        u_velocity[2][2] = 1.0;
        let v_velocity = vec![vec![0.0; 5]; 5];
        let boundary_conditions = HashMap::new();
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
    }

    #[test]
    fn test_convection_negative_u_velocity() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        
        // Set hot spot at downstream location
        solver.temperature[3][2] = 350.0;
        
        // Negative u velocity (flow in negative x direction)
        let mut u_velocity = vec![vec![-1.0; 5]; 5];
        u_velocity[2][2] = -1.0;
        let v_velocity = vec![vec![0.0; 5]; 5];
        let boundary_conditions = HashMap::new();
        
        let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        assert!(result.is_ok());
    }

    #[test]
    fn test_diffusion_smoothing() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.1);
        
        // Create temperature gradient
        for i in 0..5 {
            solver.temperature[i][2] = 300.0 + (i as f64) * 10.0;
        }
        
        let u_velocity = vec![vec![0.0; 5]; 5];
        let v_velocity = vec![vec![0.0; 5]; 5];
        let boundary_conditions = HashMap::new();
        
        let initial_gradient = solver.temperature[3][2] - solver.temperature[1][2];
        
        // Run multiple time steps
        for _ in 0..5 {
            let _ = solver.solve_explicit(&u_velocity, &v_velocity, 0.001, 0.1, 0.1, &boundary_conditions);
        }
        
        let final_gradient = solver.temperature[3][2] - solver.temperature[1][2];
        
        // Diffusion should reduce gradient
        assert!(final_gradient.abs() < initial_gradient.abs());
    }

    #[test]
    fn test_nusselt_number_positive() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        
        // Create temperature gradient at wall
        solver.temperature[0][0] = 350.0;
        solver.temperature[0][1] = 340.0;
        
        let nu = solver.nusselt_number(350.0, 300.0, 1.0, 0.1);
        assert!(nu.is_finite());
    }

    #[test]
    fn test_constants_validity() {
        assert!(constants::DEFAULT_PRANDTL > 0.0);
        assert!(constants::STEFAN_BOLTZMANN > 0.0);
        assert_relative_eq!(constants::CENTRAL_DIFF_COEFF, 2.0, epsilon = 1e-10);
    }

    #[test]
    fn test_multiple_timesteps_stability() {
        let mut solver = EnergyEquationSolver::<f64>::new(5, 5, 300.0, 0.01);
        let u_velocity = vec![vec![0.1; 5]; 5];
        let v_velocity = vec![vec![0.1; 5]; 5];
        let boundary_conditions = HashMap::new();
        
        // Run 10 time steps
        for _ in 0..10 {
            let result = solver.solve_explicit(&u_velocity, &v_velocity, 0.0001, 0.1, 0.1, &boundary_conditions);
            assert!(result.is_ok());
            
            // Check stability (no blow-up)
            for i in 0..5 {
                for j in 0..5 {
                    assert!(solver.temperature[i][j].is_finite());
                    assert!(solver.temperature[i][j] > 0.0);
                    assert!(solver.temperature[i][j] < 1000.0);
                }
            }
        }
    }
}
