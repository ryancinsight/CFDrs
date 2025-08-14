//! Pressure correction solver for SIMPLE algorithm.
//!
//! This module implements the pressure correction equation using finite volume
//! discretization with proper handling of boundary conditions.

use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use cfd_core::{BoundaryCondition, Result};
use cfd_core::constants;
use cfd_math::{SparseMatrix, SparseMatrixBuilder, LinearSolver, ConjugateGradient, LinearSolverConfig};
use nalgebra::{RealField, DVector};
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// Pressure correction equation solver
pub struct PressureCorrector<T: RealField> {
    /// Diagonal coefficients for pressure correction
    coefficients: Field2D<T>,
    /// Fluid density
    density: T,
    /// Pressure under-relaxation factor
    alpha_p: T,
}

impl<T: RealField + FromPrimitive> PressureCorrector<T> {
    /// Create new pressure corrector
    pub fn new(nx: usize, ny: usize, density: T, alpha_p: T) -> Self {
        Self {
            coefficients: Field2D::new(nx, ny, T::zero()),
            density,
            alpha_p,
        }
    }

    /// Solve pressure correction equation
    pub fn solve(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<(usize, usize), BoundaryCondition<T>>,
    ) -> Result<()> {
        // Build coefficient matrix and RHS
        let (matrix, rhs) = self.build_pressure_system(fields, grid)?;
        
        // Apply boundary conditions
        self.apply_pressure_boundary_conditions(fields, boundary_conditions)?;
        
        // Solve linear system
        let config = LinearSolverConfig::default();
        let solver = ConjugateGradient::new(config);
        let solution = solver.solve(&matrix, &rhs, None)?;
        
        // Update pressure correction field
        self.update_pressure_correction(fields, &solution)?;
        
        // Correct velocity and pressure fields
        self.correct_fields(fields, grid)?;
        
        Ok(())
    }

    /// Build pressure correction linear system
    fn build_pressure_system(
        &mut self,
        fields: &SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let dx = grid.dx.clone();
        let dy = grid.dy.clone();
        let n_inner = (fields.nx - 2) * (fields.ny - 2);
        let mut matrix_builder = SparseMatrixBuilder::new(n_inner, n_inner);
        let mut rhs = DVector::zeros(n_inner);
        
        // Calculate coefficients for each interior cell
        for i in 1..fields.nx - 1 {
            for j in 1..fields.ny - 1 {
                let row = (j - 1) * (fields.nx - 2) + (i - 1);
                
                // Coefficient based on momentum equation diagonal
                // For simplified implementation, use inverse of cell volume
                let d_u = dy.clone() / (dx.clone() * self.density.clone());
                let d_v = dx.clone() / (dy.clone() * self.density.clone());
                
                // Pressure correction coefficients
                let ae = d_u.clone() * dy.clone();
                let aw = d_u.clone() * dy.clone();
                let an = d_v.clone() * dx.clone();
                let as_ = d_v.clone() * dx.clone();
                let ap = ae.clone() + aw.clone() + an.clone() + as_.clone();
                
                // Store diagonal coefficient for velocity correction
                *self.coefficients.at_mut(i, j) = ap.clone();
                
                // Matrix entries
                matrix_builder.add_entry(row, row, ap);
                
                // East neighbor
                if i < fields.nx - 2 {
                    let col = (j - 1) * (fields.nx - 2) + i;
                    matrix_builder.add_entry(row, col, -ae);
                }
                
                // West neighbor
                if i > 1 {
                    let col = (j - 1) * (fields.nx - 2) + (i - 2);
                    matrix_builder.add_entry(row, col, -aw);
                }
                
                // North neighbor
                if j < fields.ny - 2 {
                    let col = j * (fields.nx - 2) + (i - 1);
                    matrix_builder.add_entry(row, col, -an);
                }
                
                // South neighbor
                if j > 1 {
                    let col = (j - 2) * (fields.nx - 2) + (i - 1);
                    matrix_builder.add_entry(row, col, -as_);
                }
                
                // RHS: negative of mass imbalance
                let mass_imbalance = self.calculate_mass_imbalance(fields, i, j, dx.clone(), dy.clone());
                rhs[row] = -mass_imbalance;
            }
        }
        
        Ok((matrix_builder.build()?, rhs))
    }

    /// Calculate mass imbalance for continuity equation
    fn calculate_mass_imbalance(
        &self,
        fields: &SimulationFields<T>,
        i: usize,
        j: usize,
        dx: T,
        dy: T,
    ) -> T {
        // Mass flux through cell faces
        let flux_e = self.density.clone() * fields.u_star.at(i + 1, j).x.clone() * dy.clone();
        let flux_w = self.density.clone() * fields.u_star.at(i, j).x.clone() * dy.clone();
        let flux_n = self.density.clone() * fields.u_star.at(i, j + 1).y.clone() * dx.clone();
        let flux_s = self.density.clone() * fields.u_star.at(i, j).y.clone() * dx.clone();
        
        // Net mass imbalance (should be zero for continuity)
        flux_e - flux_w + flux_n - flux_s
    }

    /// Apply boundary conditions for pressure correction
    fn apply_pressure_boundary_conditions(
        &self,
        fields: &mut SimulationFields<T>,
        boundary_conditions: &HashMap<(usize, usize), BoundaryCondition<T>>,
    ) -> Result<()> {
        // Set pressure correction to zero at boundaries (default)
        // This ensures that pressure correction doesn't affect pressure level
        
        // Bottom and top boundaries
        for i in 0..fields.nx {
            *fields.p_prime.at_mut(i, 0) = T::zero();
            *fields.p_prime.at_mut(i, fields.ny - 1) = T::zero();
        }
        
        // Left and right boundaries
        for j in 0..fields.ny {
            *fields.p_prime.at_mut(0, j) = T::zero();
            *fields.p_prime.at_mut(fields.nx - 1, j) = T::zero();
        }
        
        // Apply specific boundary conditions if needed
        for ((i, j), bc) in boundary_conditions {
            match bc {
                BoundaryCondition::PressureOutlet { pressure } => {
                    // For outlets, pressure correction can be set based on specified pressure
                    let current_pressure = fields.p.at(*i, *j).clone();
                    *fields.p_prime.at_mut(*i, *j) = pressure.clone() - current_pressure;
                },
                _ => {
                    // For other BCs, keep zero pressure correction
                    *fields.p_prime.at_mut(*i, *j) = T::zero();
                }
            }
        }
        
        Ok(())
    }

    /// Update pressure correction field with solution
    fn update_pressure_correction(
        &self,
        fields: &mut SimulationFields<T>,
        solution: &DVector<T>,
    ) -> Result<()> {
        for i in 1..fields.nx - 1 {
            for j in 1..fields.ny - 1 {
                let idx = (j - 1) * (fields.nx - 2) + (i - 1);
                *fields.p_prime.at_mut(i, j) = solution[idx].clone();
            }
        }
        
        Ok(())
    }

    /// Correct velocity and pressure fields using pressure correction
    fn correct_fields(
        &self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> Result<()> {
        let dx = grid.dx.clone();
        let dy = grid.dy.clone();
        let two = T::from_f64(constants::TWO).unwrap();
        
        // Correct velocity field
        for i in 1..fields.nx - 1 {
            for j in 1..fields.ny - 1 {
                // Get velocity correction based on pressure gradient
                let dp_dx = (fields.p_prime.at(i + 1, j).clone() - fields.p_prime.at(i - 1, j).clone()) / (two.clone() * dx.clone());
                let dp_dy = (fields.p_prime.at(i, j + 1).clone() - fields.p_prime.at(i, j - 1).clone()) / (two.clone() * dy.clone());
                
                // Apply velocity correction (simplified - should use momentum diagonal)
                let d_coeff = T::from_f64(0.1).unwrap(); // Simplified coefficient
                fields.u.at_mut(i, j).x = fields.u_star.at(i, j).x.clone() - d_coeff.clone() * dp_dx / self.density.clone();
                fields.u.at_mut(i, j).y = fields.u_star.at(i, j).y.clone() - d_coeff * dp_dy / self.density.clone();
            }
        }
        
        // Correct pressure field with under-relaxation
        for i in 0..fields.nx {
            for j in 0..fields.ny {
                let correction = self.alpha_p.clone() * fields.p_prime.at(i, j).clone();
                *fields.p.at_mut(i, j) = fields.p.at(i, j).clone() + correction;
            }
        }
        
        Ok(())
    }

    /// Calculate maximum residual for convergence checking
    pub fn calculate_residual(&self, fields: &SimulationFields<T>, grid: &StructuredGrid2D<T>) -> T {
        let dx = grid.dx.clone();
        let dy = grid.dy.clone();
        let mut max_residual = T::zero();
        
        for i in 1..fields.nx - 1 {
            for j in 1..fields.ny - 1 {
                let residual = self.calculate_mass_imbalance(fields, i, j, dx.clone(), dy.clone()).abs();
                if residual > max_residual {
                    max_residual = residual;
                }
            }
        }
        
        max_residual
    }
}