//! Momentum equation solver for SIMPLE algorithm.
//!
//! This module implements the discretization and solution of momentum equations
//! using finite volume methods with proper Rhie-Chow interpolation.

use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use crate::convection::ConvectionScheme;
use cfd_core::{BoundaryCondition, Result};
use cfd_core::constants;
use cfd_math::{SparseMatrix, SparseMatrixBuilder, LinearSolver, BiCGSTAB, LinearSolverConfig};
use nalgebra::{RealField, DVector, Vector2};
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// Momentum equation coefficients for a single cell
#[derive(Debug, Clone)]
pub struct MomentumCoefficients<T: RealField> {
    /// Central coefficient (diagonal)
    pub ap: T,
    /// East neighbor coefficient
    pub ae: T,
    /// West neighbor coefficient
    pub aw: T,
    /// North neighbor coefficient
    pub an: T,
    /// South neighbor coefficient  
    pub as_: T,
    /// Source term
    pub su: T,
}

impl<T: RealField + FromPrimitive> Default for MomentumCoefficients<T> {
    fn default() -> Self {
        Self {
            ap: T::zero(),
            ae: T::zero(),
            aw: T::zero(),
            an: T::zero(),
            as_: T::zero(),
            su: T::zero(),
        }
    }
}

/// Momentum direction for solving u or v components
#[derive(Debug, Clone, Copy)]
pub enum MomentumComponent {
    /// U-momentum (x-direction)
    U,
    /// V-momentum (y-direction)
    V,
}

/// Specialized momentum equation solver
pub struct MomentumSolver<T: RealField> {
    /// Momentum coefficients for each grid cell
    coefficients: Field2D<MomentumCoefficients<T>>,
    /// Fluid density
    density: T,
    /// Dynamic viscosity
    viscosity: T,
    /// Convection scheme strategy
    convection_scheme: Box<dyn ConvectionScheme<T>>,
}

impl<T: RealField + FromPrimitive> MomentumSolver<T> {
    /// Create new momentum solver
    pub fn new(nx: usize, ny: usize, density: T, viscosity: T, convection_scheme: Box<dyn ConvectionScheme<T>>) -> Self {
        Self {
            coefficients: Field2D::new(nx, ny, MomentumCoefficients::default()),
            density,
            viscosity,
            convection_scheme,
        }
    }

    /// Solve momentum equations for both components
    pub fn solve(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<(usize, usize), BoundaryCondition<T>>,
        dt: T,
    ) -> Result<()> {
        // Solve u-momentum
        self.solve_component(MomentumComponent::U, fields, grid, boundary_conditions, dt.clone())?;
        
        // Solve v-momentum
        self.solve_component(MomentumComponent::V, fields, grid, boundary_conditions, dt)?;
        
        Ok(())
    }

    /// Solve momentum equation for specific component
    fn solve_component(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<(usize, usize), BoundaryCondition<T>>,
        dt: T,
    ) -> Result<()> {
        // Build coefficient matrix
        self.assemble_coefficients(component, fields, grid, dt)?;
        
        // Apply boundary conditions
        self.apply_momentum_boundary_conditions(component, fields, grid, boundary_conditions)?;
        
        // Build linear system
        let (matrix, rhs) = self.build_linear_system(component, fields, grid)?;
        
        // Solve linear system
        let config = LinearSolverConfig::default();
        let solver = BiCGSTAB::new(config);
        let solution = solver.solve(&matrix, &rhs, None)?;
        
        // Update velocity field
        self.update_velocity_field(component, fields, &solution)?;
        
        Ok(())
    }

    /// Assemble discretization coefficients
    fn assemble_coefficients(
        &mut self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        dt: T,
    ) -> Result<()> {
        let dx = grid.dx.clone();
        let dy = grid.dy.clone();
        let two = T::from_f64(constants::TWO).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))?;
        
        for i in 1..fields.nx - 1 {
            for j in 1..fields.ny - 1 {
                let coeffs = self.coefficients.at_mut(i, j);
                
                // Reset coefficients
                *coeffs = MomentumCoefficients::default();
                
                // Diffusion terms
                let de = self.viscosity.clone() * dy.clone() / dx.clone();
                let dw = de.clone();
                let dn = self.viscosity.clone() * dx.clone() / dy.clone();
                let ds = dn.clone();
                
                // Convection terms - extract velocity based on component
                let (u_vel, v_vel) = match component {
                    MomentumComponent::U => {
                        let u_e = (fields.u.at(i, j).x.clone() + fields.u.at(i + 1, j).x.clone()) / two.clone();
                        let u_w = (fields.u.at(i - 1, j).x.clone() + fields.u.at(i, j).x.clone()) / two.clone();
                        let v_n = (fields.u.at(i, j).y.clone() + fields.u.at(i, j + 1).y.clone()) / two.clone();
                        let v_s = (fields.u.at(i, j - 1).y.clone() + fields.u.at(i, j).y.clone()) / two.clone();
                        ((u_e, u_w), (v_n, v_s))
                    },
                    MomentumComponent::V => {
                        let u_e = (fields.u.at(i, j).x.clone() + fields.u.at(i + 1, j).x.clone()) / two.clone();
                        let u_w = (fields.u.at(i - 1, j).x.clone() + fields.u.at(i, j).x.clone()) / two.clone();
                        let v_n = (fields.u.at(i, j).y.clone() + fields.u.at(i, j + 1).y.clone()) / two.clone();
                        let v_s = (fields.u.at(i, j - 1).y.clone() + fields.u.at(i, j).y.clone()) / two.clone();
                        ((u_e, u_w), (v_n, v_s))
                    }
                };
                
                let fe = self.density.clone() * u_vel.0 * dy.clone();
                let fw = self.density.clone() * u_vel.1 * dy.clone();
                let fn_ = self.density.clone() * v_vel.0 * dx.clone();
                let fs = self.density.clone() * v_vel.1 * dx.clone();
                
                // Apply convection scheme  
                let (ae_conv, aw_conv) = self.convection_scheme.coefficients(fe.clone(), fw.clone(), de.clone(), dw.clone());
                let (an_conv, as_conv) = self.convection_scheme.coefficients(fn_.clone(), fs.clone(), dn.clone(), ds.clone());
                
                coeffs.ae = ae_conv;
                coeffs.aw = aw_conv;
                coeffs.an = an_conv;
                coeffs.as_ = as_conv;
                
                // Time derivative (implicit Euler)
                let time_coeff = self.density.clone() * dx.clone() * dy.clone() / dt.clone();
                
                // Source term includes time derivative of old velocity
                let old_vel = match component {
                    MomentumComponent::U => fields.u.at(i, j).x.clone(),
                    MomentumComponent::V => fields.u.at(i, j).y.clone(),
                };
                coeffs.su = time_coeff * old_vel;
            }
        }
        
        Ok(())
    }

    /// Apply convection scheme to get coefficients using strategy pattern
    fn apply_convection_scheme(&self, fe: T, fw: T, de: T, dw: T) -> (T, T) {
        self.convection_scheme.coefficients(fe, fw, de, dw)
    }

    /// Apply boundary conditions for momentum equations
    fn apply_momentum_boundary_conditions(
        &mut self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
        boundary_conditions: &HashMap<(usize, usize), BoundaryCondition<T>>,
    ) -> Result<()> {
        // Implementation would apply BCs specific to momentum component
        // For brevity, implementing basic no-slip walls
        
        // Bottom and top walls (j = 0, j = ny-1)
        for i in 0..fields.nx {
            for &j in &[0, fields.ny - 1] {
                if let Some(bc) = boundary_conditions.get(&(i, j)) {
                    match bc {
                        BoundaryCondition::Wall { .. } => {
                            // No-slip: set velocity to zero
                            *fields.u_star.at_mut(i, j) = Vector2::zeros();
                        },
                        BoundaryCondition::VelocityInlet { velocity } => {
                            // Convert Vector3 to Vector2 for 2D
                            *fields.u_star.at_mut(i, j) = Vector2::new(velocity.x.clone(), velocity.y.clone());
                        },
                        BoundaryCondition::PressureOutlet { .. } => {
                            // Zero gradient
                            if j == 0 {
                                *fields.u_star.at_mut(i, j) = fields.u_star.at(i, 1).clone();
                            } else {
                                *fields.u_star.at_mut(i, j) = fields.u_star.at(i, fields.ny - 2).clone();
                            }
                        },
                        _ => {}
                    }
                }
            }
        }
        
        // Left and right walls (i = 0, i = nx-1)
        for j in 0..fields.ny {
            for &i in &[0, fields.nx - 1] {
                if let Some(bc) = boundary_conditions.get(&(i, j)) {
                    match bc {
                        BoundaryCondition::Wall { .. } => {
                            *fields.u_star.at_mut(i, j) = Vector2::zeros();
                        },
                        BoundaryCondition::VelocityInlet { velocity } => {
                            // Convert Vector3 to Vector2 for 2D
                            *fields.u_star.at_mut(i, j) = Vector2::new(velocity.x.clone(), velocity.y.clone());
                        },
                        BoundaryCondition::PressureOutlet { .. } => {
                            if i == 0 {
                                *fields.u_star.at_mut(i, j) = fields.u_star.at(1, j).clone();
                            } else {
                                *fields.u_star.at_mut(i, j) = fields.u_star.at(fields.nx - 2, j).clone();
                            }
                        },
                        _ => {}
                    }
                }
            }
        }
        
        Ok(())
    }

    /// Build linear system for momentum equation
    fn build_linear_system(
        &self,
        component: MomentumComponent,
        fields: &SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> Result<(SparseMatrix<T>, DVector<T>)> {
        let n_inner = (fields.nx - 2) * (fields.ny - 2);
        let mut matrix_builder = SparseMatrixBuilder::new(n_inner, n_inner);
        let mut rhs = DVector::zeros(n_inner);
        
        for i in 1..fields.nx - 1 {
            for j in 1..fields.ny - 1 {
                let row = (j - 1) * (fields.nx - 2) + (i - 1);
                let coeffs = self.coefficients.at(i, j);
                
                // Central coefficient
                let _ = matrix_builder.add_entry(row, row, coeffs.ap.clone());
                
                // East neighbor
                if i < fields.nx - 2 {
                    let col = (j - 1) * (fields.nx - 2) + i;
                    let _ = matrix_builder.add_entry(row, col, -coeffs.ae.clone());
                }
                
                // West neighbor
                if i > 1 {
                    let col = (j - 1) * (fields.nx - 2) + (i - 2);
                    let _ = matrix_builder.add_entry(row, col, -coeffs.aw.clone());
                }
                
                // North neighbor
                if j < fields.ny - 2 {
                    let col = j * (fields.nx - 2) + (i - 1);
                    let _ = matrix_builder.add_entry(row, col, -coeffs.an.clone());
                }
                
                // South neighbor
                if j > 1 {
                    let col = (j - 2) * (fields.nx - 2) + (i - 1);
                    let _ = matrix_builder.add_entry(row, col, -coeffs.as_.clone());
                }
                
                // Pressure gradient term
                let pressure_grad = match component {
                    MomentumComponent::U => {
                        (fields.p.at(i + 1, j).clone() - fields.p.at(i - 1, j).clone()) / T::from_f64(constants::TWO).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? / grid.dx.clone()
                    },
                    MomentumComponent::V => {
                        (fields.p.at(i, j + 1).clone() - fields.p.at(i, j - 1).clone()) / T::from_f64(constants::TWO).ok_or_else(|| cfd_core::error::Error::Numerical(cfd_core::error::NumericalErrorKind::InvalidFpOperation))? / grid.dy.clone()
                    }
                };
                
                rhs[row] = coeffs.su.clone() - pressure_grad * grid.dx.clone() * grid.dy.clone();
            }
        }
        
        Ok((matrix_builder.build()?, rhs))
    }

    /// Update velocity field with solution
    fn update_velocity_field(
        &self,
        component: MomentumComponent,
        fields: &mut SimulationFields<T>,
        solution: &DVector<T>,
    ) -> Result<()> {
        for i in 1..fields.nx - 1 {
            for j in 1..fields.ny - 1 {
                let idx = (j - 1) * (fields.nx - 2) + (i - 1);
                match component {
                    MomentumComponent::U => {
                        fields.u_star.at_mut(i, j).x = solution[idx].clone();
                    },
                    MomentumComponent::V => {
                        fields.u_star.at_mut(i, j).y = solution[idx].clone();
                    }
                }
            }
        }
        
        Ok(())
    }
}