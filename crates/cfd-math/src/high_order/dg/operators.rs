//! Core operators for Discontinuous Galerkin methods.
//!
//! This module implements the spatial discretization operators for DG methods,
//! including volume and surface integrals, weak form operators, and boundary conditions.

use super::*;
use nalgebra::{DVector, DMatrix, DMatrixSlice, DVectorSlice, DMatrixSliceMut};
use std::f64::consts::PI;
use std::ops::{Add, Sub, Mul, Div};

/// Type of boundary condition
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryCondition {
    /// Dirichlet boundary condition (fixed value)
    Dirichlet(f64),
    /// Neumann boundary condition (fixed derivative)
    Neumann(f64),
    /// Periodic boundary condition
    Periodic,
    /// Outflow boundary condition
    Outflow,
    /// Reflective boundary condition
    Reflective,
}

/// Type of numerical flux for boundary conditions
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryFlux {
    /// Central flux
    Central,
    /// Upwind flux
    Upwind,
    /// Lax-Friedrichs flux
    LaxFriedrichs,
    /// HLL flux
    HLL,
}

/// Parameters for DG operators
#[derive(Debug, Clone)]
pub struct DGOperatorParams {
    /// Type of numerical flux for volume integrals
    pub volume_flux: FluxType,
    /// Type of numerical flux for surface integrals
    pub surface_flux: FluxType,
    /// Type of limiter
    pub limiter: LimiterType,
    /// Parameters for the limiter
    pub limiter_params: LimiterParams,
    /// Parameters for the flux
    pub flux_params: FluxParams,
    /// Whether to use strong form (true) or weak form (false)
    pub strong_form: bool,
    /// Whether to use the divergence form (true) or the non-conservative form (false)
    pub divergence_form: bool,
    /// Whether to use the local Lax-Friedrichs flux
    pub use_lax_friedrichs: bool,
    /// Lax-Friedrichs parameter
    pub alpha: f64,
    /// CFL number for time stepping
    pub cfl: f64,
    /// Tolerance for time stepping
    pub tolerance: f64,
    /// Maximum number of iterations
    pub max_iter: usize,
}

impl Default for DGOperatorParams {
    fn default() -> Self {
        Self {
            volume_flux: FluxType::Central,
            surface_flux: FluxType::LaxFriedrichs,
            limiter: LimiterType::Minmod,
            limiter_params: LimiterParams::default(),
            flux_params: FluxParams::default(),
            strong_form: true,
            divergence_form: true,
            use_lax_friedrichs: true,
            alpha: 1.0,
            cfl: 0.1,
            tolerance: 1e-10,
            max_iter: 1000,
        }
    }
}

impl DGOperatorParams {
    /// Create a new set of DG operator parameters
    pub fn new() -> Self {
        Self::default()
    }
    
    /// Set the volume flux type
    pub fn with_volume_flux(mut self, flux_type: FluxType) -> Self {
        self.volume_flux = flux_type;
        self
    }
    
    /// Set the surface flux type
    pub fn with_surface_flux(mut self, flux_type: FluxType) -> Self {
        self.surface_flux = flux_type;
        self
    }
    
    /// Set the limiter type
    pub fn with_limiter(mut self, limiter: LimiterType) -> Self {
        self.limiter = limiter;
        self
    }
    
    /// Set the limiter parameters
    pub fn with_limiter_params(mut self, params: LimiterParams) -> Self {
        self.limiter_params = params;
        self
    }
    
    /// Set the flux parameters
    pub fn with_flux_params(mut self, params: FluxParams) -> Self {
        self.flux_params = params;
        self
    }
    
    /// Set whether to use the strong form
    pub fn with_strong_form(mut self, strong_form: bool) -> Self {
        self.strong_form = strong_form;
        self
    }
    
    /// Set whether to use the divergence form
    pub fn with_divergence_form(mut self, divergence_form: bool) -> Self {
        self.divergence_form = divergence_form;
        self
    }
    
    /// Set whether to use the Lax-Friedrichs flux
    pub fn with_lax_friedrichs(mut self, use_lax_friedrichs: bool) -> Self {
        self.use_lax_friedrichs = use_lax_friedrichs;
        self
    }
    
    /// Set the Lax-Friedrichs parameter
    pub fn with_alpha(mut self, alpha: f64) -> Self {
        self.alpha = alpha;
        self
    }
    
    /// Set the CFL number
    pub fn with_cfl(mut self, cfl: f64) -> Self {
        self.cfl = cfl;
        self
    }
    
    /// Set the tolerance
    pub fn with_tolerance(mut self, tolerance: f64) -> Self {
        self.tolerance = tolerance;
        self
    }
    
    /// Set the maximum number of iterations
    pub fn with_max_iter(mut self, max_iter: usize) -> Self {
        self.max_iter = max_iter;
        self
    }
}

/// Represents a DG operator for spatial discretization
pub struct DGOperator {
    /// Polynomial order
    pub order: usize,
    /// Number of components (for systems of equations)
    pub num_components: usize,
    /// Number of quadrature points
    pub num_quad: usize,
    /// Quadrature points
    pub quad_points: DVector<f64>,
    /// Quadrature weights
    pub quad_weights: DVector<f64>,
    /// Basis functions at quadrature points (num_basis × num_quad)
    pub phi: DMatrix<f64>,
    /// Basis function derivatives at quadrature points (num_basis × num_quad)
    pub dphi_dx: DMatrix<f64>,
    /// Mass matrix (num_basis × num_basis)
    pub mass_matrix: DMatrix<f64>,
    /// Stiffness matrix (num_basis × num_basis)
    pub stiffness_matrix: DMatrix<f64>,
    /// Differentiation matrix (num_basis × num_basis)
    pub diff_matrix: DMatrix<f64>,
    /// Lifting operators for boundary conditions
    pub lift_operators: Vec<DMatrix<f64>>,
    /// Parameters for the DG operator
    pub params: DGOperatorParams,
}

impl DGOperator {
    /// Create a new DG operator
    pub fn new(
        order: usize,
        num_components: usize,
        params: Option<DGOperatorParams>,
    ) -> Result<Self> {
        if order == 0 {
            return Err(DGError::InvalidOrder(order));
        }
        
        let params = params.unwrap_or_default();
        let num_basis = order + 1;
        let num_quad = order + 1; // Use (order + 1) points for exact integration
        
        // Compute quadrature points and weights
        let (quad_points, quad_weights) = gauss_lobatto_quadrature(num_quad)?;
        
        // Initialize basis functions and derivatives
        let mut phi = DMatrix::zeros(num_basis, num_quad);
        let mut dphi_dx = DMatrix::zeros(num_basis, num_quad);
        
        // Compute basis functions and derivatives at quadrature points
        for i in 0..num_basis {
            for (q, &xq) in quad_points.iter().enumerate() {
                phi[(i, q)] = legendre_poly(i, xq);
                dphi_dx[(i, q)] = legendre_poly_deriv(i, xq);
            }
        }
        
        // Compute mass matrix M_ij = ∫ φ_i(x) φ_j(x) dx
        let mut mass_matrix = DMatrix::zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                let mut m_ij = 0.0;
                for q in 0..num_quad {
                    m_ij += quad_weights[q] * phi[(i, q)] * phi[(j, q)];
                }
                mass_matrix[(i, j)] = m_ij;
            }
        }
        
        // Compute stiffness matrix K_ij = ∫ φ_i'(x) φ_j(x) dx
        let mut stiffness_matrix = DMatrix::zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                let mut k_ij = 0.0;
                for q in 0..num_quad {
                    k_ij += quad_weights[q] * dphi_dx[(i, q)] * phi[(j, q)];
                }
                stiffness_matrix[(i, j)] = k_ij;
            }
        }
        
        // Compute differentiation matrix D_ij = φ_j'(x_i)
        let mut diff_matrix = DMatrix::zeros(num_basis, num_basis);
        for i in 0..num_basis {
            for j in 0..num_basis {
                diff_matrix[(i, j)] = legendre_poly_deriv(j, quad_points[i]);
            }
        }
        
        // Compute lifting operators for boundary conditions
        let mut lift_operators = Vec::with_capacity(2);
        
        // Left boundary lifting operator
        let mut lift_l = DMatrix::zeros(num_basis, 1);
        for i in 0..num_basis {
            lift_l[(i, 0)] = phi[(i, 0)] * quad_weights[0];
        }
        lift_operators.push(lift_l);
        
        // Right boundary lifting operator
        let mut lift_r = DMatrix::zeros(num_basis, 1);
        for i in 0..num_basis {
            lift_r[(i, 0)] = phi[(i, num_quad - 1)] * quad_weights[num_quad - 1];
        }
        lift_operators.push(lift_r);
        
        Ok(Self {
            order,
            num_components,
            num_quad,
            quad_points,
            quad_weights,
            phi,
            dphi_dx,
            mass_matrix,
            stiffness_matrix,
            diff_matrix,
            lift_operators,
            params,
        })
    }
    
    /// Compute the right-hand side of the DG semi-discrete equation
    pub fn rhs<F, G>(
        &self,
        u: &DMatrix<f64>,
        flux: F,
        boundary_condition: G,
    ) -> Result<DMatrix<f64>>
    where
        F: Fn(&DVector<f64>) -> DVector<f64>,
        G: Fn(f64, &DVector<f64>, bool) -> DVector<f64>,
    {
        let num_basis = self.order + 1;
        let mut rhs = DMatrix::zeros(self.num_components, num_basis);
        
        // Compute volume integral terms
        if self.params.strong_form {
            // Strong form: M^{-1} * (K^T * F - F^*|_{-1}^1)
            
            // Compute the flux at all quadrature points
            let mut f_quad = DMatrix::zeros(self.num_components, self.num_quad);
            for q in 0..self.num_quad {
                // Evaluate the solution at the quadrature point
                let mut u_q = DVector::zeros(self.num_components);
                for i in 0..num_basis {
                    u_q += u.column(i) * self.phi[(i, q)];
                }
                
                // Compute the flux
                f_quad.column_mut(q).copy_from(&flux(&u_q));
            }
            
            // Compute K^T * F
            let mut kt_f = DMatrix::zeros(self.num_components, num_basis);
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    let mut sum = 0.0;
                    for q in 0..self.num_quad {
                        sum += self.quad_weights[q] * self.dphi_dx[(j, q)] * f_quad[(i, q)];
                    }
                    kt_f[(i, j)] = sum;
                }
            }
            
            // Compute the numerical flux at the boundaries
            let u_l = u.column(0).into_owned();
            let u_r = u.column(num_basis - 1).into_owned();
            
            // Apply boundary conditions
            let f_l = boundary_condition(-1.0, &u_l, false);
            let f_r = boundary_condition(1.0, &u_r, true);
            
            // Compute the numerical flux
            let f_num_l = self.compute_numerical_flux(-1.0, &f_l, &u_l, None, &boundary_condition);
            let f_num_r = self.compute_numerical_flux(1.0, &f_r, &u_r, None, &boundary_condition);
            
            // Compute the boundary terms
            let mut boundary_term = DMatrix::zeros(self.num_components, num_basis);
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    boundary_term[(i, j)] = 
                        self.phi[(j, self.num_quad - 1)] * f_num_r[i] - 
                        self.phi[(j, 0)] * f_num_l[i];
                }
            }
            
            // Compute the right-hand side
            rhs = kt_f - boundary_term;
            
            // Multiply by the inverse mass matrix
            let mass_lu = self.mass_matrix.lu();
            for i in 0..self.num_components {
                let mut rhs_row = rhs.row_mut(i);
                let rhs_slice = rhs_row.as_mut_slice();
                let rhs_vec = DVector::from_vec(rhs_slice.to_vec());
                let sol = mass_lu.solve(&rhs_vec);
                rhs_slice.copy_from_slice(sol.as_slice());
            }
        } else {
            // Weak form: F * dφ/dx - F^* φ|_{-1}^1
            
            // Compute the flux at all quadrature points
            let mut f_quad = DMatrix::zeros(self.num_components, self.num_quad);
            for q in 0..self.num_quad {
                // Evaluate the solution at the quadrature point
                let mut u_q = DVector::zeros(self.num_components);
                for i in 0..num_basis {
                    u_q += u.column(i) * self.phi[(i, q)];
                }
                
                // Compute the flux
                f_quad.column_mut(q).copy_from(&flux(&u_q));
            }
            
            // Compute the volume integral: F * dφ/dx
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    let mut sum = 0.0;
                    for q in 0..self.num_quad {
                        sum += self.quad_weights[q] * f_quad[(i, q)] * self.dphi_dx[(j, q)];
                    }
                    rhs[(i, j)] = -sum; // Negative sign from integration by parts
                }
            }
            
            // Compute the numerical flux at the boundaries
            let u_l = u.column(0).into_owned();
            let u_r = u.column(num_basis - 1).into_owned();
            
            // Apply boundary conditions
            let f_l = boundary_condition(-1.0, &u_l, false);
            let f_r = boundary_condition(1.0, &u_r, true);
            
            // Compute the numerical flux
            let f_num_l = self.compute_numerical_flux(-1.0, &f_l, &u_l, None, &boundary_condition);
            let f_num_r = self.compute_numerical_flux(1.0, &f_r, &u_r, None, &boundary_condition);
            
            // Add the surface integral terms
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    rhs[(i, j)] +=
                        self.phi[(j, self.num_quad - 1)] * f_num_r[i] - 
                        self.phi[(j, 0)] * f_num_l[i];
                }
            }
            
            // Multiply by the inverse mass matrix
            let mass_lu = self.mass_matrix.lu();
            for i in 0..self.num_components {
                let mut rhs_row = rhs.row_mut(i);
                let rhs_slice = rhs_row.as_mut_slice();
                let rhs_vec = DVector::from_vec(rhs_slice.to_vec());
                let sol = mass_lu.solve(&rhs_vec);
                rhs_slice.copy_from_slice(sol.as_slice());
            }
        }
        
        Ok(rhs)
    }
    
    /// Compute the numerical flux at a boundary
    fn compute_numerical_flux<G>(
        &self,
        x: f64,
        f: &DVector<f64>,
        u: &DVector<f64>,
        u_ext: Option<&DVector<f64>>,
        boundary_condition: G,
    ) -> DVector<f64>
    where
        G: Fn(f64, &DVector<f64>, bool) -> DVector<f64>,
    {
        let is_right_boundary = x > 0.0;
        
        // Get the external state (from boundary condition or neighboring element)
        let u_ext = match u_ext {
            Some(u) => u.clone(),
            None => boundary_condition(x, u, is_right_boundary),
        };
        
        // Compute the external flux
        let f_ext = boundary_condition(x, &u_ext, is_right_boundary);
        
        // Compute the numerical flux based on the selected type
        match self.params.surface_flux {
            FluxType::Central => 0.5 * (f + &f_ext),
            FluxType::Upwind => {
                // Simple upwinding based on characteristic speeds
                let a = 1.0; // Characteristic speed (simplified)
                if a >= 0.0 {
                    f.clone()
                } else {
                    f_ext
                }
            }
            FluxType::LaxFriedrichs => {
                let alpha = self.params.alpha;
                0.5 * (f + &f_ext) - 0.5 * alpha * (&u_ext - u)
            }
            _ => {
                // Default to Lax-Friedrichs
                let alpha = self.params.alpha;
                0.5 * (f + &f_ext) - 0.5 * alpha * (&u_ext - u)
            }
        }
    }
    
    /// Apply boundary conditions to the solution
    pub fn apply_boundary_conditions<G>(
        &self,
        u: &mut DMatrix<f64>,
        boundary_condition: G,
    ) where
        G: Fn(f64, &DVector<f64>, bool) -> DVector<f64>,
    {
        let num_basis = self.order + 1;
        
        // Apply boundary conditions at the left boundary (x = -1)
        let u_l = u.column(0).into_owned();
        let u_l_bc = boundary_condition(-1.0, &u_l, false);
        
        // Apply boundary conditions at the right boundary (x = 1)
        let u_r = u.column(num_basis - 1).into_owned();
        let u_r_bc = boundary_condition(1.0, &u_r, true);
        
        // Update the solution with the boundary values
        for i in 0..self.num_components {
            u[(i, 0)] = u_l_bc[i];
            u[(i, num_basis - 1)] = u_r_bc[i];
        }
    }
    
    /// Project a function onto the DG basis
    pub fn project<F>(
        &self,
        f: F,
    ) -> DMatrix<f64>
    where
        F: Fn(f64) -> DVector<f64>,
    {
        let num_basis = self.order + 1;
        let mut u = DMatrix::zeros(self.num_components, num_basis);
        
        // Compute the projection
        for q in 0..self.num_quad {
            let x_q = self.quad_points[q];
            let w_q = self.quad_weights[q];
            let f_q = f(x_q);
            
            for i in 0..num_basis {
                let phi_i = self.phi[(i, q)];
                for c in 0..self.num_components {
                    u[(c, i)] += w_q * f_q[c] * phi_i;
                }
            }
        }
        
        // Solve the mass matrix system
        let mass_lu = self.mass_matrix.lu();
        for i in 0..self.num_components {
            let mut u_row = u.row_mut(i);
            let u_slice = u_row.as_mut_slice();
            let u_vec = DVector::from_vec(u_slice.to_vec());
            let sol = mass_lu.solve(&u_vec);
            u_slice.copy_from_slice(sol.as_slice());
        }
        
        u
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_dg_operator_creation() {
        let order = 2;
        let num_components = 1;
        let params = DGOperatorParams::new()
            .with_volume_flux(FluxType::Central)
            .with_surface_flux(FluxType::LaxFriedrichs)
            .with_limiter(LimiterType::Minmod);
        
        let dg_op = DGOperator::new(order, num_components, Some(params.clone()));
        assert!(dg_op.is_ok());
        
        let dg_op = dg_op.unwrap();
        assert_eq!(dg_op.order, order);
        assert_eq!(dg_op.num_components, num_components);
        assert_eq!(dg_op.quad_points.len(), order + 1);
        assert_eq!(dg_op.quad_weights.len(), order + 1);
        assert_eq!(dg_op.phi.nrows(), order + 1);
        assert_eq!(dg_op.phi.ncols(), order + 1);
        assert_eq!(dg_op.mass_matrix.nrows(), order + 1);
        assert_eq!(dg_op.mass_matrix.ncols(), order + 1);
        assert_eq!(dg_op.stiffness_matrix.nrows(), order + 1);
        assert_eq!(dg_op.stiffness_matrix.ncols(), order + 1);
        assert_eq!(dg_op.diff_matrix.nrows(), order + 1);
        assert_eq!(dg_op.diff_matrix.ncols(), order + 1);
        assert_eq!(dg_op.lift_operators.len(), 2);
    }
    
    #[test]
    fn test_dg_operator_rhs() {
        let order = 2;
        let num_components = 1;
        let params = DGOperatorParams::new()
            .with_volume_flux(FluxType::Central)
            .with_surface_flux(FluxType::LaxFriedrichs)
            .with_limiter(LimiterType::None);
        
        let dg_op = DGOperator::new(order, num_components, Some(params)).unwrap();
        
        // Test with a constant solution u(x) = 1
        let mut u = DMatrix::zeros(1, order + 1);
        u.fill(1.0);
        
        // For a constant solution, the RHS should be zero
        let flux = |u: &DVector<f64>| u.clone();
        let bc = |x: f64, u: &DVector<f64>, _: bool| u.clone();
        
        let rhs = dg_op.rhs(&u, flux, bc).unwrap();
        
        for i in 0..=order {
            assert_relative_eq!(rhs[(0, i)], 0.0, epsilon = 1e-10);
        }
    }
    
    #[test]
    fn test_dg_operator_project() {
        let order = 3;
        let num_components = 1;
        let dg_op = DGOperator::new(order, num_components, None).unwrap();
        
        // Project a quadratic function: f(x) = x² + 2x + 1
        let f = |x: f64| DVector::from_vec(vec![x * x + 2.0 * x + 1.0]);
        
        // Project the function onto the DG basis
        let u = dg_op.project(f);
        
        // Evaluate the projection at the nodes and compare with the exact solution
        let nodes = dg_op.basis().nodes();
        for (i, &xi) in nodes.iter().enumerate() {
            let exact = f(xi)[0];
            let projected = u[(0, i)];
            assert_relative_eq!(
                projected, exact,
                epsilon = 1e-10,
                max_relative = 1e-10,
                "Mismatch at node {} (x = {}): expected {}, got {}",
                i, xi, exact, projected
            );
        }
        
        // Test that the projection is exact for polynomials of degree <= order
        // by checking the L2 error norm
        let num_quad_points = order + 2;  // Use more points for accurate integration
        let quad_points = dg_op.basis().gauss_quadrature(num_quad_points);
        
        let mut l2_error_sq = 0.0;
        let mut l2_norm_sq = 0.0;
        
        for (x, w) in quad_points {
            let exact = f(x)[0];
            let projected = dg_op.evaluate_at(&u, x)[0];
            let diff = exact - projected;
            l2_error_sq += w * diff * diff;
            l2_norm_sq += w * exact * exact;
        }
        
        let l2_error = l2_error_sq.sqrt();
        let l2_norm = l2_norm_sq.sqrt();
        let relative_error = if l2_norm > 0.0 { l2_error / l2_norm } else { l2_error };
        
        assert_relative_eq!(
            relative_error, 0.0,
            epsilon = 1e-12,
            "Projection error too large: relative L2 error = {}",
            relative_error
        );
        
        // Check that the projection is exact for polynomials of degree <= order
        for q in 0..dg_op.num_quad {
            let x_q = dg_op.quad_points[q];
            let u_q = dg_op.phi.row(q).dot(&u.row(0).transpose());
            let f_q = f(x_q)[0];
            
            assert_relative_eq!(u_q, f_q, epsilon = 1e-10);
        }
    }
}
