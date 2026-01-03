//! Core operators for Discontinuous Galerkin methods.
//!
//! This module implements the spatial discretization operators for DG methods,
//! including volume and surface integrals, weak form operators, and boundary conditions.

use super::basis::{BasisType, DGBasis};
use super::{FluxType, LimiterType, LimiterParams, FluxParams, Result, DGError};
use nalgebra::{DVector, DMatrix};

/// Type of boundary condition
#[derive(Debug, Clone, Copy, PartialEq)]
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
#[derive(Clone)]
pub struct DGOperator {
    /// Polynomial order
    pub order: usize,
    /// Number of components (for systems of equations)
    pub num_components: usize,
    /// Basis functions and matrices
    pub basis: DGBasis,
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
        Self::with_basis(order, num_components, BasisType::Orthogonal, params)
    }

    /// Create a new DG operator with given basis type
    pub fn with_basis(
        order: usize,
        num_components: usize,
        basis_type: BasisType,
        params: Option<DGOperatorParams>,
    ) -> Result<Self> {
        if order == 0 {
            return Err(DGError::InvalidOrder(order));
        }
        
        let params = params.unwrap_or_default();
        let basis = DGBasis::new(order, basis_type)?;
        let num_basis = order + 1;
        
        // Compute lifting operators for boundary conditions
        let mut lift_operators = Vec::with_capacity(2);
        
        // Left boundary lifting operator
        let mut lift_l = DMatrix::zeros(num_basis, 1);
        for i in 0..num_basis {
            lift_l[(i, 0)] = basis.phi[(i, 0)] * basis.quad_weights[0];
        }
        lift_operators.push(lift_l);
        
        // Right boundary lifting operator
        let mut lift_r = DMatrix::zeros(num_basis, 1);
        for i in 0..num_basis {
            lift_r[(i, 0)] = basis.phi[(i, basis.quad_points.len() - 1)] * basis.quad_weights[basis.quad_points.len() - 1];
        }
        lift_operators.push(lift_r);
        
        Ok(Self {
            order,
            num_components,
            basis,
            lift_operators,
            params,
        })
    }
    
    /// Compute the derivative of the solution
    pub fn compute_derivative(&self, u: &DMatrix<f64>) -> Result<DMatrix<f64>> {
        let num_basis = self.order + 1;
        let mut du_dx = DMatrix::zeros(self.num_components, num_basis);
        
        // M * c_deriv = K^T * c_orig
        // where K_ij = ∫ φ_i' φ_j dx
        // So (K^T)_ij = K_ji = ∫ φ_j' φ_i dx
        
        let kt = self.basis.stiffness_matrix.transpose();
        let mass_lu = self.basis.mass_matrix.clone().lu();
        
        for i in 0..self.num_components {
            let u_row = u.row(i).transpose();
            let rhs = &kt * u_row;
            if let Some(sol) = mass_lu.solve(&rhs) {
                du_dx.row_mut(i).copy_from(&sol.transpose());
            } else {
                return Err(DGError::NumericalError("Failed to solve mass matrix system for derivative".to_string()));
            }
        }
        
        Ok(du_dx)
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
        let num_quad = self.basis.quad_points.len();
        let mut rhs = DMatrix::zeros(self.num_components, num_basis);
        
        // Compute volume integral terms
        if self.params.strong_form {
            // Strong form: M^{-1} * (K^T * F - F^*|_{-1}^1)
            
            // Compute the flux at all quadrature points
            let mut f_quad = DMatrix::zeros(self.num_components, num_quad);
            for q in 0..num_quad {
                // Evaluate the solution at the quadrature point
                let mut u_q = DVector::zeros(self.num_components);
                for i in 0..num_basis {
                    u_q += u.column(i) * self.basis.phi[(i, q)];
                }
                
                // Compute the flux
                f_quad.column_mut(q).copy_from(&flux(&u_q));
            }
            
            // Compute K^T * F
            let mut kt_f = DMatrix::zeros(self.num_components, num_basis);
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    let mut sum = 0.0;
                    for q in 0..num_quad {
                        sum += self.basis.quad_weights[q] * self.basis.dphi_dx[(j, q)] * f_quad[(i, q)];
                    }
                    kt_f[(i, j)] = sum;
                }
            }
            
            // Compute the numerical flux at the boundaries
            // Evaluate solution at boundaries correctly for any basis
            let mut u_l = DVector::zeros(self.num_components);
            let mut u_r = DVector::zeros(self.num_components);
            for i in 0..num_basis {
                u_l += u.column(i) * self.basis.evaluate_basis(i, -1.0);
                u_r += u.column(i) * self.basis.evaluate_basis(i, 1.0);
            }
            
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
                        self.basis.evaluate_basis(j, 1.0) * f_num_r[i] - 
                        self.basis.evaluate_basis(j, -1.0) * f_num_l[i];
                }
            }
            
            // Compute the right-hand side
            rhs = kt_f - boundary_term;
            
            // Multiply by the inverse mass matrix
            let mass_lu = self.basis.mass_matrix.clone().lu();
            for i in 0..self.num_components {
                let rhs_row = rhs.row(i).transpose();
                if let Some(sol) = mass_lu.solve(&rhs_row) {
                    rhs.row_mut(i).copy_from(&sol.transpose());
                }
            }
        } else {
            // Weak form: F * dφ/dx - F^* φ|_{-1}^1
            
            // Compute the flux at all quadrature points
            let mut f_quad = DMatrix::zeros(self.num_components, num_quad);
            for q in 0..num_quad {
                // Evaluate the solution at the quadrature point
                let mut u_q = DVector::zeros(self.num_components);
                for i in 0..num_basis {
                    u_q += u.column(i) * self.basis.phi[(i, q)];
                }
                
                // Compute the flux
                f_quad.column_mut(q).copy_from(&flux(&u_q));
            }
            
            // Compute the volume integral: F * dφ/dx
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    let mut sum = 0.0;
                    for q in 0..num_quad {
                        sum += self.basis.quad_weights[q] * f_quad[(i, q)] * self.basis.dphi_dx[(j, q)];
                    }
                    rhs[(i, j)] = -sum; // Negative sign from integration by parts
                }
            }
            
            // Compute the numerical flux at the boundaries
            // Evaluate solution at boundaries correctly for any basis
            let mut u_l = DVector::zeros(self.num_components);
            let mut u_r = DVector::zeros(self.num_components);
            for i in 0..num_basis {
                u_l += u.column(i) * self.basis.evaluate_basis(i, -1.0);
                u_r += u.column(i) * self.basis.evaluate_basis(i, 1.0);
            }
            
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
                        self.basis.evaluate_basis(j, 1.0) * f_num_r[i] - 
                        self.basis.evaluate_basis(j, -1.0) * f_num_l[i];
                }
            }
            
            // Multiply by the inverse mass matrix
            let mass_lu = self.basis.mass_matrix.clone().lu();
            for i in 0..self.num_components {
                let rhs_row = rhs.row(i).transpose();
                if let Some(sol) = mass_lu.solve(&rhs_row) {
                    rhs.row_mut(i).copy_from(&sol.transpose());
                }
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
        let num_quad = self.basis.quad_points.len();
        let mut u = DMatrix::zeros(self.num_components, num_basis);
        
        // Compute the projection
        for q in 0..num_quad {
            let x_q = self.basis.quad_points[q];
            let w_q = self.basis.quad_weights[q];
            let f_q = f(x_q);
            
            for i in 0..num_basis {
                let phi_i = self.basis.phi[(i, q)];
                for c in 0..self.num_components {
                    u[(c, i)] += w_q * f_q[c] * phi_i;
                }
            }
        }
        
        // Solve the mass matrix system
        let mass_lu = self.basis.mass_matrix.clone().lu();
        for i in 0..self.num_components {
            let u_row = u.row(i).transpose();
            if let Some(sol) = mass_lu.solve(&u_row) {
                u.row_mut(i).copy_from(&sol.transpose());
            }
        }
        
        u
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_dg_operator_new() {
        let order = 2;
        let num_components = 1;
        let dg_op = DGOperator::new(order, num_components, None).unwrap();
        
        assert_eq!(dg_op.order, order);
        assert_eq!(dg_op.num_components, num_components);
        // Orthogonal basis uses order + 2 quadrature points
        assert_eq!(dg_op.basis.quad_points.len(), order + 2);
        assert_eq!(dg_op.basis.quad_weights.len(), order + 2);
        assert_eq!(dg_op.basis.phi.nrows(), order + 1);
        assert_eq!(dg_op.basis.phi.ncols(), order + 2);
        assert_eq!(dg_op.basis.mass_matrix.nrows(), order + 1);
        assert_eq!(dg_op.basis.mass_matrix.ncols(), order + 1);
        assert_eq!(dg_op.basis.stiffness_matrix.nrows(), order + 1);
        assert_eq!(dg_op.basis.stiffness_matrix.ncols(), order + 1);
        assert_eq!(dg_op.basis.diff_matrix.nrows(), order + 1);
        assert_eq!(dg_op.basis.diff_matrix.ncols(), order + 1);
    }

    #[test]
    fn test_dg_operator_compute_derivative() {
        let order = 2;
        let num_components = 1;
        let dg_op = DGOperator::new(order, num_components, None).unwrap();
        
        // u(x) = x -> u' (x) = 1
        // In Legendre basis: u(x) = P_1(x) -> u' (x) = P_0(x)
        let mut u = DMatrix::zeros(1, order + 1);
        u[(0, 1)] = 1.0;
        
        let du_dx = dg_op.compute_derivative(&u).unwrap();
        
        assert_relative_eq!(du_dx[(0, 0)], 1.0, epsilon = 1e-10);
        assert_relative_eq!(du_dx[(0, 1)], 0.0, epsilon = 1e-10);
        assert_relative_eq!(du_dx[(0, 2)], 0.0, epsilon = 1e-10);
    }

    #[test]
    fn test_dg_operator_rhs() {
        let order = 2;
        let num_components = 1;
        let dg_op = DGOperator::new(order, num_components, None).unwrap();
        
        // Test with a constant solution u(x) = 1
        let mut u = DMatrix::zeros(1, order + 1);
        u[(0, 0)] = 1.0; // P_0(x) = 1
        
        // For a constant solution, the RHS should be zero
        let flux = |u: &DVector<f64>| u.clone();
        let bc = |_: f64, u: &DVector<f64>, _: bool| u.clone();
        
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
        
        // Project f(x) = x^2
        let f = |x: f64| DVector::from_element(1, x * x);
        let u = dg_op.project(f);
        
        // Verify projection at quadrature points
        for (i, &xi) in dg_op.basis.quad_points.iter().enumerate() {
            let mut projected = 0.0;
            for j in 0..=order {
                projected += u[(0, j)] * dg_op.basis.phi[(j, i)];
            }
            assert_relative_eq!(projected, xi * xi, epsilon = 1e-10);
        }
        
        // Verify L2 error
        let mut l2_error_sq = 0.0;
        for (q, &x) in dg_op.basis.quad_points.iter().enumerate() {
            let w = dg_op.basis.quad_weights[q];
            let mut projected = 0.0;
            for j in 0..=order {
                projected += u[(0, j)] * dg_op.basis.phi[(j, q)];
            }
            let error = projected - x * x;
            l2_error_sq += w * error * error;
        }
        assert!(l2_error_sq.sqrt() < 1e-10);
    }

    #[test]
    fn test_dg_operator_numerical_flux() {
        let order = 2;
        let num_components = 1;
        let dg_op = DGOperator::new(order, num_components, None).unwrap();
        
        let u = DVector::from_element(1, 1.0);
        let f = DVector::from_element(1, 1.0);
        let bc = |_: f64, u: &DVector<f64>, _: bool| u.clone();
        
        // Central flux: f_num = 0.5 * (f_l + f_r)
        let f_num = dg_op.compute_numerical_flux(0.0, &f, &u, Some(&u), &bc);
        assert_relative_eq!(f_num[0], 1.0, epsilon = 1e-10);
    }
}
