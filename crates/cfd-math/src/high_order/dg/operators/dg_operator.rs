use super::super::basis::{BasisType, DGBasis};
use super::super::FluxType;
use super::super::{
    column_vector, matrix_solve, matrix_sub, matrix_zeros, row_vector, set_column, set_row,
    vector_add, vector_add_assign_scaled, vector_len, vector_scale, vector_sub, vector_sum,
    vector_zeros,
};
use super::params::*;
use crate::error::Result;
use cfd_core::error::{Error, ErrorContext};
use leto::{Array1, Array2};

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
    pub lift_operators: Vec<Array2<f64>>,
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
            return Err(Error::InvalidInput(format!(
                "Polynomial order must be at least 1, got {order}"
            )));
        }

        let params = params.unwrap_or_default();
        let basis = DGBasis::new(order, basis_type).context("constructing DG operator basis")?;
        let num_basis = order + 1;

        // Compute lifting operators for boundary conditions
        let mut lift_operators = Vec::with_capacity(2);

        // Left boundary lifting operator
        let mut lift_l = matrix_zeros(num_basis, 1);
        for i in 0..num_basis {
            lift_l[[i, 0]] = basis.phi[[i, 0]] * basis.quad_weights[0];
        }
        lift_operators.push(lift_l);

        // Right boundary lifting operator
        let mut lift_r = matrix_zeros(num_basis, 1);
        for i in 0..num_basis {
            lift_r[[i, 0]] = basis.phi[[i, vector_len(&basis.quad_points) - 1]]
                * basis.quad_weights[vector_len(&basis.quad_points) - 1];
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
    pub fn compute_derivative(&self, u: &Array2<f64>) -> Result<Array2<f64>> {
        let num_basis = self.order + 1;
        let mut du_dx = matrix_zeros(self.num_components, num_basis);

        // M * c_deriv = K^T * c_orig
        // where K_ij = ∫ φ_i' φ_j dx
        // So (K^T)_ij = K_ji = ∫ φ_j' φ_i dx

        for i in 0..self.num_components {
            let u_row = row_vector(u, i);
            let rhs =
                super::super::matrix_transpose_vector_mul(&self.basis.stiffness_matrix, &u_row);
            let sol = matrix_solve(&self.basis.mass_matrix, &rhs)
                .context("solving DG mass matrix for derivative coefficients")?;
            set_row(&mut du_dx, i, &sol);
        }

        Ok(du_dx)
    }

    /// Compute the right-hand side of the DG semi-discrete equation
    pub fn rhs<F, G>(&self, u: &Array2<f64>, flux: F, boundary_condition: G) -> Result<Array2<f64>>
    where
        F: Fn(&Array1<f64>) -> Array1<f64>,
        G: Fn(f64, &Array1<f64>, bool) -> Array1<f64>,
    {
        let num_basis = self.order + 1;
        let num_quad = vector_len(&self.basis.quad_points);
        let mut rhs = matrix_zeros(self.num_components, num_basis);

        // Compute volume integral terms
        if self.params.strong_form {
            // Strong form: M^{-1} * (K^T * F - F^*|_{-1}^1)

            // Compute the flux at all quadrature points
            let mut f_quad = matrix_zeros(self.num_components, num_quad);
            for q in 0..num_quad {
                // Evaluate the solution at the quadrature point
                let mut u_q = vector_zeros(self.num_components);
                for i in 0..num_basis {
                    vector_add_assign_scaled(
                        &mut u_q,
                        &column_vector(u, i),
                        self.basis.phi[[i, q]],
                    );
                }

                // Compute the flux
                set_column(&mut f_quad, q, &flux(&u_q));
            }

            // Compute K^T * F
            let mut kt_f = matrix_zeros(self.num_components, num_basis);
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    let mut sum = 0.0;
                    for q in 0..num_quad {
                        sum += self.basis.quad_weights[q]
                            * self.basis.dphi_dx[[j, q]]
                            * f_quad[[i, q]];
                    }
                    kt_f[[i, j]] = sum;
                }
            }

            // Compute the numerical flux at the boundaries
            // Evaluate solution at boundaries correctly for any basis
            let mut u_l = vector_zeros(self.num_components);
            let mut u_r = vector_zeros(self.num_components);
            for i in 0..num_basis {
                vector_add_assign_scaled(
                    &mut u_l,
                    &column_vector(u, i),
                    self.basis.evaluate_basis(i, -1.0),
                );
                vector_add_assign_scaled(
                    &mut u_r,
                    &column_vector(u, i),
                    self.basis.evaluate_basis(i, 1.0),
                );
            }

            // Apply boundary conditions
            let f_l = boundary_condition(-1.0, &u_l, false);
            let f_r = boundary_condition(1.0, &u_r, true);

            // Compute the numerical flux
            let f_num_l = self.compute_numerical_flux(-1.0, &f_l, &u_l, None, &boundary_condition);
            let f_num_r = self.compute_numerical_flux(1.0, &f_r, &u_r, None, &boundary_condition);

            // Compute the boundary terms
            let mut boundary_term = matrix_zeros(self.num_components, num_basis);
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    boundary_term[[i, j]] = self.basis.evaluate_basis(j, 1.0) * f_num_r[i]
                        - self.basis.evaluate_basis(j, -1.0) * f_num_l[i];
                }
            }

            // Compute the right-hand side
            rhs = matrix_sub(&kt_f, &boundary_term);

            // Multiply by the inverse mass matrix
            for i in 0..self.num_components {
                let rhs_row = row_vector(&rhs, i);
                let sol = matrix_solve(&self.basis.mass_matrix, &rhs_row)
                    .context("solving DG mass matrix for strong-form RHS")?;
                set_row(&mut rhs, i, &sol);
            }
        } else {
            // Weak form: F * dφ/dx - F^* φ|_{-1}^1

            // Compute the flux at all quadrature points
            let mut f_quad = matrix_zeros(self.num_components, num_quad);
            for q in 0..num_quad {
                // Evaluate the solution at the quadrature point
                let mut u_q = vector_zeros(self.num_components);
                for i in 0..num_basis {
                    vector_add_assign_scaled(
                        &mut u_q,
                        &column_vector(u, i),
                        self.basis.phi[[i, q]],
                    );
                }

                // Compute the flux
                set_column(&mut f_quad, q, &flux(&u_q));
            }

            // Compute the volume integral: F * dφ/dx
            for i in 0..self.num_components {
                for j in 0..num_basis {
                    let mut sum = 0.0;
                    for q in 0..num_quad {
                        sum += self.basis.quad_weights[q]
                            * f_quad[[i, q]]
                            * self.basis.dphi_dx[[j, q]];
                    }
                    rhs[[i, j]] = -sum; // Negative sign from integration by parts
                }
            }

            // Compute the numerical flux at the boundaries
            // Evaluate solution at boundaries correctly for any basis
            let mut u_l = vector_zeros(self.num_components);
            let mut u_r = vector_zeros(self.num_components);
            for i in 0..num_basis {
                vector_add_assign_scaled(
                    &mut u_l,
                    &column_vector(u, i),
                    self.basis.evaluate_basis(i, -1.0),
                );
                vector_add_assign_scaled(
                    &mut u_r,
                    &column_vector(u, i),
                    self.basis.evaluate_basis(i, 1.0),
                );
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
                    rhs[[i, j]] += self.basis.evaluate_basis(j, 1.0) * f_num_r[i]
                        - self.basis.evaluate_basis(j, -1.0) * f_num_l[i];
                }
            }

            // Multiply by the inverse mass matrix
            for i in 0..self.num_components {
                let rhs_row = row_vector(&rhs, i);
                let sol = matrix_solve(&self.basis.mass_matrix, &rhs_row)
                    .context("solving DG mass matrix for weak-form RHS")?;
                set_row(&mut rhs, i, &sol);
            }
        }

        Ok(rhs)
    }

    /// Compute the numerical flux at a boundary
    fn compute_numerical_flux<G>(
        &self,
        x: f64,
        f: &Array1<f64>,
        u: &Array1<f64>,
        u_ext: Option<&Array1<f64>>,
        boundary_condition: G,
    ) -> Array1<f64>
    where
        G: Fn(f64, &Array1<f64>, bool) -> Array1<f64>,
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
            FluxType::Central => vector_scale(&vector_add(f, &f_ext), 0.5),
            FluxType::Upwind => {
                // Upwind direction determined by the characteristic speed.
                // Use the local normal velocity: average of interior and exterior
                // state projections, falling back to params.alpha for sign.
                let a = 0.5 * (vector_sum(u) + vector_sum(&u_ext));
                let a = if a.abs() < f64::EPSILON {
                    self.params.alpha
                } else {
                    a
                };
                if a >= 0.0 {
                    f.clone()
                } else {
                    f_ext
                }
            }
            FluxType::LaxFriedrichs => {
                let alpha = self.params.alpha;
                vector_sub(
                    &vector_scale(&vector_add(f, &f_ext), 0.5),
                    &vector_scale(&vector_sub(&u_ext, u), 0.5 * alpha),
                )
            }
            _ => {
                // Default to Lax-Friedrichs
                let alpha = self.params.alpha;
                vector_sub(
                    &vector_scale(&vector_add(f, &f_ext), 0.5),
                    &vector_scale(&vector_sub(&u_ext, u), 0.5 * alpha),
                )
            }
        }
    }

    /// Apply boundary conditions to the solution
    pub fn apply_boundary_conditions<G>(&self, u: &mut Array2<f64>, boundary_condition: G)
    where
        G: Fn(f64, &Array1<f64>, bool) -> Array1<f64>,
    {
        let num_basis = self.order + 1;

        // Apply boundary conditions at the left boundary (x = -1)
        let u_l = column_vector(u, 0);
        let u_l_bc = boundary_condition(-1.0, &u_l, false);

        // Apply boundary conditions at the right boundary (x = 1)
        let u_r = column_vector(u, num_basis - 1);
        let u_r_bc = boundary_condition(1.0, &u_r, true);

        // Update the solution with the boundary values
        for i in 0..self.num_components {
            u[[i, 0]] = u_l_bc[i];
            u[[i, num_basis - 1]] = u_r_bc[i];
        }
    }

    /// Project a function onto the DG basis
    pub fn project<F>(&self, f: F) -> Result<Array2<f64>>
    where
        F: Fn(f64) -> Array1<f64>,
    {
        let num_basis = self.order + 1;
        let num_quad = vector_len(&self.basis.quad_points);
        let mut u = matrix_zeros(self.num_components, num_basis);

        // Compute the projection
        for q in 0..num_quad {
            let x_q = self.basis.quad_points[q];
            let w_q = self.basis.quad_weights[q];
            let f_q = f(x_q);

            for i in 0..num_basis {
                let phi_i = self.basis.phi[[i, q]];
                for c in 0..self.num_components {
                    u[[c, i]] += w_q * f_q[c] * phi_i;
                }
            }
        }

        // Solve the mass matrix system
        for i in 0..self.num_components {
            let u_row = row_vector(&u, i);
            let sol = matrix_solve(&self.basis.mass_matrix, &u_row)
                .context("solving DG mass matrix for projection")?;
            set_row(&mut u, i, &sol);
        }

        Ok(u)
    }
}
