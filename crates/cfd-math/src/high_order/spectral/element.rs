//! Spectral element data structures and operations.
//!
//! This module provides the core data structures and operations for working with
//! spectral elements, including node definitions, basis functions, and element-local
//! operations.

use super::{Result, SpectralError, compute_lgl_nodes, compute_lgl_weights, compute_derivative_matrix};
use nalgebra::{DMatrix, DVector};

/// Represents a spectral element with nodes, weights, and differentiation matrices
#[derive(Debug, Clone)]
pub struct SpectralElement {
    /// Polynomial order
    pub order: usize,
    /// Number of nodes (order + 1)
    pub num_nodes: usize,
    /// Legendre-Gauss-Lobatto nodes
    pub nodes: DVector<f64>,
    /// Quadrature weights
    pub weights: DVector<f64>,
    /// First derivative matrix
    pub derivative_matrix: DMatrix<f64>,
    /// Mass matrix (diagonal)
    pub mass_matrix: DMatrix<f64>,
    /// Stiffness matrix
    pub stiffness_matrix: DMatrix<f64>,
}

impl SpectralElement {
    /// Create a new spectral element of given order
    pub fn new(order: usize) -> Result<Self> {
        if order < 1 {
            return Err(SpectralError::InvalidOrder(order));
        }

        let num_nodes = order + 1;
        let nodes = DVector::from_vec(compute_lgl_nodes(order)?);
        let weights = DVector::from_vec(compute_lgl_weights(nodes.as_slice(), order));
        
        // Compute derivative matrix
        let derivative_matrix = compute_derivative_matrix(nodes.as_slice(), order);
        
        // Compute mass matrix (diagonal)
        let mut mass_matrix = DMatrix::zeros(num_nodes, num_nodes);
        for i in 0..num_nodes {
            mass_matrix[(i, i)] = weights[i];
        }
        
        // Compute stiffness matrix K_ij = ∫ φ'_i φ'_j dx
        let stiffness_matrix = &derivative_matrix.transpose() * &mass_matrix * &derivative_matrix;
        
        Ok(Self {
            order,
            num_nodes,
            nodes,
            weights,
            derivative_matrix,
            mass_matrix,
            stiffness_matrix,
        })
    }
    
    /// Evaluate the i-th Lagrange basis function at point x
    pub fn lagrange_basis(&self, i: usize, x: f64) -> f64 {
        let mut result = 1.0;
        let xi = self.nodes[i];
        
        for (j, &xj) in self.nodes.iter().enumerate() {
            if i != j {
                result *= (x - xj) / (xi - xj);
            }
        }
        
        result
    }
    
    /// Evaluate the derivative of the i-th Lagrange basis function at point x
    pub fn d_lagrange_basis(&self, i: usize, x: f64) -> f64 {
        let mut result = 0.0;
        let xi = self.nodes[i];
        
        for (j, &xj) in self.nodes.iter().enumerate() {
            if i == j {
                continue;
            }
            
            let mut term = 1.0 / (xi - xj);
            
            for (k, &xk) in self.nodes.iter().enumerate() {
                if k != i && k != j {
                    term *= (x - xk) / (xi - xk);
                }
            }
            
            result += term;
        }
        
        result
    }
    
    /// Interpolate a function at the element nodes
    pub fn interpolate<F>(&self, f: F) -> DVector<f64>
    where
        F: Fn(f64) -> f64,
    {
        self.nodes.map(f)
    }
    
    /// Compute the L2 error between a function and its interpolant
    pub fn l2_error<F>(&self, f: F, u: &DVector<f64>) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let mut error = 0.0;
        
        for (i, &xi) in self.nodes.iter().enumerate() {
            let diff = f(xi) - u[i];
            error += self.weights[i] * diff * diff;
        }
        
        error.sqrt()
    }
    
    /// Compute the derivative of a function represented by its nodal values
    pub fn derivative(&self, u: &DVector<f64>) -> DVector<f64> {
        &self.derivative_matrix * u
    }
    
    /// Compute the integral of a function represented by its nodal values
    pub fn integrate(&self, u: &DVector<f64>) -> f64 {
        self.weights.dot(u)
    }
}

/// A 1D spectral element mesh
#[derive(Debug, Clone)]
pub struct SpectralMesh1D {
    elements: Vec<SpectralElement>,
    num_elements: usize,
    #[allow(dead_code)]
    element_order: usize,
    global_nodes: Vec<f64>,
    element_connectivity: Vec<Vec<usize>>,
}

impl SpectralMesh1D {
    /// Create a new 1D spectral element mesh
    pub fn new(x_min: f64, x_max: f64, num_elements: usize, element_order: usize) -> Result<Self> {
        if num_elements < 1 {
            return Err(SpectralError::InvalidOrder(num_elements));
        }
        
        let element = SpectralElement::new(element_order)?;
        let num_nodes_per_element = element.num_nodes;
        
        // Create elements
        let mut elements = Vec::with_capacity(num_elements);
        let mut global_nodes = Vec::with_capacity((num_elements * (element_order + 1)) - (num_elements - 1));
        let mut element_connectivity = Vec::with_capacity(num_elements);
        
        // Generate global nodes and connectivity
        let dx = (x_max - x_min) / num_elements as f64;
        
        for e in 0..num_elements {
            let x0 = x_min + e as f64 * dx;
            let x1 = x0 + dx;
            
            // Map local nodes to global coordinates
            let local_nodes: Vec<f64> = element.nodes.iter()
                .map(|&xi| x0 + 0.5 * (xi + 1.0) * (x1 - x0))
                .collect();
            
            // Update global nodes and connectivity
            let mut connectivity = Vec::with_capacity(num_nodes_per_element);
            
            for (local_idx, &x) in local_nodes.iter().enumerate() {
                if e == 0 || local_idx > 0 {  // Skip duplicate nodes
                    connectivity.push(global_nodes.len());
                    global_nodes.push(x);
                } else {
                    // Connect to previous element's last node
                    connectivity.push(global_nodes.len() - 1);
                }
            }
            
            element_connectivity.push(connectivity);
            elements.push(element.clone());
        }
        
        Ok(Self {
            elements,
            num_elements,
            element_order,
            global_nodes,
            element_connectivity,
        })
    }
    
    /// Get the number of global nodes
    pub fn num_global_nodes(&self) -> usize {
        self.global_nodes.len()
    }
    
    /// Get the number of elements
    pub fn num_elements(&self) -> usize {
        self.num_elements
    }
    
    /// Get the global coordinates of a node
    pub fn node_coord(&self, node_idx: usize) -> f64 {
        self.global_nodes[node_idx]
    }
    
    /// Get the element containing a point
    pub fn find_element(&self, x: f64) -> Option<usize> {
        for (e, _element) in self.elements.iter().enumerate() {
            let x0 = self.node_coord(self.element_connectivity[e][0]);
            let x1 = self.node_coord(*self.element_connectivity[e].last().unwrap());
            
            // Use [x0, x1) for all but the last element, which is [x0, x1]
            let is_last = e == self.num_elements - 1;
            if x >= x0 && (x < x1 || (is_last && x <= x1 + 1e-14)) {
                return Some(e);
            }
        }
        
        None
    }
    
    /// Interpolate a function to the global nodes
    pub fn interpolate_global<F>(&self, f: F) -> DVector<f64>
    where
        F: Fn(f64) -> f64,
    {
        let mut u = DVector::zeros(self.num_global_nodes());
        
        for e in 0..self.num_elements() {
            let conn = &self.element_connectivity[e];
            
            for &global in conn {
                let x = self.node_coord(global);
                u[global] = f(x);
            }
        }
        
        u
    }
    
    /// Compute the L2 error between a function and its global interpolant
    pub fn l2_error_global<F>(&self, f: F, u: &DVector<f64>) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let mut error = 0.0;
        
        for e in 0..self.num_elements() {
            let conn = &self.element_connectivity[e];
            let element = &self.elements[e];
            
            // Extract local solution
            let mut u_local = DVector::zeros(element.num_nodes);
            for (local, &global) in conn.iter().enumerate() {
                u_local[local] = u[global];
            }
            
            // Compute element contribution to error
            let x0 = self.node_coord(conn[0]);
            let x1 = self.node_coord(*conn.last().unwrap());
            
            // Use Gauss quadrature for more accurate error computation
            let quad_points = 10;  // Number of quadrature points per element
            let dx = (x1 - x0) / f64::from(quad_points);
            
            for i in 0..=quad_points {
                let x = x0 + f64::from(i) * dx;
                let xi = 2.0 * (x - x0) / (x1 - x0) - 1.0;  // Map to reference element [-1,1]
                
                // Evaluate interpolant at quadrature point
                let mut uh = 0.0;
                for (j, &_node) in element.nodes.iter().enumerate() {
                    uh += u_local[j] * element.lagrange_basis(j, xi);
                }
                
                // Add contribution to error
                let weight = if i == 0 || i == quad_points { 0.5 } else { 1.0 };
                error += weight * dx * (f(x) - uh).powi(2);
            }
        }
        
        error.sqrt()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_spectral_element() {
        let order = 4;
        let element = SpectralElement::new(order).unwrap();
        
        // Test that the element has the correct number of nodes
        assert_eq!(element.num_nodes, order + 1);
        
        // Test that the nodes are in [-1, 1] and symmetric
        assert_relative_eq!(element.nodes[0], -1.0, epsilon = 1e-10);
        assert_relative_eq!(element.nodes[order], 1.0, epsilon = 1e-10);
        
        for i in 0..=order/2 {
            assert_relative_eq!(
                element.nodes[i], 
                -element.nodes[order - i],
                epsilon = 1e-10
            );
        }
        
        // Test that the mass matrix is diagonal
        for i in 0..element.num_nodes {
            for j in 0..element.num_nodes {
                if i != j {
                    assert_relative_eq!(element.mass_matrix[(i, j)], 0.0, epsilon = 1e-10);
                }
            }
        }
        
        // Test that the derivative of a linear function is exact
        let a = 2.0;
        let b = 3.0;
        let f = |x: f64| a * x + b;
        let df = |_: f64| a;
        
        let u = element.interpolate(f);
        let du = element.derivative(&u);
        
        for i in 0..element.num_nodes {
            assert_relative_eq!(du[i], df(element.nodes[i]), epsilon = 1e-10);
        }
    }
    
    #[test]
    fn test_spectral_mesh() {
        let x_min = 0.0;
        let x_max = 1.0;
        let num_elements = 4;
        let element_order = 3;
        
        let mesh = SpectralMesh1D::new(x_min, x_max, num_elements, element_order).unwrap();
        
        // Test number of global nodes
        let expected_nodes = num_elements * element_order + 1;
        assert_eq!(mesh.num_global_nodes(), expected_nodes);
        
        // Test node coordinates
        assert_relative_eq!(mesh.node_coord(0), x_min, epsilon = 1e-10);
        assert_relative_eq!(mesh.node_coord(mesh.num_global_nodes() - 1), x_max, epsilon = 1e-10);
        
        // Test element finding
        assert_eq!(mesh.find_element(0.0), Some(0));
        assert_eq!(mesh.find_element(0.25), Some(1));
        assert_eq!(mesh.find_element(0.5), Some(2));
        assert_eq!(mesh.find_element(0.75), Some(3));
        assert_eq!(mesh.find_element(1.0), Some(3));
        
        // Test interpolation and error computation
        let f = |x: f64| x.powi(2);
        let u = mesh.interpolate_global(f);
        
        // The error should be small for a quadratic function (exact for order >= 2)
        let error = mesh.l2_error_global(f, &u);
        assert!(error < 1e-10);
        
        // Test with a higher-order function
        let f = |x: f64| x.powi(4);
        let u = mesh.interpolate_global(f);
        let error = mesh.l2_error_global(f, &u);
        
        // The error should be small but non-zero for a quartic function with cubic elements
        assert!(error < 1e-4);
    }
}
