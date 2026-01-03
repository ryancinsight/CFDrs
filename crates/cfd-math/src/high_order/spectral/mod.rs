//! Spectral element methods for high-order numerical solutions.
//!
//! This module implements spectral element methods (SEM) for solving partial differential
//! equations with exponential convergence rates for smooth solutions.
//!
//! # Features
//! - Legendre-Gauss-Lobatto (LGL) nodes and weights
//! - Lagrange basis functions and their derivatives
//! - Spectral differentiation matrices
//! - Elemental and global assembly

mod element;
mod operators;
mod assembly;

pub use element::*;
pub use operators::*;
pub use assembly::*;

use nalgebra::DMatrix;
use std::f64::consts::PI;

/// Error types for spectral element methods
#[derive(Debug, thiserror::Error)]
pub enum SpectralError {
    /// Invalid polynomial order
    #[error("Polynomial order must be at least 1, got {0}")]
    InvalidOrder(usize),
    
    /// Node computation failed to converge
    #[error("Failed to compute nodes: {0}")]
    NodeComputation(String),
    
    /// Matrix operation error
    #[error("Matrix operation failed: {0}")]
    MatrixError(String),
}

/// Result type for spectral element operations
pub type Result<T> = std::result::Result<T, SpectralError>;

/// Legendre polynomial evaluation
fn legendre_poly(n: usize, x: f64) -> f64 {
    match n {
        0 => 1.0,
        1 => x,
        _ => {
            let mut p0 = 1.0;
            let mut p1 = x;
            
            for k in 2..=n {
                let p2 = ((2 * k - 1) as f64 * x * p1 - (k - 1) as f64 * p0) / k as f64;
                p0 = p1;
                p1 = p2;
            }
            
            p1
        }
    }
}

/// Derivative of Legendre polynomial
#[allow(dead_code)]
fn legendre_poly_deriv(n: usize, x: f64) -> f64 {
    if n == 0 {
        return 0.0;
    }
    
    // Handle endpoints to avoid division by zero
    if (x - 1.0).abs() < 1e-15 {
        return (n * (n + 1)) as f64 / 2.0;
    }
    if (x + 1.0).abs() < 1e-15 {
        let val = (n * (n + 1)) as f64 / 2.0;
        return if (n - 1) % 2 == 0 { val } else { -val };
    }
    
    (n as f64 / (1.0 - x * x)) * (legendre_poly(n - 1, x) - x * legendre_poly(n, x))
}

/// Compute Legendre-Gauss-Lobatto nodes using Newton's method
fn compute_lgl_nodes(n: usize) -> Result<Vec<f64>> {
    if n < 1 {
        return Err(SpectralError::InvalidOrder(n));
    }
    
    let mut nodes = vec![0.0; n + 1];
    nodes[0] = -1.0;
    nodes[n] = 1.0;
    
    // Initial guess for interior nodes (Chebyshev points)
    for i in 1..n {
        nodes[i] = -((i as f64 * PI) / n as f64).cos();
    }
    
    // Newton iteration to find roots of P_n'
    let max_iter = 100;
    let tol = 1e-15;
    
    for i in 1..n {
        let mut x = nodes[i];
        let mut iter = 0;
        let mut delta = 1.0;
        
        while delta > tol && iter < max_iter {
            let (p, dp) = legendre_poly_deriv_with_prev(n, x);
            
            // P''_n = (2x P'_n - n(n+1) P_n) / (1-x^2)
            let denominator = 1.0 - x * x;
            if denominator.abs() < 1e-15 {
                 return Err(SpectralError::NodeComputation(format!("Newton iteration hit endpoint at x={}", x)));
            }
            let d2p = (2.0 * x * dp - (n * (n + 1)) as f64 * p) / denominator;

            if d2p.abs() < f64::EPSILON {
                return Err(SpectralError::NodeComputation(format!(
                    "Zero second derivative at x = {x} for n = {n}"
                )));
            }
            
            let dx = dp / d2p;
            x -= dx;
            delta = dx.abs();
            iter += 1;
        }
        
        if iter >= max_iter {
            return Err(SpectralError::NodeComputation(format!(
                "Failed to converge for node {i} of {n}"
            )));
        }
        
        nodes[i] = x;
    }
    
    // Ensure symmetry
    for i in 0..=(n / 2) {
        let val = nodes[i].abs();
        nodes[i] = -val;
        nodes[n - i] = val;
    }
    
    Ok(nodes)
}

/// Compute Legendre polynomial and its derivative using recurrence
fn legendre_poly_deriv_with_prev(n: usize, x: f64) -> (f64, f64) {
    if n == 0 {
        return (1.0, 0.0);
    }
    
    let mut p_prev = 1.0;
    let mut p_curr = x;
    let mut dp_prev = 0.0;
    let mut dp_curr = 1.0;
    
    for k in 1..n {
        let p_next = ((2 * k + 1) as f64 * x * p_curr - k as f64 * p_prev) / (k + 1) as f64;
        let dp_next = dp_prev + (2 * k + 1) as f64 * p_curr;
        
        p_prev = p_curr;
        p_curr = p_next;
        dp_prev = dp_curr;
        dp_curr = dp_next;
    }
    
    (p_curr, dp_curr)
}

/// Compute Legendre-Gauss-Lobatto weights
fn compute_lgl_weights(nodes: &[f64], n: usize) -> Vec<f64> {
    let mut weights = vec![0.0; nodes.len()];
    
    for (i, &x) in nodes.iter().enumerate() {
        let pn = legendre_poly(n, x);
        weights[i] = 2.0 / (n as f64 * (n + 1) as f64) / (pn * pn);
    }
    
    weights
}

/// Compute the derivative matrix for Legendre-Gauss-Lobatto points
fn compute_derivative_matrix(nodes: &[f64], n: usize) -> DMatrix<f64> {
    let np = nodes.len();
    let mut d = DMatrix::zeros(np, np);
    
    for i in 0..np {
        let pi = legendre_poly(n, nodes[i]);
        for j in 0..np {
            if i != j {
                let pj = legendre_poly(n, nodes[j]);
                d[(i, j)] = (pi / pj) / (nodes[i] - nodes[j]);
            }
        }
    }
    
    // Diagonal entries
    d[(0, 0)] = -(n as f64 * (n + 1) as f64) / 4.0;
    d[(np - 1, np - 1)] = (n as f64 * (n + 1) as f64) / 4.0;
    
    // For interior points, the diagonal is 0.0
    // But it's more robust to enforce sum to zero property: d_ii = -sum_{j != i} d_ij
    for i in 1..np - 1 {
        let mut sum = 0.0;
        for j in 0..np {
            if i != j {
                sum += d[(i, j)];
            }
        }
        d[(i, i)] = -sum;
    }
    
    d
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use nalgebra::DVector;
    
    #[test]
    fn test_legendre_poly() {
        // Test Legendre polynomial values
        assert_eq!(legendre_poly(0, 0.5), 1.0);
        assert_eq!(legendre_poly(1, 0.5), 0.5);
        assert_relative_eq!(legendre_poly(2, 0.5), -0.125);
        assert_relative_eq!(legendre_poly(3, 0.5), -0.4375);
    }
    
    #[test]
    fn test_lgl_nodes() {
        // Test LGL nodes for n=2 (should be -1, 0, 1)
        let nodes = compute_lgl_nodes(2).unwrap();
        assert_relative_eq!(nodes[0], -1.0, epsilon = 1e-10);
        assert_relative_eq!(nodes[1], 0.0, epsilon = 1e-10);
        assert_relative_eq!(nodes[2], 1.0, epsilon = 1e-10);
        
        // Test symmetry for n=4
        let nodes = compute_lgl_nodes(4).unwrap();
        assert_relative_eq!(nodes[0], -1.0, epsilon = 1e-10);
        assert_relative_eq!(nodes[4], 1.0, epsilon = 1e-10);
        assert_relative_eq!(nodes[1], -nodes[3], epsilon = 1e-10);
        assert_relative_eq!(nodes[2], 0.0, epsilon = 1e-10);
    }
    
    #[test]
    fn test_lgl_weights() {
        let nodes = compute_lgl_nodes(2).unwrap();
        let weights = compute_lgl_weights(&nodes, 2);
        
        // Weights should be [1/3, 4/3, 1/3] for n=2
        assert_relative_eq!(weights[0], 1.0 / 3.0, epsilon = 1e-10);
        assert_relative_eq!(weights[1], 4.0 / 3.0, epsilon = 1e-10);
        assert_relative_eq!(weights[2], 1.0 / 3.0, epsilon = 1e-10);
        
        // Sum of weights should be 2.0 (integral of 1 from -1 to 1)
        assert_relative_eq!(weights.iter().sum::<f64>(), 2.0, epsilon = 1e-10);
    }
    
    #[test]
    fn test_derivative_matrix() {
        let n = 4;
        let nodes = compute_lgl_nodes(n).unwrap();
        let d = compute_derivative_matrix(&nodes, n);
        
        // Derivative of constant should be zero
        let ones = DVector::from_element(n + 1, 1.0);
        let deriv = &d * &ones;
        for &val in deriv.iter() {
            assert_relative_eq!(val, 0.0, epsilon = 1e-10);
        }
        
        // Derivative of x should be 1
        let x = DVector::from_vec(nodes.clone());
        let deriv = &d * &x;
        for &val in deriv.iter() {
            assert_relative_eq!(val, 1.0, epsilon = 1e-10);
        }
    }
}
