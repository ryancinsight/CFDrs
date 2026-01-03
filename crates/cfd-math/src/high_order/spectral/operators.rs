//! Spectral differential operators and numerical integration.
//!
//! This module provides implementations of differential operators and numerical
//! integration for spectral element methods.

use super::{SpectralElement, Result, SpectralError, compute_lgl_nodes, compute_lgl_weights};
use nalgebra::{DMatrix, DVector};

/// Spectral differentiation operator
pub struct SpectralDiffOp {
    /// Differentiation matrix
    d_matrix: DMatrix<f64>,
    /// Number of nodes
    num_nodes: usize,
}

impl SpectralDiffOp {
    /// Create a new spectral differentiation operator
    pub fn new(element: &SpectralElement) -> Self {
        Self {
            d_matrix: element.derivative_matrix.clone(),
            num_nodes: element.num_nodes,
        }
    }
    
    /// Apply the differentiation operator to a function
    pub fn apply(&self, u: &DVector<f64>) -> DVector<f64> {
        &self.d_matrix * u
    }
    
    /// Compute the gradient of a scalar field
    pub fn gradient(&self, u: &DVector<f64>) -> DVector<f64> {
        self.apply(u)
    }
    
    /// Compute the divergence of a vector field
    pub fn divergence(&self, u: &[DVector<f64>]) -> DVector<f64> {
        let mut div = DVector::zeros(self.num_nodes);
        
        for u_comp in u {
            div += &self.apply(u_comp);
        }
        
        div
    }
    
    /// Compute the curl of a 2D vector field
    pub fn curl_2d(&self, u: &[DVector<f64>; 2]) -> DVector<f64> {
        assert_eq!(u[0].len(), self.num_nodes);
        assert_eq!(u[1].len(), self.num_nodes);
        
        // In 2D, curl(u,v) = ∂v/∂x - ∂u/∂y
        let du_dy = self.apply(&u[0]);
        let dv_dx = self.apply(&u[1]);
        
        dv_dx - du_dy
    }
    
    /// Compute the Laplacian of a scalar field
    pub fn laplacian(&self, u: &DVector<f64>) -> DVector<f64> {
        self.apply(&self.apply(u))
    }
}

/// Spectral interpolation operator
pub struct SpectralInterp {
    /// Reference element
    element: SpectralElement,
}

impl SpectralInterp {
    /// Create a new spectral interpolation operator
    pub fn new(element: &SpectralElement) -> Self {
        Self {
            element: element.clone(),
        }
    }
    
    /// Interpolate a function to the spectral element nodes
    pub fn interpolate<F>(&self, f: F) -> DVector<f64>
    where
        F: Fn(f64) -> f64,
    {
        self.element.interpolate(f)
    }
    
    /// Evaluate the interpolant at a point in the reference element
    pub fn eval_at(&self, x: f64, u: &DVector<f64>) -> f64 {
        let mut result = 0.0;
        
        for i in 0..self.element.num_nodes {
            result += u[i] * self.element.lagrange_basis(i, x);
        }
        
        result
    }
    
    /// Compute the L2 projection of a function onto the spectral element space
    pub fn l2_projection<F>(&self, f: F) -> DVector<f64>
    where
        F: Fn(f64) -> f64 + Sync,
    {
        let mut mass = DMatrix::zeros(self.element.num_nodes, self.element.num_nodes);
        let mut rhs = DVector::zeros(self.element.num_nodes);
        
        // Use Gauss-Lobatto quadrature for mass matrix and right-hand side
        for i in 0..self.element.num_nodes {
            for j in 0..self.element.num_nodes {
                let mut mij = 0.0;
                
                for (k, &xk) in self.element.nodes.iter().enumerate() {
                    let phi_i = self.element.lagrange_basis(i, xk);
                    let phi_j = self.element.lagrange_basis(j, xk);
                    mij += self.element.weights[k] * phi_i * phi_j;
                }
                
                mass[(i, j)] = mij;
            }
            
            // Compute right-hand side
            let mut fi = 0.0;
            
            for (k, &xk) in self.element.nodes.iter().enumerate() {
                let phi_i = self.element.lagrange_basis(i, xk);
                fi += self.element.weights[k] * f(xk) * phi_i;
            }
            
            rhs[i] = fi;
        }
        
        // Solve the linear system M u = f
        mass.lu().solve(&rhs).unwrap_or_else(|| {
            // Fall back to interpolation if the mass matrix is singular
            self.interpolate(f)
        })
    }
}

/// Spectral quadrature for numerical integration
pub struct SpectralQuadrature {
    /// Quadrature weights
    weights: DVector<f64>,
    /// Quadrature nodes
    nodes: DVector<f64>,
    /// Number of quadrature points
    num_points: usize,
}

impl SpectralQuadrature {
    /// Create a new spectral quadrature rule
    pub fn new(weights: DVector<f64>, nodes: DVector<f64>) -> Self {
        Self {
            num_points: weights.len(),
            weights,
            nodes,
        }
    }
    
    /// Create a Gauss-Lobatto quadrature rule with n points
    pub fn gauss_lobatto(n: usize) -> Result<Self> {
        if n < 2 {
            return Err(SpectralError::InvalidOrder(n));
        }
        
        let nodes_vec = compute_lgl_nodes(n - 1)?; // n points for order 2n-3
        let weights = DVector::from_vec(compute_lgl_weights(&nodes_vec, n - 1));
        let nodes = DVector::from_vec(nodes_vec);
        
        Ok(Self {
            num_points: n,
            weights,
            nodes,
        })
    }
    
    /// Integrate a function over the reference interval [-1, 1]
    pub fn integrate<F>(&self, f: F) -> f64
    where
        F: Fn(f64) -> f64,
    {
        let mut integral = 0.0;
        
        for i in 0..self.num_points {
            let x = self.nodes[i];
            let w = self.weights[i];
            integral += w * f(x);
        }
        
        integral
    }
    
    /// Integrate the product of two functions
    pub fn integrate_product<F, G>(&self, f: F, g: G) -> f64
    where
        F: Fn(f64) -> f64,
        G: Fn(f64) -> f64,
    {
        self.integrate(|x| f(x) * g(x))
    }
    
    /// Compute the L2 inner product of two functions
    pub fn l2_inner_product<F, G>(&self, f: F, g: G) -> f64
    where
        F: Fn(f64) -> f64,
        G: Fn(f64) -> f64,
    {
        self.integrate_product(f, g)
    }
    
    /// Compute the L2 norm of a function
    pub fn l2_norm<F>(&self, f: F) -> f64
    where
        F: Fn(f64) -> f64,
    {
        self.l2_inner_product(&f, &f).sqrt()
    }
}

/// Spectral filter for stabilization
pub struct SpectralFilter {
    /// Filter coefficients
    filter_coeffs: DVector<f64>,
    /// Filter order
    #[allow(dead_code)]
    order: usize,
}

impl SpectralFilter {
    /// Create a new spectral filter
    pub fn new(order: usize, alpha: f64, p: usize) -> Self {
        let n = order + 1;
        let mut sigma = DVector::zeros(n);
        
        // Exponential filter
        for i in 0..n {
            let eta = i as f64 / order as f64;
            sigma[i] = (-alpha * eta.powi(p as i32)).exp();
        }
        
        Self {
            filter_coeffs: sigma,
            order,
        }
    }
    
    /// Apply the filter to a function
    pub fn apply(&self, u: &DVector<f64>) -> DVector<f64> {
        // In a real implementation, this would involve a forward and inverse transform
        // For simplicity, we'll just apply a pointwise multiplication here
        // A proper implementation would use FFT or similar
        
        let mut result = u.clone();
        for i in 0..result.len() {
            if i < self.filter_coeffs.len() {
                result[i] *= self.filter_coeffs[i];
            }
        }
        
        result
    }
}

/// Spectral time integration methods
pub mod time_integration {
    use super::DVector;
    
    /// Runge-Kutta time integration for spectral methods
    pub struct RungeKutta {
        /// Butcher tableau coefficients
        a: Vec<Vec<f64>>,
        /// Time step coefficients
        b: Vec<f64>,
        /// Node coefficients
        c: Vec<f64>,
        /// Number of stages
        stages: usize,
    }
    
    impl RungeKutta {
        /// Create a new Runge-Kutta method from Butcher tableau
        pub fn new(a: Vec<Vec<f64>>, b: Vec<f64>, c: Vec<f64>) -> Self {
            let stages = b.len();
            
            // Pad the Butcher tableau if necessary
            let mut a_padded = a;
            for row in &mut a_padded {
                while row.len() < stages {
                    row.push(0.0);
                }
            }
            
            Self {
                a: a_padded,
                b,
                c,
                stages,
            }
        }
        
        /// Create the classical 4th-order Runge-Kutta method
        pub fn rk4() -> Self {
            let a = vec![
                vec![],
                vec![0.5],
                vec![0.0, 0.5],
                vec![0.0, 0.0, 1.0],
            ];
            
            let b = vec![1.0/6.0, 1.0/3.0, 1.0/3.0, 1.0/6.0];
            let c = vec![0.0, 0.5, 0.5, 1.0];
            
            Self::new(a, b, c)
        }
        
        /// Perform one time step
        pub fn step<F>(
            &self,
            f: F,
            t: f64,
            dt: f64,
            u: &DVector<f64>,
        ) -> DVector<f64>
        where
            F: Fn(f64, &DVector<f64>) -> DVector<f64>,
        {
            let mut k = Vec::with_capacity(self.stages);
            let mut u_new = u.clone();
            
            // Compute stage values
            for i in 0..self.stages {
                let mut sum = DVector::zeros(u.len());
                
                for j in 0..i {
                    sum += &k[j] * self.a[i][j];
                }
                
                let t_stage = t + self.c[i] * dt;
                let u_stage = u + &sum * dt;
                
                k.push(f(t_stage, &u_stage));
            }
            
            // Combine stages
            for i in 0..self.stages {
                u_new += &k[i] * (dt * self.b[i]);
            }
            
            u_new
        }
    }
    
    /// Adams-Bashforth time integration for spectral methods
    pub struct AdamsBashforth {
        /// Coefficients
        coeffs: Vec<f64>,
        /// Order of the method
        order: usize,
    }
    
    impl AdamsBashforth {
        /// Create a new Adams-Bashforth method of given order
        pub fn new(order: usize) -> Self {
            // Coefficients for Adams-Bashforth methods
            let coeffs = match order {
                1 => vec![1.0],
                2 => vec![3.0/2.0, -1.0/2.0],
                3 => vec![23.0/12.0, -16.0/12.0, 5.0/12.0],
                4 => vec![55.0/24.0, -59.0/24.0, 37.0/24.0, -9.0/24.0],
                _ => panic!("Unsupported Adams-Bashforth order"),
            };
            
            Self { coeffs, order }
        }
        
        /// Perform one time step
        pub fn step<F>(
            &self,
            f: F,
            t: f64,
            dt: f64,
            u: &DVector<f64>,
            f_prev: &[DVector<f64>],
        ) -> (DVector<f64>, Vec<DVector<f64>>)
        where
            F: Fn(f64, &DVector<f64>) -> DVector<f64>,
        {
            let f_new = f(t, u);
            let mut sum = &f_new * self.coeffs[0];
            
            for (i, f_prev_i) in f_prev.iter().take(self.order - 1).enumerate() {
                sum += f_prev_i * self.coeffs[i + 1];
            }
            
            let u_new = u + &sum * dt;
            
            // Update history of derivatives
            let mut f_history = Vec::with_capacity(self.order);
            f_history.push(f_new);
            f_history.extend_from_slice(&f_prev[..self.order - 1]);
            
            (u_new, f_history)
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    
    #[test]
    fn test_spectral_diff_op() {
        let order = 4;
        let element = SpectralElement::new(order).unwrap();
        let diff_op = SpectralDiffOp::new(&element);
        
        // Test derivative of a linear function
        let a = 2.0;
        let b = 3.0;
        let f = |x: f64| a * x + b;
        let df = |_: f64| a;
        
        let u = element.interpolate(f);
        let du = diff_op.apply(&u);
        
        for i in 0..element.num_nodes {
            assert_relative_eq!(du[i], df(element.nodes[i]), epsilon = 1e-10);
        }
        
        // Test derivative of a quadratic function
        let c = 1.5;
        let f = |x: f64| a * x * x + b * x + c;
        let df = |x: f64| 2.0 * a * x + b;
        
        let u = element.interpolate(f);
        let du = diff_op.apply(&u);
        
        // The derivative should be exact for polynomials of degree <= order
        for i in 0..element.num_nodes {
            assert_relative_eq!(du[i], df(element.nodes[i]), epsilon = 1e-10);
        }
    }
    
    #[test]
    fn test_spectral_interp() {
        let order = 4;
        let element = SpectralElement::new(order).unwrap();
        let interp = SpectralInterp::new(&element);
        
        // Test interpolation of a polynomial
        let f = |x: f64| x.powi(3) - 2.0 * x + 1.0;
        let u = interp.interpolate(f);
        
        // The interpolation should be exact at the nodes
        for (i, &x) in element.nodes.iter().enumerate() {
            assert_relative_eq!(u[i], f(x), epsilon = 1e-10);
        }
        
        // Test evaluation at a point
        let x_test = 0.5;
        let u_test = interp.eval_at(x_test, &u);
        assert_relative_eq!(u_test, f(x_test), epsilon = 1e-10);
    }
    
    #[test]
    fn test_spectral_quadrature() {
        let order = 4;
        let quad = SpectralQuadrature::gauss_lobatto(order + 1).unwrap();
        
        // Test integration of a polynomial
        let f = |x: f64| x.powi(6);  // Degree 6
        let integral = quad.integrate(f);
        let exact = 2.0 / 7.0;  // ∫_{-1}^1 x^6 dx = 2/7
        
        // The quadrature should be exact for polynomials of degree <= 2n-3
        // For n=5 (order=4), this means degree <= 7
        assert_relative_eq!(integral, exact, epsilon = 1e-10);
    }
    
    #[test]
    fn test_runge_kutta() {
        // Test the RK4 method on a simple ODE: du/dt = -u, u(0) = 1
        // Exact solution: u(t) = exp(-t)
        let rk = time_integration::RungeKutta::rk4();
        let f = |_t: f64, u: &DVector<f64>| -u.clone();

        let t0: f64 = 0.0;
        let t_final: f64 = 1.0;
        let dt: f64 = 0.1;

        let mut t = t0;
        let mut u = DVector::from_vec(vec![1.0]);

        while t < t_final - 1e-10 {
            let dt_step = dt.min(t_final - t);
            u = rk.step(&f, t, dt_step, &u);
            t += dt_step;
        }

        let u_exact = (-t_final).exp();
        assert_relative_eq!(u[0], u_exact, epsilon = 1e-5);
    }
}
