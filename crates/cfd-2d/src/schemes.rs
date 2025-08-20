//! Numerical schemes for 2D CFD simulations
//!
//! This module provides various discretization schemes for spatial and temporal
//! derivatives in 2D computational fluid dynamics.

use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;
use cfd_core::error::Error;
use serde::{Deserialize, Serialize};

/// Named constants for scheme parameters
// Unused DEFAULT_CFL_NUMBER constant removed (CFL handled in cfd-core::constants)
const MIN_LIMITER_THRESHOLD: f64 = 1e-10;
const WENO_EPSILON: f64 = 1e-6;
const WENO_POWER: i32 = 2;

/// Spatial discretization scheme
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum SpatialScheme {
    /// First-order upwind
    FirstOrderUpwind,
    /// Second-order central differencing
    CentralDifference,
    /// Second-order upwind
    SecondOrderUpwind,
    /// QUICK (Quadratic Upstream Interpolation)
    Quick,
    /// Third-order MUSCL
    Muscl,
    /// Fifth-order WENO
    Weno5,
    /// Fourth-order explicit central difference
    FourthOrderCentral,
}

/// Flux limiter for TVD schemes
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum FluxLimiter {
    /// No limiter (unlimited)
    None,
    /// Moderate limiter (preventing oscillations)
    Moderate,
    /// Van Leer limiter
    VanLeer,
    /// Van Albada limiter
    VanAlbada,
    /// Superbee limiter
    Superbee,
    /// MC (Monotonized Central) limiter
    MC,
}

/// Time integration scheme
#[derive(Debug, Clone, Copy, PartialEq, Serialize, Deserialize)]
pub enum TimeScheme {
    /// Forward Euler (explicit)
    ForwardEuler,
    /// Backward Euler (implicit)
    BackwardEuler,
    /// Crank-Nicolson
    CrankNicolson,
    /// Second-order Runge-Kutta
    RungeKutta2,
    /// Fourth-order Runge-Kutta
    RungeKutta4,
    /// Adams-Bashforth (2nd order)
    AdamsBashforth2,
}

/// 2D grid for finite difference operations
#[derive(Debug, Clone)]
pub struct Grid2D<T: RealField> {
    /// Grid values
    pub data: DMatrix<T>,
    /// Grid spacing in x-direction
    pub dx: T,
    /// Grid spacing in y-direction
    pub dy: T,
    /// Number of ghost cells
    pub ghost_cells: usize,
}

impl<T: RealField + FromPrimitive> Grid2D<T> {
    /// Create a new 2D grid
    pub fn new(nx: usize, ny: usize, dx: T, dy: T, ghost_cells: usize) -> Self {
        let total_nx = nx + 2 * ghost_cells;
        let total_ny = ny + 2 * ghost_cells;
        
        Self {
            data: DMatrix::zeros(total_nx, total_ny),
            dx,
            dy,
            ghost_cells,
        }
    }
    
    /// Get interior dimensions (excluding ghost cells)
    pub fn interior_shape(&self) -> (usize, usize) {
        let (total_nx, total_ny) = self.data.shape();
        (
            total_nx - 2 * self.ghost_cells,
            total_ny - 2 * self.ghost_cells,
        )
    }
    
    /// Apply periodic boundary conditions
    pub fn apply_periodic_bc(&mut self) {
        let (nx, ny) = self.data.shape();
        let g = self.ghost_cells;
        
        // Periodic in x
        for j in 0..ny {
            for k in 0..g {
                self.data[(k, j)] = self.data[(nx - 2*g + k, j)].clone();
                self.data[(nx - g + k, j)] = self.data[(g + k, j)].clone();
            }
        }
        
        // Periodic in y
        for i in 0..nx {
            for k in 0..g {
                self.data[(i, k)] = self.data[(i, ny - 2*g + k)].clone();
                self.data[(i, ny - g + k)] = self.data[(i, g + k)].clone();
            }
        }
    }
    
    /// Apply Dirichlet boundary conditions
    pub fn apply_dirichlet_bc(&mut self, value: T) {
        let (nx, ny) = self.data.shape();
        let g = self.ghost_cells;
        
        // Set ghost cells to boundary value
        for j in 0..ny {
            for k in 0..g {
                self.data[(k, j)] = value.clone();
                self.data[(nx - 1 - k, j)] = value.clone();
            }
        }
        
        for i in 0..nx {
            for k in 0..g {
                self.data[(i, k)] = value.clone();
                self.data[(i, ny - 1 - k)] = value.clone();
            }
        }
    }
    
    /// Apply Neumann boundary conditions (zero gradient)
    pub fn apply_neumann_bc(&mut self) {
        let (nx, ny) = self.data.shape();
        let g = self.ghost_cells;
        
        // Copy interior values to ghost cells
        for j in 0..ny {
            for k in 0..g {
                self.data[(k, j)] = self.data[(g, j)].clone();
                self.data[(nx - 1 - k, j)] = self.data[(nx - g - 1, j)].clone();
            }
        }
        
        for i in 0..nx {
            for k in 0..g {
                self.data[(i, k)] = self.data[(i, g)].clone();
                self.data[(i, ny - 1 - k)] = self.data[(i, ny - g - 1)].clone();
            }
        }
    }
}

/// Finite difference operator for spatial derivatives
pub struct FiniteDifference<T: RealField> {
    scheme: SpatialScheme,
    limiter: Option<FluxLimiter>,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive> FiniteDifference<T> {
    /// Create a new finite difference operator
    pub fn new(scheme: SpatialScheme) -> Self {
        Self {
            scheme,
            limiter: None,
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Set flux limiter for TVD schemes
    pub fn with_limiter(mut self, limiter: FluxLimiter) -> Self {
        self.limiter = Some(limiter);
        self
    }
    
    /// Compute x-derivative using specified scheme
    /// For upwind schemes, velocity is required to determine direction
    pub fn ddx(&self, grid: &Grid2D<T>, i: usize, j: usize, u_velocity: Option<T>) -> T {
        match self.scheme {
            SpatialScheme::FirstOrderUpwind => {
                let u = u_velocity.unwrap_or_else(|| T::zero());
                self.first_order_upwind_x(grid, i, j, u)
            }
            SpatialScheme::CentralDifference => self.central_difference_x(grid, i, j),
            SpatialScheme::SecondOrderUpwind => {
                let u = u_velocity.unwrap_or_else(|| T::zero());
                self.second_order_upwind_x(grid, i, j, u)
            }
            SpatialScheme::Quick => {
                let u = u_velocity.unwrap_or_else(|| T::zero());
                self.quick_x(grid, i, j, u)
            }
            SpatialScheme::Muscl => {
                let u = u_velocity.unwrap_or_else(|| T::zero());
                self.muscl_x(grid, i, j, u)
            }
            SpatialScheme::Weno5 => self.weno5_x(grid, i, j),
            SpatialScheme::FourthOrderCentral => self.fourth_order_central_x(grid, i, j),
        }
    }
    
    /// Compute y-derivative using specified scheme
    /// For upwind schemes, velocity is required to determine direction
    pub fn ddy(&self, grid: &Grid2D<T>, i: usize, j: usize, v_velocity: Option<T>) -> T {
        match self.scheme {
            SpatialScheme::FirstOrderUpwind => {
                let v = v_velocity.unwrap_or_else(|| T::zero());
                self.first_order_upwind_y(grid, i, j, v)
            }
            SpatialScheme::CentralDifference => self.central_difference_y(grid, i, j),
            SpatialScheme::SecondOrderUpwind => {
                let v = v_velocity.unwrap_or_else(|| T::zero());
                self.second_order_upwind_y(grid, i, j, v)
            }
            SpatialScheme::Quick => {
                let v = v_velocity.unwrap_or_else(|| T::zero());
                self.quick_y(grid, i, j, v)
            }
            SpatialScheme::Muscl => {
                let v = v_velocity.unwrap_or_else(|| T::zero());
                self.muscl_y(grid, i, j, v)
            }
            SpatialScheme::Weno5 => self.weno5_y(grid, i, j),
            SpatialScheme::FourthOrderCentral => self.fourth_order_central_y(grid, i, j),
        }
    }
    
    /// First-order upwind in x-direction with velocity-dependent stencil
    fn first_order_upwind_x(&self, grid: &Grid2D<T>, i: usize, j: usize, u: T) -> T {
        if u >= T::zero() {
            // Positive velocity: backward difference
            if i > 0 {
                (grid.data[(i, j)].clone() - grid.data[(i-1, j)].clone()) / grid.dx.clone()
            } else {
                // Use forward difference at boundary
                (grid.data[(i+1, j)].clone() - grid.data[(i, j)].clone()) / grid.dx.clone()
            }
        } else {
            // Negative velocity: forward difference
            if i < grid.data.nrows() - 1 {
                (grid.data[(i+1, j)].clone() - grid.data[(i, j)].clone()) / grid.dx.clone()
            } else {
                // Use backward difference at boundary
                (grid.data[(i, j)].clone() - grid.data[(i-1, j)].clone()) / grid.dx.clone()
            }
        }
    }
    
    /// First-order upwind in y-direction with velocity-dependent stencil
    fn first_order_upwind_y(&self, grid: &Grid2D<T>, i: usize, j: usize, v: T) -> T {
        if v >= T::zero() {
            // Positive velocity: backward difference
            if j > 0 {
                (grid.data[(i, j)].clone() - grid.data[(i, j-1)].clone()) / grid.dy.clone()
            } else {
                // Use forward difference at boundary
                (grid.data[(i, j+1)].clone() - grid.data[(i, j)].clone()) / grid.dy.clone()
            }
        } else {
            // Negative velocity: forward difference
            if j < grid.data.ncols() - 1 {
                (grid.data[(i, j+1)].clone() - grid.data[(i, j)].clone()) / grid.dy.clone()
            } else {
                // Use backward difference at boundary
                (grid.data[(i, j)].clone() - grid.data[(i, j-1)].clone()) / grid.dy.clone()
            }
        }
    }
    
    /// Central difference in x-direction (velocity-independent)
    fn central_difference_x(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        if i == 0 {
            // Forward difference at left boundary
            (grid.data[(i+1, j)].clone() - grid.data[(i, j)].clone()) / grid.dx.clone()
        } else if i == grid.data.nrows() - 1 {
            // Backward difference at right boundary
            (grid.data[(i, j)].clone() - grid.data[(i-1, j)].clone()) / grid.dx.clone()
        } else {
            // Central difference in interior
            (grid.data[(i+1, j)].clone() - grid.data[(i-1, j)].clone()) / 
            (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * grid.dx.clone())
        }
    }
    
    /// Central difference in y-direction (velocity-independent)
    fn central_difference_y(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        if j == 0 {
            // Forward difference at bottom boundary
            (grid.data[(i, j+1)].clone() - grid.data[(i, j)].clone()) / grid.dy.clone()
        } else if j == grid.data.ncols() - 1 {
            // Backward difference at top boundary
            (grid.data[(i, j)].clone() - grid.data[(i, j-1)].clone()) / grid.dy.clone()
        } else {
            // Central difference in interior
            (grid.data[(i, j+1)].clone() - grid.data[(i, j-1)].clone()) / 
            (T::from_f64(2.0).unwrap_or_else(|| T::zero()) * grid.dy.clone())
        }
    }
    
    /// Second-order upwind in x-direction with velocity-dependent stencil
    fn second_order_upwind_x(&self, grid: &Grid2D<T>, i: usize, j: usize, u: T) -> T {
        let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        
        if u >= T::zero() {
            // Positive velocity: use points i-2, i-1, i
            if i >= 2 {
                (three * grid.data[(i, j)].clone() - 
                 four * grid.data[(i-1, j)].clone() + 
                 grid.data[(i-2, j)].clone()) / (two * grid.dx.clone())
            } else {
                // Fall back to first-order at boundary
                self.first_order_upwind_x(grid, i, j, u)
            }
        } else {
            // Negative velocity: use points i, i+1, i+2
            if i < grid.data.nrows() - 2 {
                (-three * grid.data[(i, j)].clone() + 
                 four * grid.data[(i+1, j)].clone() - 
                 grid.data[(i+2, j)].clone()) / (two * grid.dx.clone())
            } else {
                // Fall back to first-order at boundary
                self.first_order_upwind_x(grid, i, j, u)
            }
        }
    }
    
    /// Second-order upwind in y-direction with velocity-dependent stencil
    fn second_order_upwind_y(&self, grid: &Grid2D<T>, i: usize, j: usize, v: T) -> T {
        let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        
        if v >= T::zero() {
            // Positive velocity: use points j-2, j-1, j
            if j >= 2 {
                (three * grid.data[(i, j)].clone() - 
                 four * grid.data[(i, j-1)].clone() + 
                 grid.data[(i, j-2)].clone()) / (two * grid.dy.clone())
            } else {
                // Fall back to first-order at boundary
                self.first_order_upwind_y(grid, i, j, v)
            }
        } else {
            // Negative velocity: use points j, j+1, j+2
            if j < grid.data.ncols() - 2 {
                (-three * grid.data[(i, j)].clone() + 
                 four * grid.data[(i, j+1)].clone() - 
                 grid.data[(i, j+2)].clone()) / (two * grid.dy.clone())
            } else {
                // Fall back to first-order at boundary
                self.first_order_upwind_y(grid, i, j, v)
            }
        }
    }
    
    /// QUICK scheme in x-direction with proper velocity-dependent upwinding
    /// Reference: Leonard, B.P. (1979) "A stable and accurate convective modelling procedure"
    fn quick_x(&self, grid: &Grid2D<T>, i: usize, j: usize, u: T) -> T {
        let six_eighths = T::from_f64(0.75).unwrap_or_else(|| T::zero());
        let three_eighths = T::from_f64(0.375).unwrap_or_else(|| T::zero());
        let one_eighth = T::from_f64(0.125).unwrap_or_else(|| T::zero());
        
        if u >= T::zero() {
            // Positive flow: upstream is i-2, central is i-1, downstream is i
            if i >= 2 && i < grid.data.nrows() - 1 {
                let phi_u = grid.data[(i-2, j)].clone();
                let phi_c = grid.data[(i-1, j)].clone();
                let phi_d = grid.data[(i, j)].clone();
                
                let phi_face = six_eighths.clone() * phi_d + three_eighths.clone() * phi_c - one_eighth.clone() * phi_u;
                
                // Compute derivative using the interpolated face value
                (grid.data[(i+1, j)].clone() - phi_face) / grid.dx.clone()
            } else {
                // Fall back to second-order upwind at boundaries
                self.second_order_upwind_x(grid, i, j, u)
            }
        } else {
            // Negative flow: upstream is i+2, central is i+1, downstream is i
            if i < grid.data.nrows() - 2 && i > 0 {
                let phi_u = grid.data[(i+2, j)].clone();
                let phi_c = grid.data[(i+1, j)].clone();
                let phi_d = grid.data[(i, j)].clone();
                
                let phi_face = six_eighths * phi_d + three_eighths * phi_c - one_eighth * phi_u;
                
                // Compute derivative using the interpolated face value
                (phi_face - grid.data[(i-1, j)].clone()) / grid.dx.clone()
            } else {
                // Fall back to second-order upwind at boundaries
                self.second_order_upwind_x(grid, i, j, u)
            }
        }
    }
    
    /// QUICK scheme in y-direction with proper velocity-dependent upwinding
    fn quick_y(&self, grid: &Grid2D<T>, i: usize, j: usize, v: T) -> T {
        let six_eighths = T::from_f64(0.75).unwrap_or_else(|| T::zero());
        let three_eighths = T::from_f64(0.375).unwrap_or_else(|| T::zero());
        let one_eighth = T::from_f64(0.125).unwrap_or_else(|| T::zero());
        
        if v >= T::zero() {
            // Positive flow: upstream is j-2, central is j-1, downstream is j
            if j >= 2 && j < grid.data.ncols() - 1 {
                let phi_u = grid.data[(i, j-2)].clone();
                let phi_c = grid.data[(i, j-1)].clone();
                let phi_d = grid.data[(i, j)].clone();
                
                let phi_face = six_eighths.clone() * phi_d + three_eighths.clone() * phi_c - one_eighth.clone() * phi_u;
                
                // Compute derivative using the interpolated face value
                (grid.data[(i, j+1)].clone() - phi_face) / grid.dy.clone()
            } else {
                // Fall back to second-order upwind at boundaries
                self.second_order_upwind_y(grid, i, j, v)
            }
        } else {
            // Negative flow: upstream is j+2, central is j+1, downstream is j
            if j < grid.data.ncols() - 2 && j > 0 {
                let phi_u = grid.data[(i, j+2)].clone();
                let phi_c = grid.data[(i, j+1)].clone();
                let phi_d = grid.data[(i, j)].clone();
                
                let phi_face = six_eighths * phi_d + three_eighths * phi_c - one_eighth * phi_u;
                
                // Compute derivative using the interpolated face value
                (phi_face - grid.data[(i, j-1)].clone()) / grid.dy.clone()
            } else {
                // Fall back to second-order upwind at boundaries
                self.second_order_upwind_y(grid, i, j, v)
            }
        }
    }
    
    /// MUSCL scheme in x-direction with velocity-dependent upwinding
    fn muscl_x(&self, grid: &Grid2D<T>, i: usize, j: usize, u: T) -> T {
        if u >= T::zero() {
            // Positive velocity
            let r = self.compute_gradient_ratio_x(grid, i, j);
            let phi = self.apply_limiter(r);
            
            let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
            let base = (grid.data[(i+1, j)].clone() - grid.data[(i, j)].clone()) / grid.dx.clone();
            let correction = phi * (grid.data[(i, j)].clone() - grid.data[(i-1, j)].clone()) / grid.dx.clone();
            
            base + half * correction
        } else {
            // Negative velocity - reverse the stencil
            let r = self.compute_gradient_ratio_x_negative(grid, i, j);
            let phi = self.apply_limiter(r);
            
            let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
            let base = (grid.data[(i, j)].clone() - grid.data[(i-1, j)].clone()) / grid.dx.clone();
            let correction = phi * (grid.data[(i+1, j)].clone() - grid.data[(i, j)].clone()) / grid.dx.clone();
            
            base - half * correction
        }
    }
    
    /// MUSCL scheme in y-direction with velocity-dependent upwinding
    fn muscl_y(&self, grid: &Grid2D<T>, i: usize, j: usize, v: T) -> T {
        if v >= T::zero() {
            // Positive velocity
            let r = self.compute_gradient_ratio_y(grid, i, j);
            let phi = self.apply_limiter(r);
            
            let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
            let base = (grid.data[(i, j+1)].clone() - grid.data[(i, j)].clone()) / grid.dy.clone();
            let correction = phi * (grid.data[(i, j)].clone() - grid.data[(i, j-1)].clone()) / grid.dy.clone();
            
            base + half * correction
        } else {
            // Negative velocity - reverse the stencil
            let r = self.compute_gradient_ratio_y_negative(grid, i, j);
            let phi = self.apply_limiter(r);
            
            let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
            let base = (grid.data[(i, j)].clone() - grid.data[(i, j-1)].clone()) / grid.dy.clone();
            let correction = phi * (grid.data[(i, j+1)].clone() - grid.data[(i, j)].clone()) / grid.dy.clone();
            
            base - half * correction
        }
    }
    
    /// WENO5 in x-direction
    /// Reference: Shu, C.W. (2009) "High order weighted essentially non-oscillatory schemes"
    fn weno5_x(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        if i < 2 || i >= grid.data.nrows() - 2 {
            // Fall back to central difference at boundaries
            return self.central_difference_x(grid, i, j);
        }
        
        // WENO5 stencil values
        let v = [
            grid.data[(i-2, j)].clone(),
            grid.data[(i-1, j)].clone(),
            grid.data[(i, j)].clone(),
            grid.data[(i+1, j)].clone(),
            grid.data[(i+2, j)].clone(),
        ];
        
        // Three sub-stencils for reconstruction
        let s0 = (v[0].clone() * T::from_f64(-1.0).unwrap_or_else(|| T::zero()) + 
                  v[1].clone() * T::from_f64(3.0).unwrap_or_else(|| T::zero()) - 
                  v[2].clone() * T::from_f64(3.0).unwrap_or_else(|| T::zero()) + 
                  v[3].clone()) / T::from_f64(6.0).unwrap_or_else(|| T::zero());
                  
        let s1 = (v[1].clone() * T::from_f64(-1.0).unwrap_or_else(|| T::zero()) + 
                  v[3].clone()) / T::from_f64(2.0).unwrap_or_else(|| T::zero());
                  
        let s2 = (v[2].clone() * T::from_f64(-3.0).unwrap_or_else(|| T::zero()) + 
                  v[3].clone() * T::from_f64(3.0).unwrap_or_else(|| T::zero()) + 
                  v[4].clone() * T::from_f64(-1.0).unwrap_or_else(|| T::zero())) / T::from_f64(6.0).unwrap_or_else(|| T::zero());
        
        // Compute smoothness indicators
        let thirteen = T::from_f64(13.0).unwrap_or_else(|| T::zero());
        let twelve = T::from_f64(12.0).unwrap_or_else(|| T::zero());
        let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());
        
        let beta0 = thirteen.clone() / twelve.clone() * 
                   (v[0].clone() - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * v[1].clone() + v[2].clone()).powi(2) +
                   T::one() / four.clone() * 
                   (v[0].clone() - four.clone() * v[1].clone() + three.clone() * v[2].clone()).powi(2);
        
        let beta1 = thirteen.clone() / twelve.clone() * 
                   (v[1].clone() - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * v[2].clone() + v[3].clone()).powi(2) +
                   T::one() / four.clone() * 
                   (v[1].clone() - v[3].clone()).powi(2);
        
        let beta2 = thirteen.clone() / twelve.clone() * 
                   (v[2].clone() - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * v[3].clone() + v[4].clone()).powi(2) +
                   T::one() / four.clone() * 
                   (three.clone() * v[2].clone() - four.clone() * v[3].clone() + v[4].clone()).powi(2);
        
        // WENO weights
        let d0 = T::from_f64(0.1).unwrap_or_else(|| T::zero());
        let d1 = T::from_f64(0.6).unwrap_or_else(|| T::zero());
        let d2 = T::from_f64(0.3).unwrap_or_else(|| T::zero());
        
        let alpha0 = d0 / (T::from_f64(WENO_EPSILON).unwrap_or_else(|| T::zero()) + beta0).powi(WENO_POWER);
        let alpha1 = d1 / (T::from_f64(WENO_EPSILON).unwrap_or_else(|| T::zero()) + beta1).powi(WENO_POWER);
        let alpha2 = d2 / (T::from_f64(WENO_EPSILON).unwrap_or_else(|| T::zero()) + beta2).powi(WENO_POWER);
        
        let sum_alpha = alpha0.clone() + alpha1.clone() + alpha2.clone();
        
        let w0 = alpha0 / sum_alpha.clone();
        let w1 = alpha1 / sum_alpha.clone();
        let w2 = alpha2 / sum_alpha;
        
        // WENO reconstruction and derivative
        (w0 * s0 + w1 * s1 + w2 * s2) / grid.dx.clone()
    }
    
    /// WENO5 scheme in y-direction
    fn weno5_y(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        if j < 2 || j >= grid.data.ncols() - 2 {
            // Fall back to central difference at boundaries
            return self.central_difference_y(grid, i, j);
        }
        
        // WENO5 stencil values
        let v = [
            grid.data[(i, j-2)].clone(),
            grid.data[(i, j-1)].clone(),
            grid.data[(i, j)].clone(),
            grid.data[(i, j+1)].clone(),
            grid.data[(i, j+2)].clone(),
        ];
        
        // Three sub-stencils for reconstruction
        let s0 = (v[0].clone() * T::from_f64(-1.0).unwrap_or_else(|| T::zero()) + 
                  v[1].clone() * T::from_f64(3.0).unwrap_or_else(|| T::zero()) - 
                  v[2].clone() * T::from_f64(3.0).unwrap_or_else(|| T::zero()) + 
                  v[3].clone()) / T::from_f64(6.0).unwrap_or_else(|| T::zero());
                  
        let s1 = (v[1].clone() * T::from_f64(-1.0).unwrap_or_else(|| T::zero()) + 
                  v[3].clone()) / T::from_f64(2.0).unwrap_or_else(|| T::zero());
                  
        let s2 = (v[2].clone() * T::from_f64(-3.0).unwrap_or_else(|| T::zero()) + 
                  v[3].clone() * T::from_f64(3.0).unwrap_or_else(|| T::zero()) + 
                  v[4].clone() * T::from_f64(-1.0).unwrap_or_else(|| T::zero())) / T::from_f64(6.0).unwrap_or_else(|| T::zero());
        
        // Compute smoothness indicators
        let thirteen = T::from_f64(13.0).unwrap_or_else(|| T::zero());
        let twelve = T::from_f64(12.0).unwrap_or_else(|| T::zero());
        let three = T::from_f64(3.0).unwrap_or_else(|| T::zero());
        let four = T::from_f64(4.0).unwrap_or_else(|| T::zero());
        
        let beta0 = thirteen.clone() / twelve.clone() * 
                   (v[0].clone() - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * v[1].clone() + v[2].clone()).powi(2) +
                   T::one() / four.clone() * 
                   (v[0].clone() - four.clone() * v[1].clone() + three.clone() * v[2].clone()).powi(2);
        
        let beta1 = thirteen.clone() / twelve.clone() * 
                   (v[1].clone() - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * v[2].clone() + v[3].clone()).powi(2) +
                   T::one() / four.clone() * 
                   (v[1].clone() - v[3].clone()).powi(2);
        
        let beta2 = thirteen.clone() / twelve.clone() * 
                   (v[2].clone() - T::from_f64(2.0).unwrap_or_else(|| T::zero()) * v[3].clone() + v[4].clone()).powi(2) +
                   T::one() / four.clone() * 
                   (three.clone() * v[2].clone() - four.clone() * v[3].clone() + v[4].clone()).powi(2);
        
        // WENO weights
        let d0 = T::from_f64(0.1).unwrap_or_else(|| T::zero());
        let d1 = T::from_f64(0.6).unwrap_or_else(|| T::zero());
        let d2 = T::from_f64(0.3).unwrap_or_else(|| T::zero());
        
        let alpha0 = d0 / (T::from_f64(WENO_EPSILON).unwrap_or_else(|| T::zero()) + beta0).powi(WENO_POWER);
        let alpha1 = d1 / (T::from_f64(WENO_EPSILON).unwrap_or_else(|| T::zero()) + beta1).powi(WENO_POWER);
        let alpha2 = d2 / (T::from_f64(WENO_EPSILON).unwrap_or_else(|| T::zero()) + beta2).powi(WENO_POWER);
        
        let sum_alpha = alpha0.clone() + alpha1.clone() + alpha2.clone();
        
        let w0 = alpha0 / sum_alpha.clone();
        let w1 = alpha1 / sum_alpha.clone();
        let w2 = alpha2 / sum_alpha;
        
        // WENO reconstruction and derivative
        (w0 * s0 + w1 * s1 + w2 * s2) / grid.dy.clone()
    }
    
    /// Fourth-order explicit central difference in x-direction
    fn fourth_order_central_x(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        // Fourth-order compact scheme
        let twelve = T::from_f64(12.0).unwrap_or_else(|| T::zero());
        
        (grid.data[(i-2, j)].clone() * T::from_f64(-1.0).unwrap_or_else(|| T::zero()) +
         grid.data[(i-1, j)].clone() * T::from_f64(8.0).unwrap_or_else(|| T::zero()) -
         grid.data[(i+1, j)].clone() * T::from_f64(8.0).unwrap_or_else(|| T::zero()) +
         grid.data[(i+2, j)].clone()) / (twelve * grid.dx.clone())
    }
    
    /// Fourth-order explicit central difference in y-direction
    fn fourth_order_central_y(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let twelve = T::from_f64(12.0).unwrap_or_else(|| T::zero());
        let eight = T::from_f64(8.0).unwrap_or_else(|| T::zero());
        
        (-grid.data[(i, j+2)].clone() + 
         eight.clone() * grid.data[(i, j+1)].clone() - 
         eight * grid.data[(i, j-1)].clone() + 
         grid.data[(i, j-2)].clone()) / (twelve * grid.dy.clone())
    }
    
    /// Compute gradient ratio for flux limiting (x-direction)
    fn compute_gradient_ratio_x(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let num = grid.data[(i, j)].clone() - grid.data[(i-1, j)].clone();
        let den = grid.data[(i+1, j)].clone() - grid.data[(i, j)].clone();
        
        if den.clone().abs() < T::from_f64(MIN_LIMITER_THRESHOLD).unwrap_or_else(|| T::zero()) {
            T::zero()
        } else {
            num / den
        }
    }
    
    /// Compute gradient ratio for flux limiting (y-direction)
    fn compute_gradient_ratio_y(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let num = grid.data[(i, j)].clone() - grid.data[(i, j-1)].clone();
        let den = grid.data[(i, j+1)].clone() - grid.data[(i, j)].clone();
        
        if den.clone().abs() < T::from_f64(MIN_LIMITER_THRESHOLD).unwrap_or_else(|| T::zero()) {
            T::zero()
        } else {
            num / den
        }
    }

    /// Compute gradient ratio for flux limiting (x-direction) for negative velocity
    fn compute_gradient_ratio_x_negative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let num = grid.data[(i, j)].clone() - grid.data[(i+1, j)].clone();
        let den = grid.data[(i-1, j)].clone() - grid.data[(i, j)].clone();
        
        if den.clone().abs() < T::from_f64(MIN_LIMITER_THRESHOLD).unwrap_or_else(|| T::zero()) {
            T::zero()
        } else {
            num / den
        }
    }

    /// Compute gradient ratio for flux limiting (y-direction) for negative velocity
    fn compute_gradient_ratio_y_negative(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        let num = grid.data[(i, j)].clone() - grid.data[(i, j+1)].clone();
        let den = grid.data[(i, j-1)].clone() - grid.data[(i, j)].clone();
        
        if den.clone().abs() < T::from_f64(MIN_LIMITER_THRESHOLD).unwrap_or_else(|| T::zero()) {
            T::zero()
        } else {
            num / den
        }
    }
    
    /// Apply flux limiter
    fn apply_limiter(&self, r: T) -> T {
        match self.limiter {
            Some(FluxLimiter::None) => T::one(),
            Some(FluxLimiter::Moderate) => self.moderate_limiter(r),
            Some(FluxLimiter::VanLeer) => self.van_leer_limiter(r),
            Some(FluxLimiter::VanAlbada) => self.van_albada_limiter(r),
            Some(FluxLimiter::Superbee) => self.superbee_limiter(r),
            Some(FluxLimiter::MC) => self.mc_limiter(r),
            None => T::one(), // No limiter
        }
    }
    
    /// Moderate limiter (prevents oscillations)
    fn moderate_limiter(&self, r: T) -> T {
        T::zero().max(T::one().min(r))
    }
    
    /// Van Leer limiter
    fn van_leer_limiter(&self, r: T) -> T {
        let abs_r = r.clone().abs();
        (r.clone() + abs_r.clone()) / (T::one() + abs_r)
    }
    
    /// Van Albada limiter
    fn van_albada_limiter(&self, r: T) -> T {
        (r.clone() * r.clone() + r.clone()) / (r.clone() * r + T::one())
    }
    
    /// Superbee limiter
    fn superbee_limiter(&self, r: T) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        T::zero().max(T::one().min(two.clone() * r.clone()))
            .max(two.clone().min(r))
    }
    
    /// MC (Monotonized Central) limiter
    fn mc_limiter(&self, r: T) -> T {
        let two = T::from_f64(2.0).unwrap_or_else(|| T::zero());
        let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
        
        T::zero().max(
            ((T::one() + r.clone()) * half).min(two.clone())
                .min(two * r)
        )
    }
    
    /// Compute second derivative in x-direction
    pub fn d2dx2(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        (grid.data[(i+1, j)].clone() - 
         T::from_f64(2.0).unwrap_or_else(|| T::zero()) * grid.data[(i, j)].clone() + 
         grid.data[(i-1, j)].clone()) / grid.dx.clone().powi(2)
    }
    
    /// Compute second derivative in y-direction
    pub fn d2dy2(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        (grid.data[(i, j+1)].clone() - 
         T::from_f64(2.0).unwrap_or_else(|| T::zero()) * grid.data[(i, j)].clone() + 
         grid.data[(i, j-1)].clone()) / grid.dy.clone().powi(2)
    }
    
    /// Compute Laplacian
    pub fn laplacian(&self, grid: &Grid2D<T>, i: usize, j: usize) -> T {
        self.d2dx2(grid, i, j) + self.d2dy2(grid, i, j)
    }
}

/// Time integration for 2D problems
pub struct TimeIntegrator<T: RealField> {
    scheme: TimeScheme,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + FromPrimitive> TimeIntegrator<T> {
    /// Create a new time integrator
    pub fn new(scheme: TimeScheme) -> Self {
        Self {
            scheme,
            _phantom: std::marker::PhantomData,
        }
    }
    
    /// Advance solution in time
    pub fn advance<F>(
        &self,
        grid: &mut Grid2D<T>,
        dt: T,
        rhs: F,
    ) -> Result<(), Error>
    where
        F: Fn(&Grid2D<T>) -> Grid2D<T>,
    {
        match self.scheme {
            TimeScheme::ForwardEuler => self.forward_euler(grid, dt, rhs),
            TimeScheme::RungeKutta2 => self.runge_kutta2(grid, dt, rhs),
            TimeScheme::RungeKutta4 => self.runge_kutta4(grid, dt, rhs),
            _ => Err(Error::NotImplemented("Time scheme not yet implemented".to_string())),
        }
    }
    
    /// Forward Euler time integration
    fn forward_euler<F>(
        &self,
        grid: &mut Grid2D<T>,
        dt: T,
        rhs: F,
    ) -> Result<(), Error>
    where
        F: Fn(&Grid2D<T>) -> Grid2D<T>,
    {
        let k1 = rhs(grid);
        grid.data += k1.data * dt;
        Ok(())
    }
    
    /// Second-order Runge-Kutta
    fn runge_kutta2<F>(
        &self,
        grid: &mut Grid2D<T>,
        dt: T,
        rhs: F,
    ) -> Result<(), Error>
    where
        F: Fn(&Grid2D<T>) -> Grid2D<T>,
    {
        let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
        
        // Stage 1
        let k1 = rhs(grid);
        let mut temp = grid.clone();
        temp.data += k1.data * dt.clone() * half.clone();
        
        // Stage 2
        let k2 = rhs(&temp);
        grid.data += k2.data * dt;
        
        Ok(())
    }
    
    /// Fourth-order Runge-Kutta
    fn runge_kutta4<F>(
        &self,
        grid: &mut Grid2D<T>,
        dt: T,
        rhs: F,
    ) -> Result<(), Error>
    where
        F: Fn(&Grid2D<T>) -> Grid2D<T>,
    {
        let half = T::from_f64(0.5).unwrap_or_else(|| T::zero());
        let sixth = T::from_f64(1.0/6.0).unwrap_or_else(|| T::zero());
        let third = T::from_f64(1.0/3.0).unwrap_or_else(|| T::zero());
        
        // Stage 1
        let k1 = rhs(grid);
        let mut temp = grid.clone();
        temp.data = grid.data.clone() + k1.data.clone() * dt.clone() * half.clone();
        
        // Stage 2
        let k2 = rhs(&temp);
        temp.data = grid.data.clone() + k2.data.clone() * dt.clone() * half;
        
        // Stage 3
        let k3 = rhs(&temp);
        temp.data = grid.data.clone() + k3.data.clone() * dt.clone();
        
        // Stage 4
        let k4 = rhs(&temp);
        
        // Combine stages
        grid.data += (k1.data + k2.data * T::from_f64(2.0).unwrap_or_else(|| T::zero()) + 
                      k3.data * T::from_f64(2.0).unwrap_or_else(|| T::zero()) + k4.data) * dt * sixth;
        
        Ok(())
    }
    
    /// Compute stable time step based on CFL condition
    pub fn compute_stable_dt(
        &self,
        grid: &Grid2D<T>,
        max_velocity: T,
        cfl: T,
    ) -> T {
        let dx_term = grid.dx.clone() / max_velocity.clone();
        let dy_term = grid.dy.clone() / max_velocity;
        
        cfl * dx_term.min(dy_term)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_grid_creation() {
        let grid = Grid2D::<f64>::new(10, 10, 0.1, 0.1, 2);
        assert_eq!(grid.data.shape(), (14, 14));
        assert_eq!(grid.interior_shape(), (10, 10));
    }
    
    #[test]
    fn test_central_difference() {
        let mut grid = Grid2D::<f64>::new(5, 5, 1.0, 1.0, 2);
        
        // Set up a linear function u = x
        for i in 0..9 {
            for j in 0..9 {
                grid.data[(i, j)] = i as f64;
            }
        }
        
        let fd = FiniteDifference::new(SpatialScheme::CentralDifference);
        let ddx = fd.ddx(&grid, 4, 4, None);
        
        // For u = x, du/dx = 1
        assert!((ddx - 1.0).abs() < 1e-10);
    }
    
    #[test]
    fn test_time_integration() {
        let mut grid = Grid2D::<f64>::new(5, 5, 1.0, 1.0, 2);
        grid.data.fill(1.0);
        
        let integrator = TimeIntegrator::new(TimeScheme::ForwardEuler);
        
        // Standard decay: du/dt = -u
        let rhs = |g: &Grid2D<f64>| {
            let mut result = g.clone();
            result.data = -g.data.clone();
            result
        };
        
        integrator.advance(&mut grid, 0.1, rhs).expect("CRITICAL: Add proper error handling");
        
        // After one step: u = 1.0 * (1 - 0.1) = 0.9
        assert!((grid.data[(4, 4)] - 0.9).abs() < 1e-10);
    }
}
