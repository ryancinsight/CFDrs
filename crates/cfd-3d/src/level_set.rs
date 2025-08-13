//! Level Set Method for interface tracking in 3D multiphase flows
//!
//! The Level Set method represents interfaces as the zero level set of a
//! signed distance function, providing accurate interface tracking.

use cfd_core::Result;
use nalgebra::{Vector3, RealField};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

// Named constants for Level Set
const DEFAULT_REINITIALIZATION_INTERVAL: usize = 5;
const DEFAULT_BAND_WIDTH: f64 = 5.0;  // Width of narrow band in grid cells
const DEFAULT_CFL_NUMBER: f64 = 0.5;
const DEFAULT_TOLERANCE: f64 = 1e-6;
const DEFAULT_MAX_ITERATIONS: usize = 100;
const EPSILON_SMOOTHING: f64 = 1.5;  // Interface thickness for smoothing
const WENO_ORDER: usize = 5;  // WENO scheme order for advection

/// Level Set configuration
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LevelSetConfig {
    /// Reinitialization interval (timesteps)
    pub reinitialization_interval: usize,
    /// Narrow band width
    pub band_width: f64,
    /// CFL number for time stepping
    pub cfl_number: f64,
    /// Convergence tolerance
    pub tolerance: f64,
    /// Maximum iterations for reinitialization
    pub max_iterations: usize,
    /// Use narrow band method
    pub use_narrow_band: bool,
    /// Use WENO scheme for advection
    pub use_weno: bool,
}

impl Default for LevelSetConfig {
    fn default() -> Self {
        Self {
            reinitialization_interval: DEFAULT_REINITIALIZATION_INTERVAL,
            band_width: DEFAULT_BAND_WIDTH,
            cfl_number: DEFAULT_CFL_NUMBER,
            tolerance: DEFAULT_TOLERANCE,
            max_iterations: DEFAULT_MAX_ITERATIONS,
            use_narrow_band: false,
            use_weno: true,
        }
    }
}

/// Level Set solver for interface tracking
pub struct LevelSetSolver<T: RealField + FromPrimitive> {
    config: LevelSetConfig,
    /// Grid dimensions
    nx: usize,
    ny: usize,
    nz: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    dz: T,
    /// Level set function (signed distance)
    phi: Vec<T>,
    /// Previous level set function
    phi_old: Vec<T>,
    /// Velocity field
    velocity: Vec<Vector3<T>>,
    /// Narrow band indices (if using narrow band)
    narrow_band: Vec<usize>,
    /// Time step counter
    time_step: usize,
}

impl<T: RealField + FromPrimitive> LevelSetSolver<T> {
    /// Create a new Level Set solver
    pub fn new(
        config: LevelSetConfig,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: T,
        dy: T,
        dz: T,
    ) -> Self {
        let grid_size = nx * ny * nz;
        Self {
            config,
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            phi: vec![T::zero(); grid_size],
            phi_old: vec![T::zero(); grid_size],
            velocity: vec![Vector3::zeros(); grid_size],
            narrow_band: Vec::new(),
            time_step: 0,
        }
    }
    
    /// Initialize level set with a sphere
    pub fn initialize_sphere(&mut self, center: Vector3<T>, radius: T) {
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let x = T::from_usize(i).unwrap() * self.dx.clone();
                    let y = T::from_usize(j).unwrap() * self.dy.clone();
                    let z = T::from_usize(k).unwrap() * self.dz.clone();
                    
                    let pos = Vector3::new(x, y, z);
                    let distance = (pos - center.clone()).norm();
                    
                    let idx = self.index(i, j, k);
                    self.phi[idx] = distance - radius.clone();
                }
            }
        }
        
        if self.config.use_narrow_band {
            self.update_narrow_band();
        }
    }
    
    /// Initialize level set with a box
    pub fn initialize_box(&mut self, min_corner: Vector3<T>, max_corner: Vector3<T>) {
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let x = T::from_usize(i).unwrap() * self.dx.clone();
                    let y = T::from_usize(j).unwrap() * self.dy.clone();
                    let z = T::from_usize(k).unwrap() * self.dz.clone();
                    
                    // Signed distance to box
                    let dx_min = x.clone() - min_corner[0].clone();
                    let dx_max = max_corner[0].clone() - x.clone();
                    let dy_min = y.clone() - min_corner[1].clone();
                    let dy_max = max_corner[1].clone() - y.clone();
                    let dz_min = z.clone() - min_corner[2].clone();
                    let dz_max = max_corner[2].clone() - z;
                    
                    let dx_dist = dx_min.min(dx_max);
                    let dy_dist = dy_min.min(dy_max);
                    let dz_dist = dz_min.min(dz_max);
                    
                    let idx = self.index(i, j, k);
                    self.phi[idx] = -dx_dist.min(dy_dist).min(dz_dist);
                }
            }
        }
        
        if self.config.use_narrow_band {
            self.update_narrow_band();
        }
    }
    
    /// Convert 3D indices to linear index
    fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.nx * self.ny + j * self.nx + i
    }
    
    /// Update narrow band around interface
    fn update_narrow_band(&mut self) {
        self.narrow_band.clear();
        let band_width = T::from_f64(self.config.band_width).unwrap() 
            * self.dx.clone().max(self.dy.clone()).max(self.dz.clone());
        
        for idx in 0..self.phi.len() {
            if self.phi[idx].clone().abs() <= band_width {
                self.narrow_band.push(idx);
            }
        }
    }
    
    /// WENO5 scheme for spatial derivatives (Jiang & Peng, 2000)
    fn weno5_derivative(&self, v: &[T], direction: usize, i: usize, j: usize, k: usize) -> (T, T) {
        let idx = self.index(i, j, k);
        let eps = T::from_f64(1e-6).unwrap();
        
        // Get stencil values based on direction
        let (vm2, vm1, v0, vp1, vp2, _vp3) = match direction {
            0 => {  // x-direction
                if i < 2 || i >= self.nx - 3 {
                    // Use lower order at boundaries
                    let vm1 = if i > 0 { v[self.index(i-1, j, k)].clone() } else { v[idx].clone() };
                    let vp1 = if i < self.nx-1 { v[self.index(i+1, j, k)].clone() } else { v[idx].clone() };
                    return ((v[idx].clone() - vm1) / self.dx.clone(),
                            (vp1 - v[idx].clone()) / self.dx.clone());
                }
                (v[self.index(i-2, j, k)].clone(),
                 v[self.index(i-1, j, k)].clone(),
                 v[idx].clone(),
                 v[self.index(i+1, j, k)].clone(),
                 v[self.index(i+2, j, k)].clone(),
                 v[self.index(i.min(self.nx-1).min(i+3), j, k)].clone())
            },
            1 => {  // y-direction
                if j < 2 || j >= self.ny - 3 {
                    let vm1 = if j > 0 { v[self.index(i, j-1, k)].clone() } else { v[idx].clone() };
                    let vp1 = if j < self.ny-1 { v[self.index(i, j+1, k)].clone() } else { v[idx].clone() };
                    return ((v[idx].clone() - vm1) / self.dy.clone(),
                            (vp1 - v[idx].clone()) / self.dy.clone());
                }
                (v[self.index(i, j-2, k)].clone(),
                 v[self.index(i, j-1, k)].clone(),
                 v[idx].clone(),
                 v[self.index(i, j+1, k)].clone(),
                 v[self.index(i, j+2, k)].clone(),
                 v[self.index(i, j.min(self.ny-1).min(j+3), k)].clone())
            },
            _ => {  // z-direction
                if k < 2 || k >= self.nz - 3 {
                    let vm1 = if k > 0 { v[self.index(i, j, k-1)].clone() } else { v[idx].clone() };
                    let vp1 = if k < self.nz-1 { v[self.index(i, j, k+1)].clone() } else { v[idx].clone() };
                    return ((v[idx].clone() - vm1) / self.dz.clone(),
                            (vp1 - v[idx].clone()) / self.dz.clone());
                }
                (v[self.index(i, j, k-2)].clone(),
                 v[self.index(i, j, k-1)].clone(),
                 v[idx].clone(),
                 v[self.index(i, j, k+1)].clone(),
                 v[self.index(i, j, k+2)].clone(),
                 v[self.index(i, j, k.min(self.nz-1).min(k+3))].clone())
            }
        };
        
        // WENO5 reconstruction for negative and positive derivatives
        let two = T::from_f64(2.0).unwrap();
        let three = T::from_f64(3.0).unwrap();
        let six = T::from_f64(6.0).unwrap();
        let thirteen = T::from_f64(13.0).unwrap();
        
        // Smoothness indicators
        let beta0 = (thirteen.clone() * (vm2.clone() - two.clone() * vm1.clone() + v0.clone()).powi(2)
            + three.clone() * (vm2.clone() - T::from_f64(4.0).unwrap() * vm1.clone() + three.clone() * v0.clone()).powi(2)) / T::from_f64(12.0).unwrap();
        let beta1 = (thirteen.clone() * (vm1.clone() - two.clone() * v0.clone() + vp1.clone()).powi(2)
            + three.clone() * (vm1.clone() - vp1.clone()).powi(2)) / T::from_f64(12.0).unwrap();
        let beta2 = (thirteen.clone() * (v0.clone() - two.clone() * vp1.clone() + vp2.clone()).powi(2)
            + three.clone() * (three.clone() * v0.clone() - T::from_f64(4.0).unwrap() * vp1.clone() + vp2.clone()).powi(2)) / T::from_f64(12.0).unwrap();
        
        // Weights
        let d0 = T::from_f64(0.1).unwrap();
        let d1 = T::from_f64(0.6).unwrap();
        let d2 = T::from_f64(0.3).unwrap();
        
        let alpha0 = d0 / (eps.clone() + beta0).powi(2);
        let alpha1 = d1 / (eps.clone() + beta1).powi(2);
        let alpha2 = d2 / (eps.clone() + beta2).powi(2);
        
        let sum_alpha = alpha0.clone() + alpha1.clone() + alpha2.clone();
        
        let w0 = alpha0 / sum_alpha.clone();
        let w1 = alpha1 / sum_alpha.clone();
        let w2 = alpha2 / sum_alpha;
        
        // Compute derivatives
        let h = match direction {
            0 => self.dx.clone(),
            1 => self.dy.clone(),
            _ => self.dz.clone(),
        };
        
        // Negative-biased stencil (for positive velocities)
        let derivative_minus = (w0.clone() * (two.clone() * vm1.clone() - T::from_f64(7.0).unwrap() * v0.clone() + T::from_f64(11.0).unwrap() * vp1.clone())
            + w1.clone() * (-vm1.clone() + T::from_f64(5.0).unwrap() * v0.clone() + two.clone() * vp1.clone())
            + w2.clone() * (two.clone() * v0.clone() + T::from_f64(5.0).unwrap() * vp1.clone() - vp2.clone())) / (six.clone() * h.clone());
        
        // Positive-biased stencil (for negative velocities)
        // Need to recalculate with shifted stencil
        let derivative_plus = (w2 * (T::from_f64(11.0).unwrap() * v0.clone() - T::from_f64(7.0).unwrap() * vp1.clone() + two.clone() * vp2)
            + w1 * (two.clone() * vm1.clone() + T::from_f64(5.0).unwrap() * v0.clone() - vp1.clone())
            + w0 * (-vm2.clone() + T::from_f64(5.0).unwrap() * vm1 + two * v0)) / (six * h);
        
        (derivative_minus, derivative_plus)
    }
    
    /// Advect level set using velocity field
    pub fn advect(&mut self, dt: T) {
        // Check CFL condition for stability
        let max_velocity = self.velocity.iter()
            .map(|v| v[0].clone().abs().max(v[1].clone().abs()).max(v[2].clone().abs()))
            .fold(T::zero(), |acc, v| acc.max(v));
        
        let min_spacing = self.dx.clone().min(self.dy.clone()).min(self.dz.clone());
        let cfl = max_velocity.clone() * dt.clone() / min_spacing.clone();
        
        if cfl > T::from_f64(self.config.cfl_number).unwrap() {
            // Warning: CFL condition violated, stability may be compromised
            // In production, should adapt time step or use sub-stepping
            let recommended_dt = T::from_f64(self.config.cfl_number).unwrap() * min_spacing / max_velocity;
            eprintln!("Warning: CFL = {:?} > {:?}, recommended dt = {:?}", 
                     cfl, self.config.cfl_number, recommended_dt);
        }
        
        self.phi_old.clone_from(&self.phi);
        
        if self.config.use_weno {
            // WENO5 scheme for advection
            for k in 0..self.nz {
                for j in 0..self.ny {
                    for i in 0..self.nx {
                        let idx = self.index(i, j, k);
                        let vel = &self.velocity[idx];
                        
                        // Compute spatial derivatives using WENO5
                        let (dphi_dx_minus, dphi_dx_plus) = self.weno5_derivative(&self.phi_old, 0, i, j, k);
                        let (dphi_dy_minus, dphi_dy_plus) = self.weno5_derivative(&self.phi_old, 1, i, j, k);
                        let (dphi_dz_minus, dphi_dz_plus) = self.weno5_derivative(&self.phi_old, 2, i, j, k);
                        
                        // Upwind scheme based on velocity direction
                        let dphi_dx = if vel[0] > T::zero() { dphi_dx_minus } else { dphi_dx_plus };
                        let dphi_dy = if vel[1] > T::zero() { dphi_dy_minus } else { dphi_dy_plus };
                        let dphi_dz = if vel[2] > T::zero() { dphi_dz_minus } else { dphi_dz_plus };
                        
                        // Update level set
                        self.phi[idx] = self.phi_old[idx].clone() 
                            - dt.clone() * (vel[0].clone() * dphi_dx 
                                          + vel[1].clone() * dphi_dy 
                                          + vel[2].clone() * dphi_dz);
                    }
                }
            }
        } else {
            // Simple upwind scheme
            for k in 1..self.nz-1 {
                for j in 1..self.ny-1 {
                    for i in 1..self.nx-1 {
                        let idx = self.index(i, j, k);
                        let vel = &self.velocity[idx];
                        
                        // Upwind differences
                        let dphi_dx = if vel[0] > T::zero() {
                            (self.phi_old[idx].clone() - self.phi_old[self.index(i-1, j, k)].clone()) / self.dx.clone()
                        } else {
                            (self.phi_old[self.index(i+1, j, k)].clone() - self.phi_old[idx].clone()) / self.dx.clone()
                        };
                        
                        let dphi_dy = if vel[1] > T::zero() {
                            (self.phi_old[idx].clone() - self.phi_old[self.index(i, j-1, k)].clone()) / self.dy.clone()
                        } else {
                            (self.phi_old[self.index(i, j+1, k)].clone() - self.phi_old[idx].clone()) / self.dy.clone()
                        };
                        
                        let dphi_dz = if vel[2] > T::zero() {
                            (self.phi_old[idx].clone() - self.phi_old[self.index(i, j, k-1)].clone()) / self.dz.clone()
                        } else {
                            (self.phi_old[self.index(i, j, k+1)].clone() - self.phi_old[idx].clone()) / self.dz.clone()
                        };
                        
                        self.phi[idx] = self.phi_old[idx].clone() 
                            - dt.clone() * (vel[0].clone() * dphi_dx 
                                          + vel[1].clone() * dphi_dy 
                                          + vel[2].clone() * dphi_dz);
                    }
                }
            }
        }
        
        self.time_step += 1;
        
        // Reinitialize if needed
        if self.time_step % self.config.reinitialization_interval == 0 {
            self.reinitialize();
        }
    }
    
    /// Smooth Heaviside function for better numerical stability
    fn smooth_heaviside(&self, phi: T, epsilon: T) -> T {
        if phi < -epsilon.clone() {
            T::zero()
        } else if phi > epsilon.clone() {
            T::one()
        } else {
            let half = T::from_f64(0.5).unwrap();
            let pi = T::from_f64(std::f64::consts::PI).unwrap();
            half.clone() * (T::one() + phi.clone() / epsilon.clone() 
                + (pi.clone() * phi / epsilon).sin() / pi)
        }
    }
    
    /// Smooth sign function using smooth Heaviside
    fn smooth_sign(&self, phi: T, epsilon: T) -> T {
        let two = T::from_f64(2.0).unwrap();
        two * self.smooth_heaviside(phi, epsilon) - T::one()
    }
    
    /// Reinitialize to signed distance function
    pub fn reinitialize(&mut self) {
        let mut phi_temp = self.phi.clone();
        let epsilon = T::from_f64(EPSILON_SMOOTHING).unwrap() * self.dx.clone();
        let sign_phi = self.phi.iter()
            .map(|p| self.smooth_sign(p.clone(), epsilon.clone()))
            .collect::<Vec<_>>();
        
        let dtau = T::from_f64(0.5).unwrap() * self.dx.clone().min(self.dy.clone()).min(self.dz.clone());
        
        for _ in 0..self.config.max_iterations {
            let phi_old_temp = phi_temp.clone();
            
            for k in 1..self.nz-1 {
                for j in 1..self.ny-1 {
                    for i in 1..self.nx-1 {
                        let idx = self.index(i, j, k);
                        
                        // Godunov scheme for |∇φ|
                        let dphi_dx_plus = (phi_old_temp[self.index(i+1, j, k)].clone() - phi_old_temp[idx].clone()) / self.dx.clone();
                        let dphi_dx_minus = (phi_old_temp[idx].clone() - phi_old_temp[self.index(i-1, j, k)].clone()) / self.dx.clone();
                        let dphi_dy_plus = (phi_old_temp[self.index(i, j+1, k)].clone() - phi_old_temp[idx].clone()) / self.dy.clone();
                        let dphi_dy_minus = (phi_old_temp[idx].clone() - phi_old_temp[self.index(i, j-1, k)].clone()) / self.dy.clone();
                        let dphi_dz_plus = (phi_old_temp[self.index(i, j, k+1)].clone() - phi_old_temp[idx].clone()) / self.dz.clone();
                        let dphi_dz_minus = (phi_old_temp[idx].clone() - phi_old_temp[self.index(i, j, k-1)].clone()) / self.dz.clone();
                        
                        let a_plus = dphi_dx_plus.clone().max(T::zero());
                        let a_minus = dphi_dx_minus.clone().min(T::zero());
                        let b_plus = dphi_dy_plus.clone().max(T::zero());
                        let b_minus = dphi_dy_minus.clone().min(T::zero());
                        let c_plus = dphi_dz_plus.clone().max(T::zero());
                        let c_minus = dphi_dz_minus.clone().min(T::zero());
                        
                        let grad_phi_norm = if sign_phi[idx] > T::zero() {
                            ComplexField::sqrt(
                                a_minus.clone() * a_minus + a_plus.clone() * a_plus
                                + b_minus.clone() * b_minus + b_plus.clone() * b_plus
                                + c_minus.clone() * c_minus + c_plus.clone() * c_plus
                            )
                        } else {
                            ComplexField::sqrt(
                                a_plus.clone() * a_plus + a_minus.clone() * a_minus
                                + b_plus.clone() * b_plus + b_minus.clone() * b_minus
                                + c_plus.clone() * c_plus + c_minus.clone() * c_minus
                            )
                        };
                        
                        phi_temp[idx] = phi_old_temp[idx].clone() 
                            - dtau.clone() * sign_phi[idx].clone() * (grad_phi_norm - T::one());
                    }
                }
            }
            
            // Check convergence
            let error = phi_temp.iter()
                .zip(phi_old_temp.iter())
                .map(|(p1, p2)| (p1.clone() - p2.clone()).abs())
                .fold(T::zero(), |acc, x| acc.max(x));
            
            if error < T::from_f64(self.config.tolerance).unwrap() {
                break;
            }
        }
        
        self.phi = phi_temp;
        
        if self.config.use_narrow_band {
            self.update_narrow_band();
        }
    }
    
    /// Compute curvature of the interface
    pub fn compute_curvature(&self) -> Vec<T> {
        let mut curvature = vec![T::zero(); self.phi.len()];
        
        for k in 1..self.nz-1 {
            for j in 1..self.ny-1 {
                for i in 1..self.nx-1 {
                    let idx = self.index(i, j, k);
                    
                    // Central differences for gradients
                    let phi_x = (self.phi[self.index(i+1, j, k)].clone() - self.phi[self.index(i-1, j, k)].clone()) 
                        / (T::from_f64(2.0).unwrap() * self.dx.clone());
                    let phi_y = (self.phi[self.index(i, j+1, k)].clone() - self.phi[self.index(i, j-1, k)].clone()) 
                        / (T::from_f64(2.0).unwrap() * self.dy.clone());
                    let phi_z = (self.phi[self.index(i, j, k+1)].clone() - self.phi[self.index(i, j, k-1)].clone()) 
                        / (T::from_f64(2.0).unwrap() * self.dz.clone());
                    
                    // Second derivatives
                    let phi_xx = (self.phi[self.index(i+1, j, k)].clone() - T::from_f64(2.0).unwrap() * self.phi[idx].clone() 
                        + self.phi[self.index(i-1, j, k)].clone()) / (self.dx.clone() * self.dx.clone());
                    let phi_yy = (self.phi[self.index(i, j+1, k)].clone() - T::from_f64(2.0).unwrap() * self.phi[idx].clone() 
                        + self.phi[self.index(i, j-1, k)].clone()) / (self.dy.clone() * self.dy.clone());
                    let phi_zz = (self.phi[self.index(i, j, k+1)].clone() - T::from_f64(2.0).unwrap() * self.phi[idx].clone() 
                        + self.phi[self.index(i, j, k-1)].clone()) / (self.dz.clone() * self.dz.clone());
                    
                    // Mixed derivatives
                    let phi_xy = (self.phi[self.index(i+1, j+1, k)].clone() - self.phi[self.index(i+1, j-1, k)].clone()
                        - self.phi[self.index(i-1, j+1, k)].clone() + self.phi[self.index(i-1, j-1, k)].clone())
                        / (T::from_f64(4.0).unwrap() * self.dx.clone() * self.dy.clone());
                    let phi_xz = (self.phi[self.index(i+1, j, k+1)].clone() - self.phi[self.index(i+1, j, k-1)].clone()
                        - self.phi[self.index(i-1, j, k+1)].clone() + self.phi[self.index(i-1, j, k-1)].clone())
                        / (T::from_f64(4.0).unwrap() * self.dx.clone() * self.dz.clone());
                    let phi_yz = (self.phi[self.index(i, j+1, k+1)].clone() - self.phi[self.index(i, j+1, k-1)].clone()
                        - self.phi[self.index(i, j-1, k+1)].clone() + self.phi[self.index(i, j-1, k-1)].clone())
                        / (T::from_f64(4.0).unwrap() * self.dy.clone() * self.dz.clone());
                    
                    // Compute curvature
                    let grad_norm_sq = phi_x.clone() * phi_x.clone() + phi_y.clone() * phi_y.clone() + phi_z.clone() * phi_z.clone();
                    let grad_norm = ComplexField::sqrt(grad_norm_sq.clone() + T::from_f64(1e-10).unwrap());
                    
                    curvature[idx] = (phi_xx.clone() * (phi_y.clone() * phi_y.clone() + phi_z.clone() * phi_z.clone())
                        + phi_yy.clone() * (phi_x.clone() * phi_x.clone() + phi_z.clone() * phi_z.clone())
                        + phi_zz * (phi_x.clone() * phi_x.clone() + phi_y.clone() * phi_y.clone())
                        - T::from_f64(2.0).unwrap() * (phi_x.clone() * phi_y.clone() * phi_xy 
                            + phi_x.clone() * phi_z.clone() * phi_xz 
                            + phi_y * phi_z * phi_yz))
                        / (grad_norm.clone() * grad_norm_sq);
                }
            }
        }
        
        curvature
    }
    
    /// Get volume fraction in each cell
    pub fn get_volume_fraction(&self) -> Vec<T> {
        self.phi.iter()
            .map(|p| {
                let p = p.clone();
                let eps = T::from_f64(EPSILON_SMOOTHING).unwrap() * self.dx.clone();
                if p.clone() > eps.clone() {
                    T::zero()
                } else if p.clone() < -eps.clone() {
                    T::one()
                } else {
                    T::from_f64(0.5).unwrap() * (T::one() - p.clone() / eps.clone() 
                        - ComplexField::sin(T::from_f64(std::f64::consts::PI).unwrap() * p / eps) 
                        / T::from_f64(std::f64::consts::PI).unwrap())
                }
            })
            .collect()
    }
    
    /// Update velocity field
    pub fn set_velocity(&mut self, velocity: Vec<Vector3<T>>) {
        assert_eq!(velocity.len(), self.velocity.len());
        self.velocity = velocity;
    }
    
    /// Get the level set function
    pub fn phi(&self) -> &[T] {
        &self.phi
    }
    
    /// Main time step
    pub fn step(&mut self, dt: T) -> Result<()> {
        // Compute CFL-limited time step
        let max_vel = self.velocity.iter()
            .map(|v| v[0].clone().abs().max(v[1].clone().abs()).max(v[2].clone().abs()))
            .fold(T::zero(), |acc, v| acc.max(v));
        
        let dt_cfl = T::from_f64(self.config.cfl_number).unwrap() 
            * self.dx.clone().min(self.dy.clone()).min(self.dz.clone()) / (max_vel + T::from_f64(1e-10).unwrap());
        
        let dt_actual = dt.min(dt_cfl);
        
        // Advect level set
        self.advect(dt_actual);
        
        Ok(())
    }
}

use nalgebra::ComplexField;

#[cfg(test)]
mod tests {
    use super::*;
    
    #[test]
    fn test_level_set_creation() {
        let config = LevelSetConfig::default();
        let solver: LevelSetSolver<f64> = LevelSetSolver::new(
            config,
            10, 10, 10,
            0.1, 0.1, 0.1,
        );
        
        assert_eq!(solver.phi.len(), 1000);
    }
    
    #[test]
    fn test_sphere_initialization() {
        let config = LevelSetConfig::default();
        let mut solver: LevelSetSolver<f64> = LevelSetSolver::new(
            config,
            20, 20, 20,
            0.1, 0.1, 0.1,
        );
        
        solver.initialize_sphere(Vector3::new(1.0, 1.0, 1.0), 0.5);
        
        // Check that center is inside (negative)
        let center_idx = solver.index(10, 10, 10);
        assert!(solver.phi[center_idx] < 0.0);
        
        // Check that corners are outside (positive)
        let corner_idx = solver.index(0, 0, 0);
        assert!(solver.phi[corner_idx] > 0.0);
    }
    
    #[test]
    fn test_volume_fraction() {
        let config = LevelSetConfig::default();
        let mut solver: LevelSetSolver<f64> = LevelSetSolver::new(
            config,
            10, 10, 10,
            0.1, 0.1, 0.1,
        );
        
        solver.initialize_sphere(Vector3::new(0.5, 0.5, 0.5), 0.3);
        let vf = solver.get_volume_fraction();
        
        // Volume fraction should be between 0 and 1
        for f in &vf {
            assert!(*f >= 0.0 && *f <= 1.0);
        }
    }
}