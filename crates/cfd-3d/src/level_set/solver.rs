//! Level Set solver implementation

use super::config::LevelSetConfig;
use cfd_core::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Level Set solver for interface tracking
pub struct LevelSetSolver<T: RealField + FromPrimitive + Copy> {
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
    /// Previous timestep level set function
    phi_previous: Vec<T>,
    /// Velocity field
    velocity: Vec<Vector3<T>>,
    /// Narrow band indices (if using narrow band)
    narrow_band: Vec<usize>,
    /// Time step counter
    time_step: usize,
}

impl<T: RealField + FromPrimitive + Copy> LevelSetSolver<T> {
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
            phi_previous: vec![T::zero(); grid_size],
            velocity: vec![Vector3::zeros(); grid_size],
            narrow_band: Vec::new(),
            time_step: 0,
        }
    }

    /// Get grid index from (i, j, k) coordinates
    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.ny * self.nx + j * self.nx + i
    }

    /// Get the level set function
    pub fn phi(&self) -> &[T] {
        &self.phi
    }

    /// Get mutable reference to level set function
    pub fn phi_mut(&mut self) -> &mut [T] {
        &mut self.phi
    }

    /// Set velocity field
    pub fn set_velocity(&mut self, velocity: Vec<Vector3<T>>) {
        self.velocity = velocity;
    }

    /// Get narrow band indices
    pub fn narrow_band(&self) -> &[usize] {
        &self.narrow_band
    }

    /// Update narrow band indices based on current level set
    pub fn update_narrow_band(&mut self) {
        self.narrow_band.clear();
        let band_width = T::from_f64(self.config.band_width).unwrap_or_else(T::zero);

        for idx in 0..self.phi.len() {
            if self.phi[idx].abs() <= band_width * self.dx.min(self.dy).min(self.dz) {
                self.narrow_band.push(idx);
            }
        }
    }

    /// Advance level set by one time step
    pub fn advance(&mut self, dt: T) -> Result<()> {
        // Store previous level set
        self.phi_previous.copy_from_slice(&self.phi);

        // Advection step using upwind scheme
        self.advect_level_set(dt)?;

        // Check for reinitialization
        self.time_step += 1;
        if self.time_step % self.config.reinitialization_interval == 0 {
            self.reinitialize()?;
        }

        if self.config.use_narrow_band {
            self.update_narrow_band();
        }

        Ok(())
    }

    /// Advect level set using upwind finite difference
    fn advect_level_set(&mut self, dt: T) -> Result<()> {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;
        let dx = self.dx;
        let dy = self.dy;
        let dz = self.dz;

        // Create temporary storage for new values
        let mut phi_new = self.phi.clone();

        for k in 1..nz - 1 {
            for j in 1..ny - 1 {
                for i in 1..nx - 1 {
                    let idx = self.index(i, j, k);

                    // Get velocity at this point
                    let velocity = self.velocity[idx];
                    let u = velocity.x;
                    let v = velocity.y;
                    let w = velocity.z;

                    // Upwind differencing based on velocity direction
                    let dphi_dx = if u > T::zero() {
                        (self.phi[idx] - self.phi[self.index(i - 1, j, k)]) / dx
                    } else {
                        (self.phi[self.index(i + 1, j, k)] - self.phi[idx]) / dx
                    };

                    let dphi_dy = if v > T::zero() {
                        (self.phi[idx] - self.phi[self.index(i, j - 1, k)]) / dy
                    } else {
                        (self.phi[self.index(i, j + 1, k)] - self.phi[idx]) / dy
                    };

                    let dphi_dz = if w > T::zero() {
                        (self.phi[idx] - self.phi[self.index(i, j, k - 1)]) / dz
                    } else {
                        (self.phi[self.index(i, j, k + 1)] - self.phi[idx]) / dz
                    };

                    // Level set advection equation: ∂φ/∂t + u·∇φ = 0
                    phi_new[idx] = self.phi[idx] - dt * (u * dphi_dx + v * dphi_dy + w * dphi_dz);
                }
            }
        }

        // Update level set
        self.phi = phi_new;
        Ok(())
    }

    /// Reinitialize level set to signed distance function
    fn reinitialize(&mut self) -> Result<()> {
        let iterations = 10; // Number of reinitialization iterations
        let dtau = T::from_f64(0.5).unwrap_or_else(|| T::one()) * self.dx.min(self.dy).min(self.dz);

        for _ in 0..iterations {
            let mut phi_new = self.phi.clone();

            for k in 1..self.nz - 1 {
                for j in 1..self.ny - 1 {
                    for i in 1..self.nx - 1 {
                        let idx = self.index(i, j, k);

                        // Calculate gradient using WENO or central differences
                        let grad = self.calculate_gradient_magnitude(i, j, k);

                        // Sign function
                        let sign_phi = self.phi[idx]
                            / (self.phi[idx].abs()
                                + T::from_f64(1e-6).unwrap_or_else(|| T::zero()));

                        // Reinitialization equation: ∂φ/∂τ + S(φ₀)(|∇φ| - 1) = 0
                        phi_new[idx] = self.phi[idx] - dtau * sign_phi * (grad - T::one());
                    }
                }
            }

            self.phi = phi_new;
        }

        Ok(())
    }

    /// Calculate gradient magnitude at a grid point
    fn calculate_gradient_magnitude(&self, i: usize, j: usize, k: usize) -> T {
        let dx = self.dx;
        let dy = self.dy;
        let dz = self.dz;

        // Central differences for gradient
        let dphi_dx = (self.phi[self.index(i + 1, j, k)] - self.phi[self.index(i - 1, j, k)])
            / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dx);
        let dphi_dy = (self.phi[self.index(i, j + 1, k)] - self.phi[self.index(i, j - 1, k)])
            / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dy);
        let dphi_dz = (self.phi[self.index(i, j, k + 1)] - self.phi[self.index(i, j, k - 1)])
            / (T::from_f64(2.0).unwrap_or_else(|| T::one()) * dz);

        (dphi_dx * dphi_dx + dphi_dy * dphi_dy + dphi_dz * dphi_dz).sqrt()
    }
}
