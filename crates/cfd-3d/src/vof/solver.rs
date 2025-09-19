//! Core VOF solver implementation

use cfd_core::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

use super::advection::AdvectionMethod;
use super::config::VofConfig;
use super::initialization::Initialization;
use super::reconstruction::InterfaceReconstruction;

/// VOF solver for multiphase flow
pub struct VofSolver<T: RealField + FromPrimitive + Copy> {
    pub(crate) config: VofConfig,
    /// Grid dimensions
    pub(crate) nx: usize,
    pub(crate) ny: usize,
    pub(crate) nz: usize,
    /// Grid spacing
    pub(crate) dx: T,
    pub(crate) dy: T,
    pub(crate) dz: T,
    /// Volume fraction field (0 = phase 1, 1 = phase 2)
    pub(crate) alpha: Vec<T>,
    /// Previous timestep volume fraction
    pub(crate) alpha_previous: Vec<T>,
    /// Velocity field
    pub(crate) velocity: Vec<Vector3<T>>,
    /// Interface normal vectors
    pub(crate) normals: Vec<Vector3<T>>,
    /// Interface curvature
    pub(crate) curvature: Vec<T>,
}

impl<T: RealField + FromPrimitive + Copy> VofSolver<T> {
    /// Create VOF solver instance
    pub fn create(config: VofConfig, nx: usize, ny: usize, nz: usize, dx: T, dy: T, dz: T) -> Self {
        let grid_size = nx * ny * nz;
        Self {
            config,
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            alpha: vec![T::zero(); grid_size],
            alpha_previous: vec![T::zero(); grid_size],
            velocity: vec![Vector3::zeros(); grid_size],
            normals: vec![Vector3::zeros(); grid_size],
            curvature: vec![T::zero(); grid_size],
        }
    }

    /// Convert 3D indices to linear index
    pub(crate) fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.ny * self.nx + j * self.nx + i
    }

    /// Set velocity field
    pub fn set_velocity_field(&mut self, velocity: Vec<Vector3<T>>) -> Result<()> {
        if velocity.len() != self.velocity.len() {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: self.velocity.len(),
                actual: velocity.len(),
            });
        }
        self.velocity = velocity;
        Ok(())
    }

    /// Get volume fraction field
    pub fn get_volume_fraction(&self) -> &[T] {
        &self.alpha
    }

    /// Get interface normals
    pub fn get_normals(&self) -> &[Vector3<T>] {
        &self.normals
    }

    /// Get interface curvature
    pub fn get_curvature(&self) -> &[T] {
        &self.curvature
    }

    /// Main time step
    pub fn advance(&mut self, dt: T) -> Result<()> {
        // Store previous values
        self.alpha_previous.copy_from_slice(&self.alpha);

        // Reconstruct interface
        self.reconstruct_interface();

        // Advect volume fraction
        let advection = AdvectionMethod::create(&self.config);
        advection.advect(self, dt)?;

        // Apply compression if enabled
        if self.config.enable_compression {
            advection.apply_compression(self, dt)?;
        }

        Ok(())
    }

    /// Reconstruct interface normals and curvature
    pub fn reconstruct_interface(&mut self) {
        let reconstruction = InterfaceReconstruction::create(&self.config);
        reconstruction.reconstruct(self);
    }

    /// Initialize the volume fraction field
    pub fn initialize(&mut self, init: Initialization<T>) {
        init.apply(self);
        self.reconstruct_interface();
    }

    /// Calculate time step based on CFL condition
    pub fn calculate_timestep(&self) -> T {
        let mut max_velocity = T::zero();
        for vel in &self.velocity {
            let speed = vel.norm();
            if speed > max_velocity {
                max_velocity = speed;
            }
        }

        if max_velocity > T::zero() {
            let dx_min = self.dx.min(self.dy).min(self.dz);
            T::from_f64(self.config.cfl_number).unwrap_or(T::from_f64(0.3).unwrap_or(T::one()))
                * dx_min
                / max_velocity
        } else {
            T::from_f64(1e-3).unwrap_or(T::one())
        }
    }

    /// Check mass conservation
    pub fn total_volume(&self) -> T {
        let mut total = T::zero();
        let cell_volume = self.dx * self.dy * self.dz;
        for &alpha in &self.alpha {
            total += alpha * cell_volume;
        }
        total
    }
}
