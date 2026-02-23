//! # Volume of Fluid (VOF) Method for 3D Multiphase Flows
//!
//! This module implements the Volume of Fluid method for accurate tracking
//! of interfaces in multiphase flows with guaranteed volume conservation.
//!
//! ## Mathematical Foundation
//!
//! ### VOF Equation
//! The volume fraction α satisfies:
//!
//! ```math
//! ∂α/∂t + ∇·(α u) = 0
//! ```
//!
//! where α = 1 in fluid 1, α = 0 in fluid 2, and 0 < α < 1 at interfaces.
//!
//! ### Interface Reconstruction
//! **PLIC (Piecewise Linear Interface Construction)**:
//! - Reconstructs linear interface within each cell
//! - Ensures exact volume conservation
//! - Maintains interface sharpness
//!
//! ### Geometric Advection
//! **Theorem (Volume Conservation)**: Geometric advection preserves
//! the total volume of each phase to machine precision:
//!
//! ```math
//! ∫_{Ω} α dΩ = constant  (for each phase)
//! ```
//!
//! ## Algorithm Overview
//!
//! 1. **Interface Reconstruction**: PLIC method to reconstruct interface geometry
//! 2. **Geometric Advection**: Flux calculation through cell faces
//! 3. **Volume Update**: Conservative volume fraction update
//! 4. **Surface Tension**: Continuum Surface Force (CSF) model
//! 5. **Compression**: Interface compression to prevent smearing
//!
//! ## Interface Reconstruction
//!
//! ### PLIC Algorithm (Youngs, 1982)
//! 1. **Normal Calculation**: ∇α / |∇α| using mixed Youngs-centered scheme
//! 2. **Plane Equation**: Find linear interface that matches volume fraction
//! 3. **Root Finding**: Solve nonlinear equation for interface position
//!
//! ### Geometric Flux Calculation
//! For each cell face, compute the volume flux using polygon clipping:
//! - **Face Geometry**: Compute intersection of interface plane with cell face
//! - **Flux Volume**: Calculate volume of intersected polygon
//! - **Conservation**: Exact volume conservation per time step
//!
//! ## Surface Tension Implementation
//!
//! ### Continuum Surface Force (CSF) Model (Brackbill et al., 1992)
//! ```math
//! F_σ = σ κ ∇α
//! ```
//!
//! where κ is the interface curvature computed from the divergence of the normal:
//! ```math
//! κ = ∇·(∇α / |∇α|)
//! ```
//!
//! ## Accuracy and Conservation
//!
//! **Theorem (VOF Conservation)**: The geometric VOF method is exactly
//! conservative for incompressible flows and preserves phase volumes
//! to machine precision, independent of time step size.
//!
//! **Interface Sharpness**: PLIC reconstruction maintains sub-cell
//! interface resolution without artificial smearing.
//!
//! ## Implementation Features
//!
//! - **PLIC Reconstruction**: Exact volume matching with Newton iteration
//! - **Geometric Advection**: Polygon clipping for exact conservation
//! - **Cache Blocking**: Optimized memory access patterns for 3D grids
//! - **SIMD Operations**: Vectorized normal and curvature calculations
//!
//! ## References
//!
//! - Hirt, C.W. & Nichols, B.D. (1981). "Volume of fluid (VOF) method for the dynamics of free boundaries"
//! - Youngs, D.L. (1982). "Time-dependent multi-material flow with large fluid distortion"
//! - Brackbill, J.U. et al. (1992). "A continuum method for modeling surface tension"
//! - Scardovelli, R. & Zaleski, S. (1999). *Direct Numerical Simulation of Free-Surface and Interfacial Flow*

use cfd_core::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

use super::advection::AdvectionMethod;
use super::config::VofConfig;
use super::initialization::Initialization;
use super::reconstruction::InterfaceReconstruction;

/// VOF solver for multiphase flow
pub struct VofSolver<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> {
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

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> VofSolver<T> {
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
    #[inline]
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

    /// Set volume fraction field
    pub fn set_volume_fraction(&mut self, alpha: Vec<T>) -> Result<()> {
        if alpha.len() != self.alpha.len() {
            return Err(cfd_core::error::Error::DimensionMismatch {
                expected: self.alpha.len(),
                actual: alpha.len(),
            });
        }
        self.alpha = alpha;
        Ok(())
    }

    /// Get interface normals
    pub fn get_normals(&self) -> &[Vector3<T>] {
        &self.normals
    }

    /// Get interface curvature
    pub fn get_curvature(&self) -> &[T] {
        &self.curvature
    }

    /// Copy boundary values from alpha to alpha_previous
    /// This is used to initialize the boundary conditions in the temporary buffer
    /// without copying the entire field, which provides a performance optimization.
    pub fn copy_boundaries(&mut self) {
        // Copy Z boundaries (k=0 and k=nz-1)
        // These are contiguous blocks of memory (assuming Z is slowest index)
        let plane_size = self.nx * self.ny;
        let total_size = self.alpha.len();

        // k=0 plane
        if total_size >= plane_size {
            self.alpha_previous[0..plane_size].copy_from_slice(&self.alpha[0..plane_size]);
        }

        // k=nz-1 plane
        if total_size >= plane_size {
            let start = total_size - plane_size;
            self.alpha_previous[start..total_size].copy_from_slice(&self.alpha[start..total_size]);
        }

        // Copy Y boundaries (j=0 and j=ny-1) for interior k
        if self.nz > 2 {
            for k in 1..self.nz - 1 {
                // j=0 row
                let idx_j0 = self.index(0, 0, k);
                let range_j0 = idx_j0..idx_j0 + self.nx;
                self.alpha_previous[range_j0.clone()].copy_from_slice(&self.alpha[range_j0]);

                // j=ny-1 row
                let idx_jend = self.index(0, self.ny - 1, k);
                let range_jend = idx_jend..idx_jend + self.nx;
                self.alpha_previous[range_jend.clone()].copy_from_slice(&self.alpha[range_jend]);
            }
        }

        // Copy X boundaries (i=0 and i=nx-1) for interior k, interior j
        if self.nz > 2 && self.ny > 2 {
            for k in 1..self.nz - 1 {
                for j in 1..self.ny - 1 {
                    let idx_i0 = self.index(0, j, k);
                    self.alpha_previous[idx_i0] = self.alpha[idx_i0];

                    let idx_iend = self.index(self.nx - 1, j, k);
                    self.alpha_previous[idx_iend] = self.alpha[idx_iend];
                }
            }
        }
    }

    /// Main time step
    pub fn advance(&mut self, dt: T) -> Result<()> {
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
            let dx_min = <T as num_traits::Float>::min(<T as num_traits::Float>::min(self.dx, self.dy), self.dz);
            <T as FromPrimitive>::from_f64(self.config.cfl_number).unwrap_or(<T as FromPrimitive>::from_f64(0.3).unwrap_or(T::one()))
                * dx_min
                / max_velocity
        } else {
            <T as FromPrimitive>::from_f64(1e-3).unwrap_or(T::one())
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
