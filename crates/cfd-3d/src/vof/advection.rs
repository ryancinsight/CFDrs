//! VOF advection methods

use cfd_core::error::Result;
use nalgebra::{Vector3, RealField};
use num_traits::FromPrimitive;
use super::solver::VofSolver;
use super::config::{VofConfig, VOF_INTERFACE_LOWER, VOF_INTERFACE_UPPER};

/// Advection method for VOF
pub struct AdvectionMethod {
    use_geometric: bool,
}

impl AdvectionMethod {
    /// Create advection method based on configuration
    pub fn create(config: &VofConfig) -> Self {
        Self {
            use_geometric: config.use_geometric_advection,
        }
    }

    /// Advect volume fraction field
    pub fn advect<T: RealField + FromPrimitive + Copy>(
        &self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
        if self.use_geometric {
            self.geometric_advection(solver, dt)
        } else {
            self.algebraic_advection(solver, dt)
        }
    }

    /// Geometric advection (more accurate, conserves volume)
    fn geometric_advection<T: RealField + FromPrimitive + Copy>(
        &self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
        let mut alpha_temp = solver.alpha.clone();

        for k in 1..solver.nz - 1 {
            for j in 1..solver.ny - 1 {
                for i in 1..solver.nx - 1 {
                    let idx = solver.index(i, j, k);
                    
                    // Calculate fluxes through all faces
                    let flux_x_minus = self.calculate_flux(solver, i - 1, j, k, i, j, k, 0, dt)?;
                    let flux_x_plus = self.calculate_flux(solver, i, j, k, i + 1, j, k, 0, dt)?;
                    let flux_y_minus = self.calculate_flux(solver, i, j - 1, k, i, j, k, 1, dt)?;
                    let flux_y_plus = self.calculate_flux(solver, i, j, k, i, j + 1, k, 1, dt)?;
                    let flux_z_minus = self.calculate_flux(solver, i, j, k - 1, i, j, k, 2, dt)?;
                    let flux_z_plus = self.calculate_flux(solver, i, j, k, i, j, k + 1, 2, dt)?;

                    // Update volume fraction
                    let cell_volume = solver.dx * solver.dy * solver.dz;
                    let net_flux = flux_x_plus - flux_x_minus
                        + flux_y_plus - flux_y_minus
                        + flux_z_plus - flux_z_minus;
                    
                    alpha_temp[idx] = solver.alpha[idx] - net_flux / cell_volume;
                    
                    // Bound volume fraction
                    alpha_temp[idx] = alpha_temp[idx].max(T::zero()).min(T::one());
                }
            }
        }

        solver.alpha = alpha_temp;
        Ok(())
    }

    /// Calculate flux through a face
    fn calculate_flux<T: RealField + FromPrimitive + Copy>(
        &self,
        solver: &VofSolver<T>,
        i1: usize,
        j1: usize,
        k1: usize,
        i2: usize,
        j2: usize,
        k2: usize,
        direction: usize, // 0=x, 1=y, 2=z
        dt: T,
    ) -> Result<T> {
        let idx1 = solver.index(i1, j1, k1);
        let idx2 = solver.index(i2, j2, k2);
        
        // Get velocity at face (simple average)
        let velocity = (solver.velocity[idx1] + solver.velocity[idx2]) * T::from_f64(0.5).unwrap_or(T::zero());
        
        // Get velocity component normal to face
        let u_normal = match direction {
            0 => velocity.x,
            1 => velocity.y,
            2 => velocity.z,
            _ => return Err(cfd_core::Error::InvalidInput("Invalid direction".to_string())),
        };
        
        // Upwind scheme for flux calculation
        let alpha_upwind = if u_normal > T::zero() {
            solver.alpha[idx1]
        } else {
            solver.alpha[idx2]
        };
        
        // Face area
        let face_area = match direction {
            0 => solver.dy * solver.dz,
            1 => solver.dx * solver.dz,
            2 => solver.dx * solver.dy,
            _ => return Err(cfd_core::Error::InvalidInput("Invalid direction".to_string())),
        };
        
        Ok(alpha_upwind * u_normal * face_area * dt)
    }

    /// Algebraic advection (simpler but less accurate)
    fn algebraic_advection<T: RealField + FromPrimitive + Copy>(
        &self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
        let mut alpha_temp = solver.alpha.clone();
        let two = T::from_f64(2.0).unwrap_or(T::one() + T::one());

        for k in 1..solver.nz - 1 {
            for j in 1..solver.ny - 1 {
                for i in 1..solver.nx - 1 {
                    let idx = solver.index(i, j, k);
                    let vel = solver.velocity[idx];

                    // Central differences for spatial derivatives
                    let dalpha_dx = (solver.alpha[solver.index(i + 1, j, k)]
                        - solver.alpha[solver.index(i - 1, j, k)])
                        / (two * solver.dx);
                    let dalpha_dy = (solver.alpha[solver.index(i, j + 1, k)]
                        - solver.alpha[solver.index(i, j - 1, k)])
                        / (two * solver.dy);
                    let dalpha_dz = (solver.alpha[solver.index(i, j, k + 1)]
                        - solver.alpha[solver.index(i, j, k - 1)])
                        / (two * solver.dz);

                    // Advection equation: ∂α/∂t + u·∇α = 0
                    alpha_temp[idx] = solver.alpha[idx]
                        - dt * (vel.x * dalpha_dx + vel.y * dalpha_dy + vel.z * dalpha_dz);

                    // Bound volume fraction
                    alpha_temp[idx] = alpha_temp[idx].max(T::zero()).min(T::one());
                }
            }
        }

        solver.alpha = alpha_temp;
        Ok(())
    }

    /// Apply artificial compression to sharpen interface
    pub fn apply_compression<T: RealField + FromPrimitive + Copy>(
        &self,
        solver: &mut VofSolver<T>,
        dt: T,
    ) -> Result<()> {
        let mut alpha_temp = solver.alpha.clone();

        for k in 1..solver.nz - 1 {
            for j in 1..solver.ny - 1 {
                for i in 1..solver.nx - 1 {
                    let idx = solver.index(i, j, k);
                    
                    // Only apply compression near interface
                    let alpha = solver.alpha[idx];
                    if alpha > T::from_f64(VOF_INTERFACE_LOWER).unwrap_or(T::zero())
                        && alpha < T::from_f64(VOF_INTERFACE_UPPER).unwrap_or(T::one())
                    {
                        // Get interface normal
                        let normal = solver.normals[idx];
                        
                        if normal.norm() > T::zero() {
                            // Compression velocity proportional to normal
                            let compression_factor = T::from_f64(0.5).unwrap_or(T::zero());
                            let u_compression = normal * compression_factor;
                            
                            // Apply compression flux
                            let two = T::from_f64(2.0).unwrap_or(T::one() + T::one());
                            let dalpha_dx = (solver.alpha[solver.index(i + 1, j, k)]
                                - solver.alpha[solver.index(i - 1, j, k)])
                                / (two * solver.dx);
                            let dalpha_dy = (solver.alpha[solver.index(i, j + 1, k)]
                                - solver.alpha[solver.index(i, j - 1, k)])
                                / (two * solver.dy);
                            let dalpha_dz = (solver.alpha[solver.index(i, j, k + 1)]
                                - solver.alpha[solver.index(i, j, k - 1)])
                                / (two * solver.dz);
                            
                            let compression_term = u_compression.x * dalpha_dx
                                + u_compression.y * dalpha_dy
                                + u_compression.z * dalpha_dz;
                            
                            alpha_temp[idx] = alpha - dt * compression_term * alpha * (T::one() - alpha);
                            
                            // Bound volume fraction
                            alpha_temp[idx] = alpha_temp[idx].max(T::zero()).min(T::one());
                        }
                    }
                }
            }
        }

        solver.alpha = alpha_temp;
        Ok(())
    }
}