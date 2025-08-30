//! Accelerated solvers using SIMD/GPU when available
//!
//! Provides unified interface that automatically uses best available acceleration

use crate::error::{Error, Result};
use crate::fields::Field2D;
use crate::grid::StructuredGrid2D;
use nalgebra::RealField;
use std::sync::Arc;

/// Acceleration backend selection
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum Backend {
    /// GPU acceleration via wgpu
    Gpu,
    /// CPU SIMD acceleration
    Simd,
    /// Standard CPU (fallback)
    Cpu,
}

/// Unified accelerated Poisson solver
pub struct AcceleratedPoissonSolver {
    backend: Backend,
    #[cfg(feature = "gpu")]
    gpu_solver: Option<Arc<cfd_core::compute::gpu::GpuPoissonSolver>>,
}

impl AcceleratedPoissonSolver {
    /// Create solver with automatic backend selection
    pub fn new(nx: usize, ny: usize, dx: f32, dy: f32) -> Result<Self> {
        // Try GPU first if feature enabled
        #[cfg(feature = "gpu")]
        {
            if let Ok(gpu_context) = cfd_core::compute::gpu::GpuContext::create() {
                if let Ok(gpu_solver) = cfd_core::compute::gpu::GpuPoissonSolver::new(
                    gpu_context.device.clone(),
                    gpu_context.queue.clone(),
                    nx,
                    ny,
                    dx,
                    dy,
                ) {
                    return Ok(Self {
                        backend: Backend::Gpu,
                        gpu_solver: Some(Arc::new(gpu_solver)),
                    });
                }
            }
        }

        // Check for SIMD capability
        let simd_cap = cfd_math::simd::SimdCapability::detect();
        let backend = match simd_cap {
            cfd_math::simd::SimdCapability::Avx2
            | cfd_math::simd::SimdCapability::Sse42
            | cfd_math::simd::SimdCapability::Neon => Backend::Simd,
            _ => Backend::Cpu,
        };

        Ok(Self {
            backend,
            #[cfg(feature = "gpu")]
            gpu_solver: None,
        })
    }

    /// Solve Poisson equation with automatic acceleration
    pub fn solve(
        &self,
        phi: &mut Field2D<f32>,
        source: &Field2D<f32>,
        iterations: usize,
        omega: f32,
    ) -> Result<f32> {
        match self.backend {
            #[cfg(feature = "gpu")]
            Backend::Gpu => {
                if let Some(ref gpu_solver) = self.gpu_solver {
                    // Use GPU solver
                    let residual = gpu_solver.solve_jacobi(
                        phi.as_mut_slice(),
                        source.as_slice(),
                        iterations,
                        omega,
                    )?;
                    Ok(residual)
                } else {
                    self.solve_simd(phi, source, iterations, omega)
                }
            }
            Backend::Simd => self.solve_simd(phi, source, iterations, omega),
            Backend::Cpu => self.solve_cpu(phi, source, iterations, omega),
            #[cfg(not(feature = "gpu"))]
            Backend::Gpu => {
                // GPU not available, fall back to CPU
                tracing::warn!("GPU backend requested but not available, using CPU");
                self.solve_cpu(phi, source, iterations, omega)
            }
        }
    }

    /// SIMD-accelerated solver
    fn solve_simd(
        &self,
        phi: &mut Field2D<f32>,
        source: &Field2D<f32>,
        iterations: usize,
        omega: f32,
    ) -> Result<f32> {
        let nx = phi.nx();
        let ny = phi.ny();
        let dx = 1.0 / (nx as f32 - 1.0);
        let dy = 1.0 / (ny as f32 - 1.0);

        for _ in 0..iterations {
            super::simd_kernels::gauss_seidel_simd(
                phi.as_mut_slice(),
                source.as_slice(),
                nx,
                ny,
                dx,
                dy,
                omega,
            )?;
        }

        // Calculate residual
        let mut residual_field = vec![0.0f32; nx * ny];
        super::simd_kernels::calculate_residual_simd(
            phi.as_slice(),
            source.as_slice(),
            &mut residual_field,
            nx,
            ny,
            dx,
            dy,
        )
        .map_err(|e| Error::from(format!("Failed to calculate residual: {:?}", e)))
    }

    /// Standard CPU solver (fallback)
    fn solve_cpu(
        &self,
        phi: &mut Field2D<f32>,
        source: &Field2D<f32>,
        iterations: usize,
        omega: f32,
    ) -> Result<f32> {
        let nx = phi.nx();
        let ny = phi.ny();
        let dx = 1.0 / (nx as f32 - 1.0);
        let dy = 1.0 / (ny as f32 - 1.0);

        let dx2 = dx * dx;
        let dy2 = dy * dy;
        let factor = omega / (2.0 * (1.0 / dx2 + 1.0 / dy2));

        let mut max_residual = 0.0f32;

        for _ in 0..iterations {
            // Red-Black Gauss-Seidel
            for color in 0..2 {
                for i in 1..nx - 1 {
                    for j in 1..ny - 1 {
                        if (i + j) % 2 == color {
                            let residual = (phi[(i - 1, j)] + phi[(i + 1, j)]) / dx2
                                + (phi[(i, j - 1)] + phi[(i, j + 1)]) / dy2
                                - source[(i, j)];

                            phi[(i, j)] = (1.0 - omega) * phi[(i, j)] + factor * residual;
                        }
                    }
                }
            }

            // Calculate max residual for convergence check
            max_residual = 0.0;
            for i in 1..nx - 1 {
                for j in 1..ny - 1 {
                    let laplacian = (phi[(i - 1, j)] - 2.0 * phi[(i, j)] + phi[(i + 1, j)]) / dx2
                        + (phi[(i, j - 1)] - 2.0 * phi[(i, j)] + phi[(i, j + 1)]) / dy2;

                    let residual = (laplacian - source[(i, j)]).abs();
                    max_residual = max_residual.max(residual);
                }
            }
        }

        Ok(max_residual)
    }

    /// Get active backend
    pub fn backend(&self) -> Backend {
        self.backend
    }
}

/// Accelerated velocity-pressure solver
pub struct AcceleratedNavierStokesSolver {
    poisson_solver: AcceleratedPoissonSolver,
    backend: Backend,
}

impl AcceleratedNavierStokesSolver {
    /// Create new accelerated Navier-Stokes solver
    pub fn new(grid: &StructuredGrid2D<f32>) -> Result<Self> {
        let nx = grid.nx;
        let ny = grid.ny;
        let (dx, dy) = grid.spacing();

        let poisson_solver = AcceleratedPoissonSolver::new(nx, ny, dx, dy)?;
        let backend = poisson_solver.backend();

        Ok(Self {
            poisson_solver,
            backend,
        })
    }

    /// Solve pressure Poisson equation
    pub fn solve_pressure(
        &self,
        pressure: &mut Field2D<f32>,
        divergence: &Field2D<f32>,
        iterations: usize,
    ) -> Result<f32> {
        self.poisson_solver.solve(
            pressure,
            divergence,
            iterations,
            cfd_core::constants::numerical::relaxation::SOR_OMEGA_DEFAULT as f32,
        )
    }

    /// Calculate velocity divergence with acceleration
    pub fn calculate_divergence(
        &self,
        u: &Field2D<f32>,
        v: &Field2D<f32>,
        divergence: &mut Field2D<f32>,
    ) -> Result<()> {
        let nx = u.nx();
        let ny = u.ny();
        let dx = 1.0 / (nx as f32 - 1.0);
        let dy = 1.0 / (ny as f32 - 1.0);

        match self.backend {
            Backend::Simd => super::simd_kernels::calculate_divergence_simd(
                u.as_slice(),
                v.as_slice(),
                divergence.as_mut_slice(),
                nx,
                ny,
                dx,
                dy,
            )
            .map_err(|e| Error::from(format!("Failed to calculate divergence: {:?}", e))),
            _ => {
                // CPU fallback
                for i in 1..nx - 1 {
                    for j in 1..ny - 1 {
                        let dudx = (u[(i + 1, j)] - u[(i - 1, j)]) / (2.0 * dx);
                        let dvdy = (v[(i, j + 1)] - v[(i, j - 1)]) / (2.0 * dy);
                        divergence[(i, j)] = dudx + dvdy;
                    }
                }
                Ok(())
            }
        }
    }

    /// Get active backend
    pub fn backend(&self) -> Backend {
        self.backend
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_backend_selection() {
        let solver = AcceleratedPoissonSolver::new(100, 100, 0.01, 0.01).unwrap();

        // Should select SIMD or GPU if available
        match solver.backend() {
            Backend::Gpu => println!("Using GPU acceleration"),
            Backend::Simd => println!("Using SIMD acceleration"),
            Backend::Cpu => println!("Using CPU fallback"),
        }
    }

    #[test]
    fn test_poisson_solve() {
        let nx = 50;
        let ny = 50;
        let mut phi = Field2D::zeros(nx, ny);
        let mut source = Field2D::zeros(nx, ny);

        // Set source term
        source[(nx / 2, ny / 2)] = 1.0;

        let solver = AcceleratedPoissonSolver::new(nx, ny, 0.02, 0.02).unwrap();
        let residual = solver.solve(&mut phi, &source, 100, 1.0).unwrap();

        assert!(residual < 1.0); // Should converge
    }
}
