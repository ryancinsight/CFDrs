//! Constants for GPU compute operations

/// Workgroup dimensions
pub mod workgroup {
    /// Default 2D workgroup size (8x8 threads)
    pub const SIZE_2D: u32 = 8;

    /// Default 1D workgroup size for reductions
    pub const SIZE_1D: u32 = 256;

    /// Maximum workgroup size (hardware dependent, conservative default)
    pub const MAX_SIZE: u32 = 1024;
}

/// Computational complexity estimates (FLOPs per grid point)
pub mod complexity {
    /// Advection kernel FLOPs per point (upwind scheme)
    pub const ADVECTION_UPWIND: usize = 10;

    /// Diffusion kernel FLOPs per point (central difference)
    pub const DIFFUSION_CENTRAL: usize = 15;

    /// Pressure solver FLOPs per iteration per point (Jacobi)
    pub const PRESSURE_JACOBI: usize = 20;

    /// Velocity update FLOPs per point
    pub const VELOCITY_UPDATE: usize = 8;
}

/// Numerical parameters
pub mod numerical {
    /// Default relaxation factor for SOR (Successive Over-Relaxation)
    pub const SOR_OMEGA_DEFAULT: f32 = 1.2;

    /// Optimal SOR relaxation factor range
    pub const SOR_OMEGA_MIN: f32 = 1.0;
    /// Maximum omega value for SOR (Successive Over-Relaxation) iteration
    pub const SOR_OMEGA_MAX: f32 = 1.95;

    /// Default convergence tolerance
    pub const CONVERGENCE_TOLERANCE: f32 = 1e-6;

    /// Maximum iterations for iterative solvers
    pub const MAX_ITERATIONS: u32 = 1000;
}

/// Memory alignment requirements
pub mod alignment {
    /// Alignment for buffer data (16 bytes for SIMD)
    pub const BUFFER_ALIGNMENT: usize = 16;

    /// Uniform buffer alignment requirement (256 bytes typical)
    pub const UNIFORM_ALIGNMENT: usize = 256;
}
