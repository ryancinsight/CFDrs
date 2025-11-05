//! GPU acceleration for Smagorinsky LES model
//!
//! Provides GPU-accelerated computation of SGS viscosity using
//! WebGPU compute shaders.

#[cfg(feature = "gpu")]
use cfd_core::compute::gpu::turbulence_compute::GpuTurbulenceCompute;
use nalgebra::DMatrix;

/// GPU-accelerated SGS viscosity computation
#[cfg(feature = "gpu")]
pub fn compute_sgs_viscosity_gpu(
    gpu_compute: &mut GpuTurbulenceCompute,
    velocity_u: &DMatrix<f64>,
    velocity_v: &DMatrix<f64>,
    nx: usize,
    ny: usize,
    dx: f64,
    dy: f64,
    smagorinsky_constant: f64,
) -> cfd_core::error::Result<DMatrix<f64>> {
    // Convert f64 matrices to f32 for GPU computation
    let velocity_u_f32: Vec<f32> = velocity_u.iter().map(|&x| x as f32).collect();
    let velocity_v_f32: Vec<f32> = velocity_v.iter().map(|&x| x as f32).collect();

    // Compute SGS viscosity on GPU
    let sgs_buffer = gpu_compute.compute_smagorinsky_sgs(
        &velocity_u_f32,
        &velocity_v_f32,
        nx,
        ny,
        dx as f32,
        dy as f32,
        smagorinsky_constant as f32,
    )?;

    // Read back results and convert to f64
    let sgs_f32 = gpu_compute.read_buffer(&sgs_buffer)?;
    let mut sgs_viscosity = DMatrix::zeros(nx, ny);

    // Convert back to f64 matrix
    for (i, &val) in sgs_f32.iter().enumerate() {
        let j = i / nx;
        let ii = i % nx;
        sgs_viscosity[(ii, j)] = val as f64;
    }

    Ok(sgs_viscosity)
}

/// Check if GPU acceleration is available
#[cfg(feature = "gpu")]
pub fn gpu_acceleration_available() -> bool {
    GpuTurbulenceCompute::new().is_ok()
}

/// Fallback CPU implementation when GPU is not available
#[cfg(not(feature = "gpu"))]
pub fn compute_sgs_viscosity_gpu(
    _gpu_compute: &mut (),
    _velocity_u: &DMatrix<f64>,
    _velocity_v: &DMatrix<f64>,
    _nx: usize,
    _ny: usize,
    _dx: f64,
    _dy: f64,
    _smagorinsky_constant: f64,
) -> cfd_core::error::Result<DMatrix<f64>> {
    // This should never be called when GPU feature is not enabled
    Err(cfd_core::error::Error::UnsupportedOperation("GPU acceleration not available".to_string()))
}

/// Check if GPU acceleration is available (always false without GPU feature)
#[cfg(not(feature = "gpu"))]
pub fn gpu_acceleration_available() -> bool {
    false
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_gpu_availability() {
        let available = gpu_acceleration_available();

        #[cfg(feature = "gpu")]
        {
            // With GPU feature, availability depends on hardware
            // We can't assert much here without specific hardware
            let _ = available; // Just acknowledge the result
        }

        #[cfg(not(feature = "gpu"))]
        assert!(!available);
    }

    #[cfg(feature = "gpu")]
    #[test]
    fn test_gpu_sgs_computation() {
        if !gpu_acceleration_available() {
            return; // Skip test if GPU not available
        }

        let mut gpu_compute = GpuTurbulenceCompute::new().unwrap();
        let nx = 16;
        let ny = 16;

        // Create test velocity field
        let mut velocity_u = DMatrix::zeros(nx, ny);
        let mut velocity_v = DMatrix::zeros(nx, ny);

        // Simple shear flow
        for i in 0..nx {
            for j in 0..ny {
                velocity_u[(i, j)] = (j as f64) * 0.1;
                velocity_v[(i, j)] = 0.0;
            }
        }

        // Compute SGS viscosity on GPU
        let result = compute_sgs_viscosity_gpu(
            &mut gpu_compute,
            &velocity_u,
            &velocity_v,
            nx,
            ny,
            0.1,
            0.1,
            0.1, // Smagorinsky constant
        );

        assert!(result.is_ok());
        let sgs_viscosity = result.unwrap();

        // Check dimensions
        assert_eq!(sgs_viscosity.nrows(), nx);
        assert_eq!(sgs_viscosity.ncols(), ny);

        // Check that SGS viscosity is non-negative
        assert!(sgs_viscosity.iter().all(|&v| v >= 0.0));
    }
}
