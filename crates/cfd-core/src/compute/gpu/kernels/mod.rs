//! GPU compute kernels for CFD operations

use hephaestus_wgpu::DispatchGrid;
pub mod advection;
pub mod diffusion;
pub mod laplacian;
pub mod pressure;
pub mod turbulence;
pub mod velocity;

pub use advection::{AdvectionConfig, GpuAdvectionKernel};
pub use diffusion::{DiffusionConfig, GpuDiffusionKernel};
pub use laplacian::Laplacian2DKernel;
pub use pressure::{GpuPressureKernel, PressureConfig};
pub use velocity::{GpuVelocityKernel, VelocityConfig};

use crate::error::{Error, Result};

fn validate_field_len(expected: usize, actual: usize) -> Result<()> {
    if actual == expected {
        Ok(())
    } else {
        Err(Error::DimensionMismatch { expected, actual })
    }
}

fn validate_finite_field(name: &str, values: &[f32]) -> Result<()> {
    if let Some(index) = values.iter().position(|value| !value.is_finite()) {
        Err(Error::PhysicsViolation(format!(
            "{name} field is non-finite at index {index}"
        )))
    } else {
        Ok(())
    }
}

fn dispatch_grid_3d(dimensions: [u32; 4], workgroup: [usize; 3]) -> Result<DispatchGrid> {
    let [nx, ny, nz, _] = dimensions;
    Ok(DispatchGrid::covering_domain(
        [
            usize::try_from(nx).expect("invariant: u32 fits usize"),
            usize::try_from(ny).expect("invariant: u32 fits usize"),
            usize::try_from(nz).expect("invariant: u32 fits usize"),
        ],
        workgroup,
    )?)
}
