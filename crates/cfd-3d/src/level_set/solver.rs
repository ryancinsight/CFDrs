//! Level set solver implementation

use super::config::LevelSetConfig;
use cfd_core::error::{Error, Result, NumericalErrorKind};
use nalgebra::{DMatrix, RealField};
use num_traits::FromPrimitive;
/// Level set solver
pub struct LevelSetSolver<T: RealField> {
    /// Configuration
    pub config: LevelSetConfig<T>,
    /// Level set function values
    pub phi: DMatrix<T>,
    /// Previous phi for time stepping
    pub phi_previous: DMatrix<T>,
    /// Grid dimensions
    pub nx: usize,
    pub ny: usize,
    pub nz: usize,
}
impl<T: RealField + FromPrimitive> LevelSetSolver<T> {
    /// Create new solver
    pub fn new(config: LevelSetConfig<T>, nx: usize, ny: usize, nz: usize) -> Self {
        let size = nx * ny * nz;
        Self {
            config,
            phi: DMatrix::zeros(size, 1),
            phi_previous: DMatrix::zeros(size, 1),
            nx,
            ny,
            nz,
        }
    }
    /// Initialize level set function
    pub fn initialize<F>(&mut self, init_fn: F) -> Result<()>
    where
        F: Fn(T, T, T) -> T,
    {
        for k in 0..self.nz {
            for j in 0..self.ny {
                for i in 0..self.nx {
                    let x = T::from_usize(i).ok_or_else(|| {
                        Error::Numerical(NumericalErrorKind::ConversionFailed {
                            from_type: "usize",
                            to_type: std::any::type_name::<T>(),
                            value: i.to_string(),
                        })
                    })? * self.config.dx;
                    let y = T::from_usize(j).ok_or_else(|| {
                            value: j.to_string(),
                    })? * self.config.dy;
                    let z = T::from_usize(k).ok_or_else(|| {
                            value: k.to_string(),
                    })? * self.config.dz;
                    
                    let idx = self.linear_index(i, j, k);
                    self.phi[(idx, 0)] = init_fn(x, y, z);
                }
            }
        self.phi_previous = self.phi.clone();
        Ok(())
    /// Get linear index from 3D indices
    #[inline]
    pub fn linear_index(&self, i: usize, j: usize, k: usize) -> usize {
        i + j * self.nx + k * self.nx * self.ny
    /// Time step the level set equation
    }

    pub fn step(&mut self, dt: T) -> Result<()> {
        // Actual stepping would be implemented here


    }

}
}
}
}
