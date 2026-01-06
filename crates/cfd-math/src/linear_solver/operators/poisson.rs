use crate::linear_solver::traits::LinearOperator;
use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};

/// 2D Laplacian operator using finite differences.
pub struct LaplacianOperator2D<T: RealField + Copy> {
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
}

impl<T: RealField + Copy> LaplacianOperator2D<T> {
    /// Create a new 2D Laplacian operator
    pub fn new(nx: usize, ny: usize, dx: T, dy: T) -> Self {
        Self { nx, ny, dx, dy }
    }

    fn apply_boundary_conditions(&self, y: &mut [T]) {
        for j in 0..self.ny {
            let left_idx = j * self.nx;
            let right_idx = j * self.nx + (self.nx - 1);
            if self.nx > 1 {
                y[left_idx] = y[left_idx + 1];
                y[right_idx] = y[right_idx - 1];
            }
        }
        for i in 0..self.nx {
            let bottom_idx = i;
            let top_idx = (self.ny - 1) * self.nx + i;
            if self.ny > 1 {
                y[bottom_idx] = y[bottom_idx + self.nx];
                y[top_idx] = y[top_idx - self.nx];
            }
        }
    }
}

impl<T: RealField + Copy + num_traits::FromPrimitive> LinearOperator<T> for LaplacianOperator2D<T> {
    fn apply(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<()> {
        if x.len() != self.size() || y.len() != self.size() {
            return Err(Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let dx2_inv = T::one() / (self.dx * self.dx);
        let dy2_inv = T::one() / (self.dy * self.dy);
        let two = T::from_f64(2.0).unwrap();

        let x_slice = x.as_slice();
        let y_slice = y.as_mut_slice();

        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x_slice[idx];

                let d2x_dx2 = (x_slice[idx - 1] + x_slice[idx + 1] - two * center) * dx2_inv;
                let d2x_dy2 = (x_slice[idx - self.nx] + x_slice[idx + self.nx] - two * center) * dy2_inv;

                y_slice[idx] = -(d2x_dx2 + d2x_dy2);
            }
        }

        self.apply_boundary_conditions(y_slice);
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny
    }

    fn is_symmetric(&self) -> bool {
        true
    }
}

/// 3D Poisson operator using finite differences.
pub struct PoissonOperator3D<T: RealField + Copy> {
    nx: usize,
    ny: usize,
    nz: usize,
    dx: T,
    dy: T,
    dz: T,
}

impl<T: RealField + Copy> PoissonOperator3D<T> {
    /// Create a new 3D Poisson operator
    pub fn new(nx: usize, ny: usize, nz: usize, dx: T, dy: T, dz: T) -> Self {
        Self { nx, ny, nz, dx, dy, dz }
    }
}

impl<T: RealField + Copy + num_traits::FromPrimitive> LinearOperator<T> for PoissonOperator3D<T> {
    fn apply(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<()> {
        let n = self.size();
        if x.len() != n || y.len() != n {
            return Err(Error::InvalidConfiguration("Dimensions mismatch".into()));
        }

        let dx2_inv = T::one() / (self.dx * self.dx);
        let dy2_inv = T::one() / (self.dy * self.dy);
        let dz2_inv = T::one() / (self.dz * self.dz);
        let two = T::from_f64(2.0).unwrap();

        let x_s = x.as_slice();
        let y_s = y.as_mut_slice();

        for k in 1..(self.nz - 1) {
            for j in 1..(self.ny - 1) {
                for i in 1..(self.nx - 1) {
                    let idx = (k * self.ny + j) * self.nx + i;
                    let center = x_s[idx];

                    let d2x = (x_s[idx - 1] + x_s[idx + 1] - two * center) * dx2_inv;
                    let d2y = (x_s[idx - self.nx] + x_s[idx + self.nx] - two * center) * dy2_inv;
                    let d2z = (x_s[idx - self.nx * self.ny] + x_s[idx + self.nx * self.ny] - two * center) * dz2_inv;

                    y_s[idx] = -(d2x + d2y + d2z);
                }
            }
        }
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny * self.nz
    }

    fn is_symmetric(&self) -> bool {
        true
    }
}
