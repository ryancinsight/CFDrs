use crate::linear_solver::traits::LinearOperator;
use cfd_core::error::{Error, Result};
use nalgebra::{DVector, RealField};

/// 1D Momentum operator (advection-diffusion).
pub struct MomentumOperator1D<T: RealField + Copy> {
    n: usize,
    dx: T,
    u: T, // Velocity
    nu: T, // Viscosity
}

impl<T: RealField + Copy> MomentumOperator1D<T> {
    /// Create a new 1D momentum operator
    pub fn new(n: usize, dx: T, u: T, nu: T) -> Self {
        Self { n, dx, u, nu }
    }
}

impl<T: RealField + Copy + num_traits::FromPrimitive> LinearOperator<T> for MomentumOperator1D<T> {
    fn apply(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<()> {
        if x.len() != self.n || y.len() != self.n {
            return Err(Error::InvalidConfiguration("Dimensions mismatch".into()));
        }

        let dx_inv = T::one() / self.dx;
        let dx2_inv = T::one() / (self.dx * self.dx);
        let two = T::from_f64(2.0).unwrap();
        let half = T::from_f64(0.5).unwrap();

        let x_s = x.as_slice();
        let y_s = y.as_mut_slice();

        for i in 1..(self.n - 1) {
            // Central difference for advection and diffusion
            let advection = self.u * (x_s[i + 1] - x_s[i - 1]) * half * dx_inv;
            let diffusion = self.nu * (x_s[i - 1] + x_s[i + 1] - two * x_s[i]) * dx2_inv;
            y_s[i] = advection - diffusion;
        }

        Ok(())
    }

    fn size(&self) -> usize {
        self.n
    }
}

/// 2D Momentum operator.
pub struct MomentumOperator2D<T: RealField + Copy> {
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    u: T,
    v: T,
    nu: T,
}

impl<T: RealField + Copy> MomentumOperator2D<T> {
    /// Create a new 2D momentum operator
    pub fn new(nx: usize, ny: usize, dx: T, dy: T, u: T, v: T, nu: T) -> Self {
        Self { nx, ny, dx, dy, u, v, nu }
    }
}

impl<T: RealField + Copy + num_traits::FromPrimitive> LinearOperator<T> for MomentumOperator2D<T> {
    fn apply(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<()> {
        let n = self.size();
        if x.len() != n || y.len() != n {
            return Err(Error::InvalidConfiguration("Dimensions mismatch".into()));
        }

        let dx_inv = T::one() / self.dx;
        let dy_inv = T::one() / self.dy;
        let dx2_inv = T::one() / (self.dx * self.dx);
        let dy2_inv = T::one() / (self.dy * self.dy);
        let two = T::from_f64(2.0).unwrap();
        let half = T::from_f64(0.5).unwrap();

        let x_s = x.as_slice();
        let y_s = y.as_mut_slice();

        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x_s[idx];

                let adv_x = self.u * (x_s[idx + 1] - x_s[idx - 1]) * half * dx_inv;
                let adv_y = self.v * (x_s[idx + self.nx] - x_s[idx - self.nx]) * half * dy_inv;

                let diff_x = self.nu * (x_s[idx - 1] + x_s[idx + 1] - two * center) * dx2_inv;
                let diff_y = self.nu * (x_s[idx - self.nx] + x_s[idx + self.nx] - two * center) * dy2_inv;

                y_s[idx] = adv_x + adv_y - (diff_x + diff_y);
            }
        }
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny
    }
}

/// 2D Energy operator.
pub struct EnergyOperator2D<T: RealField + Copy> {
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    u: T,
    v: T,
    alpha: T, // Thermal diffusivity
}

impl<T: RealField + Copy> EnergyOperator2D<T> {
    /// Create a new 2D energy operator
    pub fn new(nx: usize, ny: usize, dx: T, dy: T, u: T, v: T, alpha: T) -> Self {
        Self { nx, ny, dx, dy, u, v, alpha }
    }
}

impl<T: RealField + Copy + num_traits::FromPrimitive> LinearOperator<T> for EnergyOperator2D<T> {
    fn apply(&self, x: &DVector<T>, y: &mut DVector<T>) -> Result<()> {
        let n = self.size();
        if x.len() != n || y.len() != n {
            return Err(Error::InvalidConfiguration("Dimensions mismatch".into()));
        }

        let dx_inv = T::one() / self.dx;
        let dy_inv = T::one() / self.dy;
        let dx2_inv = T::one() / (self.dx * self.dx);
        let dy2_inv = T::one() / (self.dy * self.dy);
        let two = T::from_f64(2.0).unwrap();
        let half = T::from_f64(0.5).unwrap();

        let x_s = x.as_slice();
        let y_s = y.as_mut_slice();

        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x_s[idx];

                let adv_x = self.u * (x_s[idx + 1] - x_s[idx - 1]) * half * dx_inv;
                let adv_y = self.v * (x_s[idx + self.nx] - x_s[idx - self.nx]) * half * dy_inv;

                let diff_x = self.alpha * (x_s[idx - 1] + x_s[idx + 1] - two * center) * dx2_inv;
                let diff_y = self.alpha * (x_s[idx - self.nx] + x_s[idx + self.nx] - two * center) * dy2_inv;

                y_s[idx] = adv_x + adv_y - (diff_x + diff_y);
            }
        }
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny
    }
}
