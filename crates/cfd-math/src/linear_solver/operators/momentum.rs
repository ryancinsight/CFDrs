use crate::linear_solver::traits::LinearOperator;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

/// 1D Momentum operator (advection-diffusion).
pub struct MomentumOperator1D<T: RealField + Copy> {
    n: usize,
    dx: T,
    u: T,  // Velocity
    nu: T, // Viscosity
}

impl<T: RealField + Copy> MomentumOperator1D<T> {
    /// Create a new 1D momentum operator
    pub fn new(n: usize, dx: T, u: T, nu: T) -> Self {
        Self { n, dx, u, nu }
    }
}

impl<T: RealField + Copy + FloatElement> LinearOperator<T> for MomentumOperator1D<T> {
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> Result<()> {
        if x.shape()[0] != self.n || y.shape()[0] != self.n {
            return Err(Error::InvalidConfiguration("Dimensions mismatch".into()));
        }

        let dx_inv = <T as NumericElement>::ONE / self.dx;
        let dx2_inv = <T as NumericElement>::ONE / (self.dx * self.dx);
        let two: T = from_f64(2.0);
        let half: T = from_f64(0.5);

        for i in 1..(self.n - 1) {
            // Central difference for advection and diffusion
            let advection = self.u * (x[i + 1] - x[i - 1]) * half * dx_inv;
            let diffusion = self.nu * (x[i - 1] + x[i + 1] - two * x[i]) * dx2_inv;
            y[i] = advection - diffusion;
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
        Self {
            nx,
            ny,
            dx,
            dy,
            u,
            v,
            nu,
        }
    }
}

impl<T: RealField + Copy + FloatElement> LinearOperator<T> for MomentumOperator2D<T> {
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> Result<()> {
        let n = self.size();
        if x.shape()[0] != n || y.shape()[0] != n {
            return Err(Error::InvalidConfiguration("Dimensions mismatch".into()));
        }

        let dx_inv = <T as NumericElement>::ONE / self.dx;
        let dy_inv = <T as NumericElement>::ONE / self.dy;
        let dx2_inv = <T as NumericElement>::ONE / (self.dx * self.dx);
        let dy2_inv = <T as NumericElement>::ONE / (self.dy * self.dy);
        let two: T = from_f64(2.0);
        let half: T = from_f64(0.5);

        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x[idx];

                let adv_x = self.u * (x[idx + 1] - x[idx - 1]) * half * dx_inv;
                let adv_y = self.v * (x[idx + self.nx] - x[idx - self.nx]) * half * dy_inv;

                let diff_x = self.nu * (x[idx - 1] + x[idx + 1] - two * center) * dx2_inv;
                let diff_y =
                    self.nu * (x[idx - self.nx] + x[idx + self.nx] - two * center) * dy2_inv;

                y[idx] = adv_x + adv_y - (diff_x + diff_y);
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
        Self {
            nx,
            ny,
            dx,
            dy,
            u,
            v,
            alpha,
        }
    }
}

impl<T: RealField + Copy + FloatElement> LinearOperator<T> for EnergyOperator2D<T> {
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> Result<()> {
        let n = self.size();
        if x.shape()[0] != n || y.shape()[0] != n {
            return Err(Error::InvalidConfiguration("Dimensions mismatch".into()));
        }

        let dx_inv = <T as NumericElement>::ONE / self.dx;
        let dy_inv = <T as NumericElement>::ONE / self.dy;
        let dx2_inv = <T as NumericElement>::ONE / (self.dx * self.dx);
        let dy2_inv = <T as NumericElement>::ONE / (self.dy * self.dy);
        let two: T = from_f64(2.0);
        let half: T = from_f64(0.5);

        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x[idx];

                let adv_x = self.u * (x[idx + 1] - x[idx - 1]) * half * dx_inv;
                let adv_y = self.v * (x[idx + self.nx] - x[idx - self.nx]) * half * dy_inv;

                let diff_x = self.alpha * (x[idx - 1] + x[idx + 1] - two * center) * dx2_inv;
                let diff_y =
                    self.alpha * (x[idx - self.nx] + x[idx + self.nx] - two * center) * dy2_inv;

                y[idx] = adv_x + adv_y - (diff_x + diff_y);
            }
        }
        Ok(())
    }

    fn size(&self) -> usize {
        self.nx * self.ny
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use leto::Array1;

    #[test]
    fn momentum_1d_applies_central_advection_diffusion() -> Result<()> {
        let operator = MomentumOperator1D::new(3, 1.0, 2.0, 0.5);
        let x = Array1::from_shape_vec([3], vec![1.0, 2.0, 4.0]).expect("valid shape");
        let mut y = Array1::zeros([3]);

        operator.apply(&x, &mut y)?;

        assert_eq!(y[1], 2.5);
        Ok(())
    }

    #[test]
    fn momentum_2d_applies_center_advection_diffusion() -> Result<()> {
        let operator = MomentumOperator2D::new(3, 3, 1.0, 1.0, 1.0, 2.0, 0.5);
        let x = Array1::from_shape_vec([9], (0..9).map(f64::from).collect()).expect("valid shape");
        let mut y = Array1::zeros([9]);

        operator.apply(&x, &mut y)?;

        assert_eq!(y[4], 7.0);
        Ok(())
    }

    #[test]
    fn energy_2d_applies_center_advection_diffusion() -> Result<()> {
        let operator = EnergyOperator2D::new(3, 3, 1.0, 1.0, 1.0, 2.0, 0.5);
        let x = Array1::from_shape_vec([9], (0..9).map(f64::from).collect()).expect("valid shape");
        let mut y = Array1::zeros([9]);

        operator.apply(&x, &mut y)?;

        assert_eq!(y[4], 7.0);
        Ok(())
    }
}
