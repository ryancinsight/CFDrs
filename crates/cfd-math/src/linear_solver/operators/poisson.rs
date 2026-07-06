use crate::linear_solver::traits::LinearOperator;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::Array1;

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

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

    fn apply_boundary_conditions(&self, y: &mut Array1<T>) {
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

impl<T: RealField + Copy + FloatElement> LinearOperator<T> for LaplacianOperator2D<T> {
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> Result<()> {
        if x.shape()[0] != self.size() || y.shape()[0] != self.size() {
            return Err(Error::InvalidConfiguration(
                "Vector dimensions don't match operator size".to_string(),
            ));
        }

        let dx2_inv = <T as NumericElement>::ONE / (self.dx * self.dx);
        let dy2_inv = <T as NumericElement>::ONE / (self.dy * self.dy);
        let two: T = from_f64(2.0);

        for j in 1..(self.ny - 1) {
            for i in 1..(self.nx - 1) {
                let idx = j * self.nx + i;
                let center = x[idx];

                let d2x_dx2 = (x[idx - 1] + x[idx + 1] - two * center) * dx2_inv;
                let d2x_dy2 = (x[idx - self.nx] + x[idx + self.nx] - two * center) * dy2_inv;

                y[idx] = -(d2x_dx2 + d2x_dy2);
            }
        }

        self.apply_boundary_conditions(y);
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
        Self {
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
        }
    }
}

impl<T: RealField + Copy + FloatElement> LinearOperator<T> for PoissonOperator3D<T> {
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> Result<()> {
        let n = self.size();
        if x.shape()[0] != n || y.shape()[0] != n {
            return Err(Error::InvalidConfiguration("Dimensions mismatch".into()));
        }

        let dx2_inv = <T as NumericElement>::ONE / (self.dx * self.dx);
        let dy2_inv = <T as NumericElement>::ONE / (self.dy * self.dy);
        let dz2_inv = <T as NumericElement>::ONE / (self.dz * self.dz);
        let two: T = from_f64(2.0);

        for k in 1..(self.nz - 1) {
            for j in 1..(self.ny - 1) {
                for i in 1..(self.nx - 1) {
                    let idx = (k * self.ny + j) * self.nx + i;
                    let center = x[idx];

                    let d2x = (x[idx - 1] + x[idx + 1] - two * center) * dx2_inv;
                    let d2y = (x[idx - self.nx] + x[idx + self.nx] - two * center) * dy2_inv;
                    let d2z = (x[idx - self.nx * self.ny] + x[idx + self.nx * self.ny]
                        - two * center)
                        * dz2_inv;

                    y[idx] = -(d2x + d2y + d2z);
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

#[cfg(test)]
mod tests {
    use super::*;
    use leto::Array1;

    #[test]
    fn laplacian_center_impulse_matches_five_point_stencil() -> Result<()> {
        let operator = LaplacianOperator2D::new(3, 3, 1.0, 1.0);
        let x = Array1::from_shape_vec([9], vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
            .expect("valid shape");
        let mut y = Array1::zeros([9]);

        operator.apply(&x, &mut y)?;

        assert_eq!(y[4], 4.0);
        Ok(())
    }

    #[test]
    fn poisson_center_impulse_matches_seven_point_stencil() -> Result<()> {
        let operator = PoissonOperator3D::new(3, 3, 3, 1.0, 1.0, 1.0);
        let mut values = vec![0.0; 27];
        values[13] = 1.0;
        let x = Array1::from_shape_vec([27], values).expect("valid shape");
        let mut y = Array1::zeros([27]);

        operator.apply(&x, &mut y)?;

        assert_eq!(y[13], 6.0);
        Ok(())
    }
}
