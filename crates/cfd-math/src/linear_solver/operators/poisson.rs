use crate::linear_solver::traits::LinearOperator;
use aequitas::systems::si::quantities::Length;
use cfd_core::error::{Error, Result};
use eunomia::{FloatElement, NumericElement, RealField};
use leto::{Array1, BoundaryCondition, Laplacian2D, LaplacianPolarity};
use leto_ops::{laplacian_2d_into, RealScalar};

#[inline]
fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

/// Provider-backed negative two-dimensional Laplacian operator.
pub struct LaplacianOperator2D<T> {
    stencil: Laplacian2D<T>,
}

impl<T: RealScalar> LaplacianOperator2D<T> {
    /// Create a typed negative-Laplacian operator over a uniform Cartesian grid.
    ///
    /// # Errors
    ///
    /// Returns an invalid-configuration error when the Leto provider rejects
    /// the grid dimensions or physical spacing.
    pub fn new(
        nx: usize,
        ny: usize,
        dx: Length<T>,
        dy: Length<T>,
        boundary: BoundaryCondition,
    ) -> Result<Self> {
        let stencil = Laplacian2D::new(nx, ny, dx, dy, boundary)
            .map_err(|error| Error::InvalidConfiguration(error.to_string()))?
            .with_polarity(LaplacianPolarity::NegativeLaplacian);
        Ok(Self { stencil })
    }
}

impl<T: RealField + RealScalar> LinearOperator<T> for LaplacianOperator2D<T> {
    fn apply(&self, x: &Array1<T>, y: &mut Array1<T>) -> Result<()> {
        laplacian_2d_into(&self.stencil, &x.view(), &mut y.view_mut())
            .map_err(|error| Error::InvalidConfiguration(error.to_string()))
    }

    fn size(&self) -> usize {
        self.stencil.len()
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
    use aequitas::systems::si::units::Meter;
    use leto::Array1;

    #[test]
    fn laplacian_center_impulse_matches_five_point_stencil() -> Result<()> {
        let operator = LaplacianOperator2D::new(
            3,
            3,
            Length::from_unit::<Meter>(1.0),
            Length::from_unit::<Meter>(1.0),
            BoundaryCondition::Dirichlet,
        )?;
        let x = Array1::from_shape_vec([9], vec![0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0])
            .expect("valid shape");
        let mut y = Array1::zeros([9]);

        operator.apply(&x, &mut y)?;

        assert_eq!(y[4], 4.0);
        Ok(())
    }

    #[test]
    fn laplacian_neumann_quadratic_matches_closed_form_at_every_point() -> Result<()> {
        let nx = 4;
        let ny = 4;
        let dx = 0.5f32;
        let dy = 0.25f32;
        let values = (0..ny)
            .flat_map(|j| {
                (0..nx).map(move |i| {
                    let x = f32::from(u16::try_from(i).expect("small test index")) * dx;
                    let y = f32::from(u16::try_from(j).expect("small test index")) * dy;
                    x * x + 3.0 * y * y
                })
            })
            .collect();
        let input = Array1::from_shape_vec([nx * ny], values).expect("valid shape");
        let mut output = Array1::zeros([nx * ny]);
        let operator = LaplacianOperator2D::new(
            nx,
            ny,
            Length::from_unit::<Meter>(dx),
            Length::from_unit::<Meter>(dy),
            BoundaryCondition::Neumann,
        )?;

        operator.apply(&input, &mut output)?;

        assert_eq!(
            output.iter().copied().collect::<Vec<_>>(),
            vec![-8.0; nx * ny]
        );
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
