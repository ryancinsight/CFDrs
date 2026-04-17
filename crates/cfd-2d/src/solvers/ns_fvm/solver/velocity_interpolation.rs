//! Staggered-grid velocity interpolation for the SIMPLE solver.
//!
//! # Algorithm
//! The SIMPLE solver stores `u` on vertical faces and `v` on horizontal faces.
//! To evaluate velocity at an arbitrary point `(x, y)`, this module performs
//! axis-aligned linear interpolation on the corresponding staggered lattice.
//! Uniform x coordinates use direct index arithmetic; stretched y grids use a
//! monotone bracket search over the physical face/center coordinates.
//!
//! # Theorem - Linear Reproduction
//! For any affine field `f(x, y) = a + b x + c y`, the interpolated value on
//! the staggered lattice equals the exact analytic value at every point inside
//! the convex hull of the sampled coordinates.
//!
//! **Proof sketch**: 1D linear interpolation reproduces affine functions exactly.
//! Bilinear interpolation is the tensor product of two exact 1D interpolants, so
//! an affine field remains unchanged under interpolation. The staggered layout
//! changes sample locations, not the interpolation degree.

use super::NavierStokesSolver2D;
use crate::solvers::cell_tracking::physics::VelocityFieldInterpolator;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

impl<T: RealField + Copy + Float + FromPrimitive> NavierStokesSolver2D<T> {
    #[inline]
    fn uniform_bracket(query: T, origin: T, step: T, count: usize) -> (usize, T) {
        let zero = T::zero();
        let one = T::one();

        if count <= 1 || step <= zero {
            return (0, zero);
        }

        let last = origin + step * T::from_usize(count - 1).expect("grid size fits in T");
        if query <= origin {
            return (0, zero);
        }
        if query >= last {
            return (count - 2, one);
        }

        let scaled = (query - origin) / step;
        let mut index: usize = Float::floor(scaled).to_usize().unwrap_or(0);
        if index >= count - 1 {
            index = count - 2;
        }

        let left = origin + step * T::from_usize(index).expect("grid index fits in T");
        let mut fraction = (query - left) / step;
        if fraction < zero {
            fraction = zero;
        }
        if fraction > one {
            fraction = one;
        }
        (index, fraction)
    }

    #[inline]
    fn monotonic_bracket<F>(query: T, count: usize, coordinate: F) -> (usize, T)
    where
        F: Fn(usize) -> T,
    {
        let zero = T::zero();
        let one = T::one();

        if count <= 1 {
            return (0, zero);
        }

        let first = coordinate(0);
        if query <= first {
            return (0, zero);
        }

        let last_index = count - 1;
        let last = coordinate(last_index);
        if query >= last {
            return (count - 2, one);
        }

        let mut low = 0usize;
        let mut high = last_index;
        while high - low > 1 {
            let mid = low + (high - low) / 2;
            if coordinate(mid) <= query {
                low = mid;
            } else {
                high = mid;
            }
        }

        let left = coordinate(low);
        let right = coordinate(low + 1);
        let denom = right - left;
        if Float::abs(denom) <= T::from_f64(1.0e-30).expect("tolerance fits in T") {
            (low, zero)
        } else {
            let mut fraction = (query - left) / denom;
            if fraction < zero {
                fraction = zero;
            }
            if fraction > one {
                fraction = one;
            }
            (low, fraction)
        }
    }

    #[inline]
    fn interpolate_u_component(&self, x: T, y: T) -> T {
        let x_count = self.field.u.rows();
        let y_count = self.field.u.cols();
        let zero = T::zero();
        let one = T::one();
        let half = T::from_f64(0.5).expect("0.5 fits in T");

        if x_count == 0 || y_count == 0 {
            return zero;
        }
        if x_count == 1 && y_count == 1 {
            return self.field.u[(0, 0)];
        }

        if x_count == 1 {
            let (j, ty) = if self.grid.y_faces.is_some() {
                Self::monotonic_bracket(y, y_count, |jj| self.grid.y_center(jj))
            } else {
                Self::uniform_bracket(y, self.grid.dy * half, self.grid.dy, y_count)
            };

            let u0 = self.field.u[(0, j)];
            if j + 1 >= y_count {
                return u0;
            }
            let u1 = self.field.u[(0, j + 1)];
            return u0 * (one - ty) + u1 * ty;
        }

        let (i, tx) = Self::uniform_bracket(x, T::zero(), self.grid.dx, x_count);

        if y_count == 1 {
            let u0 = self.field.u[(i, 0)];
            let u1 = self.field.u[(i + 1, 0)];
            return u0 * (one - tx) + u1 * tx;
        }

        let (j, ty) = if self.grid.y_faces.is_some() {
            Self::monotonic_bracket(y, y_count, |jj| self.grid.y_center(jj))
        } else {
            Self::uniform_bracket(y, self.grid.dy * half, self.grid.dy, y_count)
        };

        let u00 = self.field.u[(i, j)];
        let u10 = self.field.u[(i + 1, j)];
        let u01 = self.field.u[(i, j + 1)];
        let u11 = self.field.u[(i + 1, j + 1)];
        let u_bottom = u00 * (one - tx) + u10 * tx;
        let u_top = u01 * (one - tx) + u11 * tx;
        u_bottom * (one - ty) + u_top * ty
    }

    #[inline]
    fn interpolate_v_component(&self, x: T, y: T) -> T {
        let x_count = self.field.v.rows();
        let y_count = self.field.v.cols();
        let zero = T::zero();
        let one = T::one();
        let half = T::from_f64(0.5).expect("0.5 fits in T");

        if x_count == 0 || y_count == 0 {
            return zero;
        }
        if x_count == 1 && y_count == 1 {
            return self.field.v[(0, 0)];
        }

        if x_count == 1 {
            let (j, ty) = if self.grid.y_faces.is_some() {
                Self::monotonic_bracket(y, y_count, |jj| self.grid.y_v_face(jj))
            } else {
                Self::uniform_bracket(y, T::zero(), self.grid.dy, y_count)
            };

            let v0 = self.field.v[(0, j)];
            if j + 1 >= y_count {
                return v0;
            }
            let v1 = self.field.v[(0, j + 1)];
            return v0 * (one - ty) + v1 * ty;
        }

        let (i, tx) = Self::uniform_bracket(x, self.grid.dx * half, self.grid.dx, x_count);

        if y_count == 1 {
            let v0 = self.field.v[(i, 0)];
            let v1 = self.field.v[(i + 1, 0)];
            return v0 * (one - tx) + v1 * tx;
        }

        let (j, ty) = if self.grid.y_faces.is_some() {
            Self::monotonic_bracket(y, y_count, |jj| self.grid.y_v_face(jj))
        } else {
            Self::uniform_bracket(y, T::zero(), self.grid.dy, y_count)
        };

        let v00 = self.field.v[(i, j)];
        let v10 = self.field.v[(i + 1, j)];
        let v01 = self.field.v[(i, j + 1)];
        let v11 = self.field.v[(i + 1, j + 1)];
        let v_bottom = v00 * (one - tx) + v10 * tx;
        let v_top = v01 * (one - tx) + v11 * tx;
        v_bottom * (one - ty) + v_top * ty
    }
}

impl<T: RealField + Copy + Float + FromPrimitive> VelocityFieldInterpolator for NavierStokesSolver2D<T> {
    fn velocity_at(&self, x: f64, y: f64) -> (f64, f64) {
        let x_t = T::from_f64(x).unwrap_or(T::zero());
        let y_t = T::from_f64(y).unwrap_or(T::zero());

        let u = self.interpolate_u_component(x_t, y_t);
        let v = self.interpolate_v_component(x_t, y_t);

        (u.to_f64().unwrap_or(0.0), v.to_f64().unwrap_or(0.0))
    }

    fn is_fluid(&self, x: f64, y: f64) -> bool {
        let lx_f64 = self.grid.lx.to_f64().unwrap_or(0.0);
        let ly_f64 = self.grid.ly.to_f64().unwrap_or(0.0);
        x >= 0.0 && x <= lx_f64 && y >= 0.0 && y <= ly_f64
    }

    fn bounds(&self) -> (f64, f64, f64, f64) {
        let lx_f64 = self.grid.lx.to_f64().unwrap_or(0.0);
        let ly_f64 = self.grid.ly.to_f64().unwrap_or(0.0);
        (0.0, lx_f64, 0.0, ly_f64)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use cfd_core::physics::fluid::BloodModel;
    use cfd_math::pressure_velocity::SIMPLEConfig;
    use crate::solvers::ns_fvm::StaggeredGrid2D;

    #[test]
    fn velocity_at_reproduces_affine_staggered_field() {
        let grid = StaggeredGrid2D::<f64>::new_stretched_y(4, 3, 2.0, vec![0.0, 0.15, 0.55, 1.0]);
        let blood = BloodModel::Newtonian(1.0e-3);
        let config = SIMPLEConfig::new(1, 1.0e-6, 0.7, 0.3, 1.0, 1, 1);
        let mut solver = NavierStokesSolver2D::new(grid.clone(), blood, 1000.0, config);

        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                let x = grid.x_u_face(i);
                let y = grid.y_center(j);
                solver.field.u[(i, j)] = 1.25 + 2.0 * x + 3.0 * y;
            }
        }
        for i in 0..grid.nx {
            for j in 0..=grid.ny {
                let x = grid.x_center(i);
                let y = grid.y_v_face(j);
                solver.field.v[(i, j)] = -0.75 - x + 4.0 * y;
            }
        }

        let x = 0.73;
        let y = 0.41;
        let (u, v) = solver.velocity_at(x, y);

        assert_relative_eq!(u, 1.25 + 2.0 * x + 3.0 * y, epsilon = 1e-12);
        assert_relative_eq!(v, -0.75 - x + 4.0 * y, epsilon = 1e-12);
    }

    #[test]
    fn velocity_at_clamps_outside_domain_to_edge_samples() {
        let grid = StaggeredGrid2D::<f64>::new(4, 4, 2.0, 1.0);
        let blood = BloodModel::Newtonian(1.0e-3);
        let config = SIMPLEConfig::default();
        let mut solver = NavierStokesSolver2D::new(grid.clone(), blood, 1000.0, config);

        for i in 0..=grid.nx {
            for j in 0..grid.ny {
                let x = grid.x_u_face(i);
                let y = grid.y_center(j);
                solver.field.u[(i, j)] = x + y;
            }
        }
        for i in 0..grid.nx {
            for j in 0..=grid.ny {
                let x = grid.x_center(i);
                let y = grid.y_v_face(j);
                solver.field.v[(i, j)] = x - y;
            }
        }

        let (u, v) = solver.velocity_at(-0.2, 1.4);
        assert!(u.is_finite() && v.is_finite());
        let (u_edge, v_edge) = solver.velocity_at(0.0, 1.0);
        assert_relative_eq!(u, u_edge, epsilon = 1e-12);
        assert_relative_eq!(v, v_edge, epsilon = 1e-12);
    }
}
