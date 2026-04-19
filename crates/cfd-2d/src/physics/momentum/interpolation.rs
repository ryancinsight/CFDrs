//! Interpolation schemes including Rhie-Chow.
//!
//! # Theorem
//! The momentum discretization must conserve linear momentum globally and locally.
//!
//! **Proof sketch**:
//! By integrating the Navier-Stokes momentum equation over a control volume $\Omega$,
//! Gauss's divergence theorem converts the convective and diffusive volume integrals
//! into surface fluxes. The finite volume method ensures that the flux leaving one
//! cell exactly equals the flux entering the adjacent cell. Thus, in the absence of
//! external forces and boundary fluxes, the total momentum $\int_\Omega \rho \mathbf{u} dV$
//! is exactly conserved to machine precision.

use crate::fields::{Field2D, SimulationFields};
use crate::grid::array2d::Array2D;
use nalgebra::RealField;
use num_traits::FromPrimitive;

fn harmonic_face_coefficient<T: RealField + Copy + FromPrimitive>(
    ap_left: T,
    ap_right: T,
    volume: T,
) -> T {
    let tiny = T::from_f64(1e-30).unwrap_or_else(T::default_epsilon);
    let ap_left = if ap_left > tiny { ap_left } else { tiny };
    let ap_right = if ap_right > tiny { ap_right } else { tiny };
    let d_left = volume / ap_left;
    let d_right = volume / ap_right;
    let denom = d_left + d_right;
    if denom <= tiny {
        T::zero()
    } else {
        (T::one() + T::one()) * d_left * d_right / denom
    }
}

/// Rhie-Chow interpolation for face velocities using full momentum diagonals.
///
/// The face coefficient uses harmonic averaging of the adjacent momentum diagonals,
/// matching the coefficient-aware Rhie-Chow formulation used by the pressure-velocity
/// coupling solvers elsewhere in the crate.
pub fn rhie_chow_interpolation<T: RealField + Copy + FromPrimitive>(
    fields: &SimulationFields<T>,
    ap_u: &Field2D<T>,
    ap_v: &Field2D<T>,
    dx: T,
    dy: T,
) -> (Array2D<T>, Array2D<T>) {
    let (nx, ny) = fields.u.dimensions();
    let two = T::one() + T::one();
    let volume = dx * dy;

    let mut u_face = Array2D::new(nx + 1, ny, T::zero());
    let mut v_face = Array2D::new(nx, ny + 1, T::zero());

    for j in 0..ny {
        for i in 0..=nx {
            if i == 0 {
                u_face[(i, j)] = fields.u.at(0, j);
            } else if i == nx {
                u_face[(i, j)] = fields.u.at(nx - 1, j);
            } else {
                let left = i - 1;
                let right = i;
                let u_bar = (fields.u.at(left, j) + fields.u.at(right, j)) / two;
                let d_face = harmonic_face_coefficient(ap_u.at(left, j), ap_u.at(right, j), volume);

                let dp_dx_left = if left > 0 && left < nx - 1 {
                    (fields.p.at(left + 1, j) - fields.p.at(left - 1, j)) / (two * dx)
                } else if left == 0 {
                    (fields.p.at(left + 1, j) - fields.p.at(left, j)) / dx
                } else {
                    (fields.p.at(left, j) - fields.p.at(left - 1, j)) / dx
                };
                let dp_dx_right = if right > 0 && right < nx - 1 {
                    (fields.p.at(right + 1, j) - fields.p.at(right - 1, j)) / (two * dx)
                } else if right == nx - 1 {
                    (fields.p.at(right, j) - fields.p.at(right - 1, j)) / dx
                } else {
                    (fields.p.at(right + 1, j) - fields.p.at(right, j)) / dx
                };
                let dp_dx_cells = (dp_dx_left + dp_dx_right) / two;
                let dp_dx_face = (fields.p.at(right, j) - fields.p.at(left, j)) / dx;

                u_face[(i, j)] = u_bar + d_face * (dp_dx_cells - dp_dx_face);
            }
        }
    }

    for i in 0..nx {
        for j in 0..=ny {
            if j == 0 {
                v_face[(i, j)] = fields.v.at(i, 0);
            } else if j == ny {
                v_face[(i, j)] = fields.v.at(i, ny - 1);
            } else {
                let bottom = j - 1;
                let top = j;
                let v_bar = (fields.v.at(i, bottom) + fields.v.at(i, top)) / two;
                let d_face = harmonic_face_coefficient(ap_v.at(i, bottom), ap_v.at(i, top), volume);

                let dp_dy_bottom = if bottom > 0 && bottom < ny - 1 {
                    (fields.p.at(i, bottom + 1) - fields.p.at(i, bottom - 1)) / (two * dy)
                } else if bottom == 0 {
                    (fields.p.at(i, bottom + 1) - fields.p.at(i, bottom)) / dy
                } else {
                    (fields.p.at(i, bottom) - fields.p.at(i, bottom - 1)) / dy
                };
                let dp_dy_top = if top > 0 && top < ny - 1 {
                    (fields.p.at(i, top + 1) - fields.p.at(i, top - 1)) / (two * dy)
                } else if top == ny - 1 {
                    (fields.p.at(i, top) - fields.p.at(i, top - 1)) / dy
                } else {
                    (fields.p.at(i, top + 1) - fields.p.at(i, top)) / dy
                };
                let dp_dy_cells = (dp_dy_bottom + dp_dy_top) / two;
                let dp_dy_face = (fields.p.at(i, top) - fields.p.at(i, bottom)) / dy;

                v_face[(i, j)] = v_bar + d_face * (dp_dy_cells - dp_dy_face);
            }
        }
    }

    (u_face, v_face)
}

#[cfg(test)]
mod tests {
    use super::rhie_chow_interpolation;
    use crate::fields::{Field2D, SimulationFields};
    use crate::grid::StructuredGrid2D;
    use crate::pressure_velocity::RhieChowInterpolation;
    use approx::assert_relative_eq;
    use nalgebra::Vector2;

    #[test]
    fn coefficient_aware_interpolation_matches_exact_reference() {
        let grid = StructuredGrid2D::new(4, 4, 0.0, 1.0, 0.0, 1.0).unwrap();
        let mut fields = SimulationFields::new(4, 4);
        let mut u_vec = Field2D::new(4, 4, Vector2::new(0.0, 0.0));
        let ap_u = Field2D::new(4, 4, 3.0_f64);
        let ap_v = Field2D::new(4, 4, 5.0_f64);

        for i in 0..4 {
            for j in 0..4 {
                let x = i as f64;
                let y = j as f64;
                fields.u.set(i, j, 1.0 + 0.2 * x - 0.1 * y);
                fields.v.set(i, j, -0.5 + 0.3 * x + 0.15 * y);
                fields.p.set(i, j, x * x * x + 0.25 * y * y);
                u_vec.set(i, j, Vector2::new(fields.u.at(i, j), fields.v.at(i, j)));
            }
        }

        let (u_face, v_face) = rhie_chow_interpolation(&fields, &ap_u, &ap_v, grid.dx, grid.dy);

        let mut exact = RhieChowInterpolation::new(&grid);
        exact.update_u_coefficients(&ap_u);
        exact.update_v_coefficients(&ap_v);

        for j in 0..grid.ny {
            for i in 1..grid.nx {
                let reference =
                    exact.face_velocity_x(&u_vec, &fields.p, grid.dx, grid.dy, None, i - 1, j);
                assert_relative_eq!(u_face[(i, j)], reference, epsilon = 1e-12);
            }
        }

        for i in 0..grid.nx {
            for j in 1..grid.ny {
                let reference =
                    exact.face_velocity_y(&u_vec, &fields.p, grid.dx, grid.dy, None, i, j - 1);
                assert_relative_eq!(v_face[(i, j)], reference, epsilon = 1e-12);
            }
        }

        assert_relative_eq!(u_face[(0, 1)], fields.u.at(0, 1), epsilon = 1e-12);
        assert_relative_eq!(u_face[(4, 1)], fields.u.at(3, 1), epsilon = 1e-12);
        assert_relative_eq!(v_face[(1, 0)], fields.v.at(1, 0), epsilon = 1e-12);
        assert_relative_eq!(v_face[(1, 4)], fields.v.at(1, 3), epsilon = 1e-12);
    }
}
