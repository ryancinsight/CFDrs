//! Interpolation schemes including Rhie-Chow

use crate::fields::SimulationFields;
use nalgebra::RealField;

/// Rhie-Chow interpolation for face velocities
/// Prevents checkerboard pressure oscillations
pub fn rhie_chow_interpolation<T: RealField + Copy>(
    fields: &SimulationFields<T>,
    dx: T,
    dy: T,
    dt: T,
) -> (Vec<Vec<T>>, Vec<Vec<T>>) {
    let ny = fields.u.rows();
    let nx = fields.u.cols();

    // Face velocities
    let mut u_face = vec![vec![T::zero(); nx + 1]; ny];
    let mut v_face = vec![vec![T::zero(); nx]; ny + 1];

    // U-velocity at vertical faces
    for j in 0..ny {
        for i in 0..=nx {
            if i == 0 {
                // West boundary
                u_face[j][i] = fields.u.at(0, j);
            } else if i == nx {
                // East boundary
                u_face[j][i] = fields.u.at(nx - 1, j);
            } else {
                // Interior faces - linear interpolation with pressure correction
                let u_avg = (fields.u.at(i - 1, j) + fields.u.at(i, j)) / (T::one() + T::one());

                // Rhie-Chow pressure correction term
                let dp_dx = (fields.p.at(i, j) - fields.p.at(i - 1, j)) / dx;
                let correction = dt * dp_dx / fields.density.at(i, j);

                u_face[j][i] = u_avg - correction;
            }
        }
    }

    // V-velocity at horizontal faces
    for j in 0..=ny {
        for i in 0..nx {
            if j == 0 {
                // South boundary
                v_face[j][i] = fields.v.at(i, 0);
            } else if j == ny {
                // North boundary
                v_face[j][i] = fields.v.at(i, ny - 1);
            } else {
                // Interior faces
                let v_avg = (fields.v.at(i, j - 1) + fields.v.at(i, j)) / (T::one() + T::one());

                // Rhie-Chow pressure correction
                let dp_dy = (fields.p.at(i, j) - fields.p.at(i, j - 1)) / dy;
                let correction = dt * dp_dy / fields.density.at(i, j);

                v_face[j][i] = v_avg - correction;
            }
        }
    }

    (u_face, v_face)
}
