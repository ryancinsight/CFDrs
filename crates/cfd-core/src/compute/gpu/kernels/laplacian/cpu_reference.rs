//! CPU reference implementation of the 2D Laplacian operator.
//!
//! Provides a numerically deterministic fallback for validation and small grids.

use super::types::BoundaryType;

/// CPU reference implementation of the 2D Laplacian with full boundary condition support.
///
/// Computes in f64 intermediate precision and casts back to f32 for accuracy.
pub(super) fn execute_cpu_reference(
    field: &[f32],
    nx: usize,
    ny: usize,
    dx: f32,
    dy: f32,
    bc: BoundaryType,
    result: &mut [f32],
) {
    let dx_inv2_64 = 1.0f64 / (f64::from(dx) * f64::from(dx));
    let dy_inv2_64 = 1.0f64 / (f64::from(dy) * f64::from(dy));

    for y in 0..ny {
        for x in 0..nx {
            let idx = y * nx + x;
            let mut laplacian = 0.0f64;

            if x > 0 && x < nx - 1 {
                let left = field[y * nx + (x - 1)];
                let center = field[idx];
                let right = field[y * nx + (x + 1)];
                laplacian +=
                    (f64::from(left) - 2.0f64 * f64::from(center) + f64::from(right)) * dx_inv2_64;
            } else if x == 0 {
                match bc {
                    BoundaryType::Dirichlet => {
                        let center = field[idx];
                        laplacian += (-2.0f64 * f64::from(center)) * dx_inv2_64;
                    }
                    BoundaryType::Neumann => {
                        let center = field[idx];
                        if nx >= 4 {
                            let u0 = center;
                            let u1 = field[y * nx + 1];
                            let u2 = field[y * nx + 2];
                            let u3 = field[y * nx + 3];
                            let d2x = f64::from(2.0 * u0 - 5.0 * u1 + 4.0 * u2 - u3) * dx_inv2_64;
                            laplacian += d2x;
                        } else {
                            let right = field[y * nx + (x + 1)];
                            laplacian += (f64::from(right) - 2.0f64 * f64::from(center)
                                + f64::from(right))
                                * dx_inv2_64;
                        }
                    }
                    BoundaryType::Periodic => {
                        let left = field[y * nx + (nx - 2)];
                        let center = field[idx];
                        let right = field[y * nx + (x + 1)];
                        laplacian += (f64::from(left) - 2.0f64 * f64::from(center)
                            + f64::from(right))
                            * dx_inv2_64;
                    }
                }
            } else if x == nx - 1 {
                match bc {
                    BoundaryType::Dirichlet => {
                        let center = field[idx];
                        laplacian += (-2.0f64 * f64::from(center)) * dx_inv2_64;
                    }
                    BoundaryType::Neumann => {
                        let center = field[idx];
                        if nx >= 4 {
                            let u0 = center;
                            let u1 = field[y * nx + (nx - 2)];
                            let u2 = field[y * nx + (nx - 3)];
                            let u3 = field[y * nx + (nx - 4)];
                            let d2x = f64::from(2.0 * u0 - 5.0 * u1 + 4.0 * u2 - u3) * dx_inv2_64;
                            laplacian += d2x;
                        } else {
                            let left = field[y * nx + (x - 1)];
                            laplacian += (f64::from(left) - 2.0f64 * f64::from(center)
                                + f64::from(left))
                                * dx_inv2_64;
                        }
                    }
                    BoundaryType::Periodic => {
                        let left = field[y * nx + (x - 1)];
                        let center = field[idx];
                        let right = field[y * nx + 1];
                        laplacian += (f64::from(left) - 2.0f64 * f64::from(center)
                            + f64::from(right))
                            * dx_inv2_64;
                    }
                }
            }

            if y > 0 && y < ny - 1 {
                let bottom = field[(y - 1) * nx + x];
                let center = field[idx];
                let top = field[(y + 1) * nx + x];
                laplacian +=
                    (f64::from(bottom) - 2.0f64 * f64::from(center) + f64::from(top)) * dy_inv2_64;
            } else if y == 0 {
                match bc {
                    BoundaryType::Dirichlet => {
                        let center = field[idx];
                        laplacian += (-2.0f64 * f64::from(center)) * dy_inv2_64;
                    }
                    BoundaryType::Neumann => {
                        let center = field[idx];
                        if ny >= 4 {
                            let u0 = center;
                            let u1 = field[nx + x];
                            let u2 = field[2 * nx + x];
                            let u3 = field[3 * nx + x];
                            let d2y = f64::from(2.0 * u0 - 5.0 * u1 + 4.0 * u2 - u3) * dy_inv2_64;
                            laplacian += d2y;
                        } else {
                            let top = field[(y + 1) * nx + x];
                            laplacian += (f64::from(top) - 2.0f64 * f64::from(center)
                                + f64::from(top))
                                * dy_inv2_64;
                        }
                    }
                    BoundaryType::Periodic => {
                        let bottom = field[(ny - 2) * nx + x];
                        let center = field[idx];
                        let top = field[(y + 1) * nx + x];
                        laplacian += (f64::from(bottom) - 2.0f64 * f64::from(center)
                            + f64::from(top))
                            * dy_inv2_64;
                    }
                }
            } else if y == ny - 1 {
                match bc {
                    BoundaryType::Dirichlet => {
                        let center = field[idx];
                        laplacian += (-2.0f64 * f64::from(center)) * dy_inv2_64;
                    }
                    BoundaryType::Neumann => {
                        let center = field[idx];
                        if ny >= 4 {
                            let u0 = center;
                            let u1 = field[(ny - 2) * nx + x];
                            let u2 = field[(ny - 3) * nx + x];
                            let u3 = field[(ny - 4) * nx + x];
                            let d2y = f64::from(2.0 * u0 - 5.0 * u1 + 4.0 * u2 - u3) * dy_inv2_64;
                            laplacian += d2y;
                        } else {
                            let bottom = field[(y - 1) * nx + x];
                            laplacian += (f64::from(bottom) - 2.0f64 * f64::from(center)
                                + f64::from(bottom))
                                * dy_inv2_64;
                        }
                    }
                    BoundaryType::Periodic => {
                        let bottom = field[(y - 1) * nx + x];
                        let center = field[idx];
                        let top = field[nx + x];
                        laplacian += (f64::from(bottom) - 2.0f64 * f64::from(center)
                            + f64::from(top))
                            * dy_inv2_64;
                    }
                }
            }

            result[idx] = laplacian as f32;
        }
    }
}
