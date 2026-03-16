//! Pressure correction Poisson solver.

use crate::error::Error;
use crate::solvers::ns_fvm::solver::NavierStokesSolver2D;
use crate::solvers::ns_fvm::array2d::Array2D;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

impl<T: RealField + Copy + Float + FromPrimitive> NavierStokesSolver2D<T> {
    /// Solves the pressure-correction Poisson equation.
    ///
    /// Formulates the continuity residual from the intermediate velocity fields 
    /// and solves for the pressure correction `p'` iteratively via SOR.
    pub fn solve_pressure_correction(&mut self) -> Result<(), Error> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let rho = self.density;
        let zero = T::zero();
        let tiny = T::from_f64(1e-30).unwrap_or(zero);

        let mut d_u = Array2D::new(nx + 1, ny, zero);
        let mut d_v = Array2D::new(nx, ny + 1, zero);

        for i in 1..nx {
            for j in 0..ny {
                let a = self.a_p_u[(i, j)];
                if a > tiny {
                    d_u[(i, j)] = self.grid.dy_at(j) / a;
                }
            }
        }
        for j in 0..ny {
            d_u[(nx, j)] = d_u[(nx - 1, j)];
        }
        for i in 0..nx {
            for j in 1..ny {
                let a = self.a_p_v[(i, j)];
                if a > tiny {
                    d_v[(i, j)] = dx / a;
                }
            }
        }
        if ny > 0 {
            for i in 0..nx {
                d_v[(i, ny)] = d_v[(i, ny - 1)];
            }
        }

        let mut p_prime = Array2D::new(nx, ny, zero);
        let mut b = Array2D::new(nx, ny, zero);

        for i in 0..nx {
            for j in 0..ny {
                if !self.field.mask[(i, j)] {
                    continue;
                }
                let dy_j = self.grid.dy_at(j);
                b[(i, j)] = rho
                    * ((self.field.u[(i, j)] - self.field.u[(i + 1, j)]) * dy_j
                        + (self.field.v[(i, j)] - self.field.v[(i, j + 1)]) * dx);
            }
        }

        for _ in 0..200 {
            for i in 0..nx {
                for j in 0..ny {
                    if !self.field.mask[(i, j)] {
                        continue;
                    }
                    let dy_j = self.grid.dy_at(j);
                    let a_e = if i + 1 < nx {
                        if self.field.mask[(i + 1, j)] {
                            rho * d_u[(i + 1, j)] * dy_j
                        } else {
                            zero
                        }
                    } else {
                        rho * d_u[(i + 1, j)] * dy_j
                    };
                    let a_w = if i > 0 {
                        if self.field.mask[(i - 1, j)] {
                            rho * d_u[(i, j)] * dy_j
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    let a_n = if j + 1 < ny {
                        if self.field.mask[(i, j + 1)] {
                            rho * d_v[(i, j + 1)] * dx
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    let a_s = if j > 0 {
                        if self.field.mask[(i, j - 1)] {
                            rho * d_v[(i, j)] * dx
                        } else {
                            zero
                        }
                    } else {
                        zero
                    };
                    let a_p = a_e + a_w + a_n + a_s;
                    if a_p < tiny {
                        continue;
                    }
                    let pe = if i + 1 < nx { p_prime[(i + 1, j)] } else { zero };
                    let pw = if i > 0 {
                        p_prime[(i - 1, j)]
                    } else {
                        p_prime[(i, j)]
                    };
                    let pn = if j + 1 < ny {
                        p_prime[(i, j + 1)]
                    } else {
                        p_prime[(i, j)]
                    };
                    let ps = if j > 0 {
                        p_prime[(i, j - 1)]
                    } else {
                        p_prime[(i, j)]
                    };
                    p_prime[(i, j)] = (a_e * pe + a_w * pw + a_n * pn + a_s * ps + b[(i, j)]) / a_p;
                }
            }
        }

        for i in 1..=nx {
            for j in 0..ny {
                if i < nx {
                    if !self.field.mask[(i, j)] && !self.field.mask[(i - 1, j)] {
                        continue;
                    }
                } else if !self.field.mask[(nx - 1, j)] {
                    continue;
                }
                let dp = if i > 0 && i < nx {
                    p_prime[(i - 1, j)] - p_prime[(i, j)]
                } else if i == nx {
                    p_prime[(i - 1, j)]
                } else {
                    zero
                };
                self.field.u[(i, j)] += d_u[(i, j)] * dp;
            }
        }

        for i in 0..nx {
            for j in 1..=ny {
                if j < ny {
                    if !self.field.mask[(i, j)] && !self.field.mask[(i, j - 1)] {
                        continue;
                    }
                } else if !self.field.mask[(i, ny - 1)] {
                    continue;
                }
                let dp = if j > 0 && j < ny {
                    p_prime[(i, j - 1)] - p_prime[(i, j)]
                } else {
                    zero
                };
                self.field.v[(i, j)] += d_v[(i, j)] * dp;
            }
        }

        let alpha_p = self.config.alpha_p;
        for i in 0..nx {
            for j in 0..ny {
                self.field.p[(i, j)] += alpha_p * p_prime[(i, j)];
            }
        }

        Ok(())
    }
}
