//! Pressure corrector step for PISO algorithm

use crate::fields::{Field2D, SimulationFields};
use crate::grid::StructuredGrid2D;
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_core::error::Result;
use eunomia::{FloatElement, NumericElement};
use leto::geometry::Vector2;

// Named constants
const ONE: f64 = 1.0;

/// Pressure corrector for PISO algorithm
pub struct PressureCorrector<T: Cfd2dScalar + Copy> {
    /// Grid dimensions
    nx: usize,
    ny: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    /// Number of corrector steps
    num_correctors: usize,
    /// Under-relaxation factor for pressure
    pressure_relaxation: T,
}

impl<T: Cfd2dScalar + Copy + FloatElement> PressureCorrector<T> {
    /// Create new pressure corrector
    pub fn new(grid: &StructuredGrid2D<T>, num_correctors: usize, pressure_relaxation: T) -> Self {
        Self {
            nx: grid.nx,
            ny: grid.ny,
            dx: grid.dx,
            dy: grid.dy,
            num_correctors,
            pressure_relaxation,
        }
    }

    /// Perform pressure correction steps
    pub fn correct(&self, fields: &mut SimulationFields<T>, dt: T) -> Result<()> {
        for corrector in 0..self.num_correctors {
            // Solve pressure correction equation
            let p_prime = self.solve_pressure_correction(fields, dt)?;

            // Correct pressure field
            self.correct_pressure(fields, &p_prime);

            // Correct velocity field
            self.correct_velocity(fields, &p_prime, dt);

            // Update face fluxes for next corrector (if any)
            if corrector < self.num_correctors - 1 {
                self.update_face_fluxes(fields, dt);
            }
        }

        Ok(())
    }

    /// Solve pressure correction equation according to Issa (1986)
    /// Reference: Issa, R.I. (1986). "Solution of the implicitly discretised fluid flow equations by operator-splitting"
    /// Journal of Computational Physics, 62(1), 40-65.
    fn solve_pressure_correction(&self, fields: &SimulationFields<T>, dt: T) -> Result<Field2D<T>> {
        let mut p_prime = Field2D::new(self.nx, self.ny, scalar::zero::<T>());
        let mut residual = scalar::from_f64::<T>(ONE);
        let tolerance = scalar::from_f64(
            cfd_core::physics::constants::numerical::solver::CONVERGENCE_TOLERANCE,
        );
        let max_iter = cfd_core::physics::constants::numerical::solver::MAX_ITERATIONS_OUTER;
        let mut iter = 0;

        // Calculate H(u) operator for PISO neighbor correction
        let h_operator = self.calculate_h_operator(fields);

        while residual > tolerance && iter < max_iter {
            residual = scalar::zero::<T>();

            // Gauss-Seidel iteration with H(u) correction
            for i in 1..self.nx - 1 {
                for j in 1..self.ny - 1 {
                    // Calculate mass imbalance including H(u) terms
                    let mass_imbalance =
                        self.calculate_mass_imbalance_with_h(fields, &h_operator, i, j);

                    // Calculate coefficients
                    let ae = fields.density.at(i, j) * self.dy * dt / self.dx;
                    let aw = fields.density.at(i, j) * self.dy * dt / self.dx;
                    let an = fields.density.at(i, j) * self.dx * dt / self.dy;
                    let as_ = fields.density.at(i, j) * self.dx * dt / self.dy;
                    let ap = ae + aw + an + as_;

                    // Update pressure correction
                    let p_current = p_prime.at(i, j);
                    let p_updated = (ae * p_prime.at(i + 1, j)
                        + aw * p_prime.at(i - 1, j)
                        + an * p_prime.at(i, j + 1)
                        + as_ * p_prime.at(i, j - 1)
                        - mass_imbalance)
                        / ap;

                    if let Some(p) = p_prime.at_mut(i, j) {
                        *p = p_updated;
                    }

                    // Calculate residual
                    let diff = <T as NumericElement>::abs(p_updated - p_current);
                    if diff > residual {
                        residual = diff;
                    }
                }
            }

            iter += 1;
        }

        Ok(p_prime)
    }

    /// Calculate H(u) operator for PISO neighbor correction
    /// H(u) = -`sum(A_nb` * `u_nb`) where `A_nb` are momentum equation coefficients
    fn calculate_h_operator(&self, fields: &SimulationFields<T>) -> Field2D<Vector2<T>> {
        let mut h_field = Field2D::new(self.nx, self.ny, Vector2::zeros());

        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                // Calculate neighbor contributions to H(u)
                let u_e = Vector2::new(fields.u.at(i + 1, j), fields.v.at(i + 1, j));
                let u_w = Vector2::new(fields.u.at(i - 1, j), fields.v.at(i - 1, j));
                let u_n = Vector2::new(fields.u.at(i, j + 1), fields.v.at(i, j + 1));
                let u_s = Vector2::new(fields.u.at(i, j - 1), fields.v.at(i, j - 1));

                // Momentum equation coefficients with proper discretization
                // Include both diffusion and convection terms
                let visc = fields.viscosity.at(i, j);
                let density = fields.density.at(i, j);

                let two_t = scalar::one::<T>() + scalar::one::<T>();

                // Face velocities (using upwind values)
                let u_e_face = (fields.u.at(i, j) + fields.u.at(i + 1, j)) / two_t;
                let u_w_face = (fields.u.at(i - 1, j) + fields.u.at(i, j)) / two_t;
                let v_n_face = (fields.v.at(i, j) + fields.v.at(i, j + 1)) / two_t;
                let v_s_face = (fields.v.at(i, j - 1) + fields.v.at(i, j)) / two_t;

                // Convective fluxes
                let f_e = density * u_e_face * self.dy;
                let f_w = density * u_w_face * self.dy;
                let f_n = density * v_n_face * self.dx;
                let f_s = density * v_s_face * self.dx;

                // Diffusion conductances
                let d_e = visc * self.dy / self.dx;
                let d_w = visc * self.dy / self.dx;
                let d_n = visc * self.dx / self.dy;
                let d_s = visc * self.dx / self.dy;

                // Hybrid differencing scheme coefficients
                let half = scalar::one::<T>() / (scalar::one::<T>() + scalar::one::<T>());
                let zero = scalar::zero::<T>();
                let one = scalar::one::<T>();
                let ae =
                    d_e * <T as NumericElement>::max_scalar(
                        zero,
                        one - half * <T as NumericElement>::abs(f_e) / d_e,
                    ) + <T as NumericElement>::max_scalar(-f_e, zero);
                let aw =
                    d_w * <T as NumericElement>::max_scalar(
                        zero,
                        one - half * <T as NumericElement>::abs(f_w) / d_w,
                    ) + <T as NumericElement>::max_scalar(f_w, zero);
                let an =
                    d_n * <T as NumericElement>::max_scalar(
                        zero,
                        one - half * <T as NumericElement>::abs(f_n) / d_n,
                    ) + <T as NumericElement>::max_scalar(-f_n, zero);
                let as_ =
                    d_s * <T as NumericElement>::max_scalar(
                        zero,
                        one - half * <T as NumericElement>::abs(f_s) / d_s,
                    ) + <T as NumericElement>::max_scalar(f_s, zero);

                // H(u) = -sum(A_nb * u_nb)
                let h_u = -(u_e * ae + u_w * aw + u_n * an + u_s * as_);

                if let Some(h) = h_field.at_mut(i, j) {
                    *h = h_u;
                }
            }
        }

        h_field
    }

    /// Calculate mass imbalance including H(u) correction terms
    fn calculate_mass_imbalance_with_h(
        &self,
        fields: &SimulationFields<T>,
        h_operator: &Field2D<Vector2<T>>,
        i: usize,
        j: usize,
    ) -> T {
        // Standard mass imbalance
        let mass_imbalance = self.calculate_mass_imbalance(fields, i, j);

        // Add H(u) correction terms
        let two_t = scalar::one::<T>() + scalar::one::<T>();
        let h_correction_x =
            (h_operator.at(i + 1, j)[0] - h_operator.at(i - 1, j)[0]) / (two_t * self.dx);
        let h_correction_y =
            (h_operator.at(i, j + 1)[1] - h_operator.at(i, j - 1)[1]) / (two_t * self.dy);

        mass_imbalance + h_correction_x + h_correction_y
    }

    /// Calculate mass imbalance at a cell
    fn calculate_mass_imbalance(&self, fields: &SimulationFields<T>, i: usize, j: usize) -> T {
        let rho = fields.density.at(i, j);

        // Face mass fluxes
        let me = rho * fields.u.at(i + 1, j) * self.dy;
        let mw = rho * fields.u.at(i, j) * self.dy;
        let mn = rho * fields.v.at(i, j + 1) * self.dx;
        let ms = rho * fields.v.at(i, j) * self.dx;

        // Mass imbalance (should be zero for continuity)
        me - mw + mn - ms
    }

    /// Correct pressure field
    fn correct_pressure(&self, fields: &mut SimulationFields<T>, p_prime: &Field2D<T>) {
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let p_correction = self.pressure_relaxation * p_prime.at(i, j);
                let p_val = fields.p.at(i, j);
                if let Some(p) = fields.p.at_mut(i, j) {
                    *p = p_val + p_correction;
                }
            }
        }
    }

    /// Correct velocity field based on pressure correction
    fn correct_velocity(&self, fields: &mut SimulationFields<T>, p_prime: &Field2D<T>, dt: T) {
        // Correct u-velocity
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let dp_dx = (p_prime.at(i, j) - p_prime.at(i - 1, j)) / self.dx;
                let u_correction = -dt * dp_dx / fields.density.at(i, j);
                let u_val = fields.u.at(i, j);
                if let Some(u) = fields.u.at_mut(i, j) {
                    *u = u_val + u_correction;
                }
            }
        }

        // Correct v-velocity
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let dp_dy = (p_prime.at(i, j) - p_prime.at(i, j - 1)) / self.dy;
                let v_correction = -dt * dp_dy / fields.density.at(i, j);
                let v_val = fields.v.at(i, j);
                if let Some(v) = fields.v.at_mut(i, j) {
                    *v = v_val + v_correction;
                }
            }
        }
    }

    /// Update face fluxes using the full Rhie-Chow interpolation.
    ///
    /// # Theorem (Rhie & Chow, 1983, Eq. 13)
    /// The face velocity includes both cell-centred and face pressure gradients:
    ///
    ///   u_f = ū_f + d_f · [(∇p)_cells − (∇p)_face]
    ///
    /// where:
    /// - `ū_f` = arithmetic mean of adjacent cell velocities
    /// - `d_f` = Δt / ρ_f, the transient pressure-momentum coupling coefficient
    /// - `(∇p)_cells` = average of the two adjacent cell-centred pressure gradients
    /// - `(∇p)_face` = direct face pressure gradient
    ///
    /// Omitting `(∇p)_cells` reduces this to plain pressure interpolation which does
    /// NOT suppress checkerboard oscillations. Both terms are mathematically mandatory.
    fn update_face_fluxes(&self, fields: &mut SimulationFields<T>, dt: T) {
        let tiny = scalar::from_f64::<T>(1e-30);
        for i in 1..self.nx - 1 {
            for j in 1..self.ny - 1 {
                let two_t = scalar::one::<T>() + scalar::one::<T>();

                // Transient Rhie-Chow coefficient: for the current uniform-cell
                // transient discretization, d_f reduces to Δt/ρ_f.
                let rho_face_x = (fields.density.at(i, j) + fields.density.at(i + 1, j)) / two_t;
                let rho_face_y = (fields.density.at(i, j) + fields.density.at(i, j + 1)) / two_t;
                let d_u = dt / <T as NumericElement>::max_scalar(rho_face_x, tiny);
                let d_v = dt / <T as NumericElement>::max_scalar(rho_face_y, tiny);

                // ─── U-velocity on east face ──────────────────────────────────────
                let u_bar = (fields.u.at(i, j) + fields.u.at(i + 1, j)) / two_t;

                // Cell-centred ∇p_x at cell i and i+1
                let dp_dx_p = if i > 0 && i < self.nx - 1 {
                    (fields.p.at(i + 1, j) - fields.p.at(i - 1, j)) / (two_t * self.dx)
                } else {
                    (fields.p.at(i + 1, j) - fields.p.at(i, j)) / self.dx
                };
                let dp_dx_e = if i + 1 < self.nx - 1 {
                    (fields.p.at(i + 2, j) - fields.p.at(i, j)) / (two_t * self.dx)
                } else {
                    (fields.p.at(i + 1, j) - fields.p.at(i, j)) / self.dx
                };
                let dp_dx_cells = (dp_dx_p + dp_dx_e) / two_t;
                let dp_dx_face = (fields.p.at(i + 1, j) - fields.p.at(i, j)) / self.dx;

                // Full Rhie-Chow: ū_f + d_f · [(∇p)_cells − (∇p)_face]
                let u_corrected = u_bar + d_u * (dp_dx_cells - dp_dx_face);

                // ─── V-velocity on north face ─────────────────────────────────────
                let v_bar = (fields.v.at(i, j) + fields.v.at(i, j + 1)) / two_t;

                let dp_dy_p = if j > 0 && j < self.ny - 1 {
                    (fields.p.at(i, j + 1) - fields.p.at(i, j - 1)) / (two_t * self.dy)
                } else {
                    (fields.p.at(i, j + 1) - fields.p.at(i, j)) / self.dy
                };
                let dp_dy_n = if j + 1 < self.ny - 1 {
                    (fields.p.at(i, j + 2) - fields.p.at(i, j)) / (two_t * self.dy)
                } else {
                    (fields.p.at(i, j + 1) - fields.p.at(i, j)) / self.dy
                };
                let dp_dy_cells = (dp_dy_p + dp_dy_n) / two_t;
                let dp_dy_face = (fields.p.at(i, j + 1) - fields.p.at(i, j)) / self.dy;

                let v_corrected = v_bar + d_v * (dp_dy_cells - dp_dy_face);

                if let Some(u) = fields.u.at_mut(i, j) {
                    *u = u_corrected;
                }
                if let Some(v) = fields.v.at_mut(i, j) {
                    *v = v_corrected;
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::SimulationFields;
    use crate::grid::StructuredGrid2D;

    fn make_grid(n: usize) -> StructuredGrid2D<f64> {
        StructuredGrid2D::new(n, n, 0.0, 1.0, 0.0, 1.0).unwrap()
    }

    #[test]
    fn corrector_creation_with_valid_config() {
        let grid = make_grid(8);
        let corrector = PressureCorrector::new(&grid, 2, 0.3);

        assert_eq!(corrector.nx, 8);
        assert_eq!(corrector.ny, 8);
        assert_eq!(corrector.num_correctors, 2);
        assert!((corrector.pressure_relaxation - 0.3).abs() < 1e-15);
    }

    #[test]
    fn single_correction_step_does_not_panic() {
        let grid = make_grid(8);
        let corrector = PressureCorrector::new(&grid, 1, 0.5);
        let mut fields: SimulationFields<f64> = SimulationFields::new(8, 8);

        // Set mild initial velocity for a non-trivial correction
        for i in 0..8 {
            for j in 0..8 {
                if let Some(u) = fields.u.at_mut(i, j) {
                    *u = 0.01 * (i as f64);
                }
                if let Some(v) = fields.v.at_mut(i, j) {
                    *v = -0.01 * (j as f64);
                }
            }
        }

        let dt = 0.001;
        let result = corrector.correct(&mut fields, dt);
        assert!(result.is_ok());

        // Verify all fields remain finite
        for i in 0..8 {
            for j in 0..8 {
                assert!(fields.u.at(i, j).is_finite(), "u[{i}][{j}] is not finite");
                assert!(fields.v.at(i, j).is_finite(), "v[{i}][{j}] is not finite");
                assert!(fields.p.at(i, j).is_finite(), "p[{i}][{j}] is not finite");
            }
        }
    }

    #[test]
    fn two_corrector_steps_remain_finite_after_face_update() {
        let grid = make_grid(8);
        let corrector = PressureCorrector::new(&grid, 2, 0.5);
        let mut fields: SimulationFields<f64> = SimulationFields::new(8, 8);

        for i in 0..8 {
            for j in 0..8 {
                if let Some(u) = fields.u.at_mut(i, j) {
                    *u = 0.01 * (i as f64);
                }
                if let Some(v) = fields.v.at_mut(i, j) {
                    *v = -0.01 * (j as f64);
                }
                if let Some(rho) = fields.density.at_mut(i, j) {
                    *rho = 1.0 + 0.001 * ((i + j) as f64);
                }
            }
        }

        let result = corrector.correct(&mut fields, 0.001);
        assert!(result.is_ok());

        for i in 0..8 {
            for j in 0..8 {
                assert!(fields.u.at(i, j).is_finite(), "u[{i}][{j}] is not finite");
                assert!(fields.v.at(i, j).is_finite(), "v[{i}][{j}] is not finite");
                assert!(fields.p.at(i, j).is_finite(), "p[{i}][{j}] is not finite");
            }
        }
    }
}
