use crate::scalar;
use crate::scalar::Cfd2dScalar;
use crate::solvers::ns_fvm::solver::NavierStokesSolver2D;
use eunomia::{FloatElement, NumericElement};

impl<T: Cfd2dScalar + Copy + FloatElement> NavierStokesSolver2D<T> {
    /// Compute separate L2-norm residuals for convergence assessment.
    ///
    /// Returns (res_continuity_rms, res_continuity_l1, res_max_pointwise):
    /// - `res_continuity`: RMS mass imbalance across all cells
    /// - `res_continuity_l1`: mean absolute mass imbalance across all cells
    /// - `res_max_pointwise`: L-infinity norm (maximum pointwise continuity error)
    ///
    /// The separate residuals allow distinguishing between:
    /// - Oscillating pressure (high res_max but moderate res_continuity)
    /// - Globally poor convergence (both high)
    /// - Localized divergence (high res_max, low res_continuity)
    pub fn compute_residuals(&self) -> (T, T, T) {
        self.compute_residuals_inner()
    }

    pub(super) fn compute_residuals_inner(&self) -> (T, T, T) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let rho = self.density;

        let mut cont_sum: T = scalar::zero();
        let mut cont_l1_sum: T = scalar::zero();
        let mut max_imb: T = scalar::zero();
        let mut n_fluid = 0_usize;

        for i in 0..nx {
            for j in 0..ny {
                if !self.field.mask[(i, j)] {
                    continue;
                }
                n_fluid += 1;
                let dy_j = self.grid.dy_at(j);

                // Continuity residual: div(rho * u) per cell.
                let imb = rho
                    * ((self.field.u[(i + 1, j)] - self.field.u[(i, j)]) * dy_j
                        + (self.field.v[(i, j + 1)] - self.field.v[(i, j)]) * dx);
                cont_sum += imb * imb;

                // L-infinity: track max pointwise imbalance.
                let abs_imb = <T as NumericElement>::abs(imb);
                if abs_imb > max_imb {
                    max_imb = abs_imb;
                }

                // Continuity L1 norm over fluid cells.
                cont_l1_sum += <T as NumericElement>::abs(imb);
            }
        }

        let n: T = scalar::from_usize(n_fluid.max(1));
        let res_cont = <T as NumericElement>::sqrt(cont_sum / n);
        let res_cont_l1 = cont_l1_sum / n;
        (res_cont, res_cont_l1, max_imb)
    }

    /// Check for NaN/Inf in the velocity field (divergence guard).
    pub fn check_divergence(&self) -> bool {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        for i in 0..=nx {
            for j in 0..ny {
                let u = self.field.u[(i, j)];
                if !<T as NumericElement>::is_finite(u) {
                    return true;
                }
            }
        }
        for i in 0..nx {
            for j in 0..=ny {
                let v = self.field.v[(i, j)];
                if !<T as NumericElement>::is_finite(v) {
                    return true;
                }
            }
        }
        false
    }
}
