//! Mass-flux and face-interpolation solvers.

use crate::scalar;
use crate::scalar::Cfd2dScalar;
use crate::solvers::ns_fvm::solver::NavierStokesSolver2D;
use eunomia::{FloatElement, NumericElement};

impl<T: Cfd2dScalar + eunomia::RealField + Copy + FloatElement> NavierStokesSolver2D<T> {
    /// Applies the mass flux correction from the solved pressure field.
    ///
    /// Updates face velocities (`u` and `v`) based on the pressure correction gradients
    /// to satisfy continuity, and updates the cell-centred pressure `p`.
    pub fn apply_mass_flux_correction(&mut self) {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let zero: T = scalar::zero();

        let mut q_in = zero;
        for j in 0..ny {
            if self.field.mask[(0, j)] {
                let dy_j = self.grid.dy_at(j);
                q_in += self.field.u[(0, j)] * dy_j;
            }
        }

        let mut q_out = zero;
        for j in 0..ny {
            if self.field.mask[(nx - 1, j)] {
                let dy_j = self.grid.dy_at(j);
                q_out += self.field.u[(nx, j)] * dy_j;
            }
        }

        let tiny: T = scalar::from_f64(1e-30);
        if <T as NumericElement>::abs(q_out) > tiny && <T as NumericElement>::abs(q_in) > tiny {
            let scale = q_in / q_out;
            for j in 0..ny {
                if self.field.mask[(nx - 1, j)] {
                    self.field.u[(nx, j)] *= scale;
                }
            }
        }
    }
}
