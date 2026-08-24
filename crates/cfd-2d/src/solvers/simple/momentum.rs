use super::algorithm::SimpleAlgorithm;
use super::STAGNANT_CELL_AP_THRESHOLD;
use crate::fields::SimulationFields;
use crate::grid::StructuredGrid2D;
use crate::physics::momentum::{MomentumComponent, MomentumSolver};
use crate::scalar;
use cfd_core::error::Result;
use cfd_core::CfdScalar;
use eunomia::{FloatElement, NumericElement, RealField as EunomiaRealField};

impl<T: CfdScalar + EunomiaRealField + Copy + std::fmt::Debug + FloatElement> SimpleAlgorithm<T> {
    pub(crate) fn predict_momentum(
        &mut self,
        momentum_solver: &mut MomentumSolver<T>,
        fields: &mut SimulationFields<T>,
        dt: T,
    ) -> Result<()> {
        momentum_solver.set_velocity_relaxation(self.velocity_relaxation);
        momentum_solver.solve(MomentumComponent::U, fields, dt)?;
        momentum_solver.solve(MomentumComponent::V, fields, dt)?;
        Ok(())
    }

    pub(crate) fn compute_d_coefficients(
        &mut self,
        momentum_solver: &MomentumSolver<T>,
        grid: &StructuredGrid2D<T>,
    ) {
        let (ap_u, _, ap_v, _) = momentum_solver.get_ap_coefficients();
        let dx = grid.dx;
        let dy = grid.dy;
        let min_ap = <T as FloatElement>::from_f64(STAGNANT_CELL_AP_THRESHOLD);

        let d_u = self.d_u.as_mut().expect("expected value");
        let d_v = self.d_v.as_mut().expect("expected value");

        for j in 0..grid.ny {
            for i in 0..grid.nx {
                // d = V / |a_P| (Patankar 1980): the velocity response to a
                // unit pressure gradient, δu = −(V/a_P)∇p′. Both the p′
                // Poisson coefficients and the Rhie–Chow interpolation must
                // use this same volume-based coefficient to remain mutually
                // consistent. Prior to 2026-08-20 d_u = dy/|a_P| omitted the
                // dx factor, leaving the p′ equation inconsistent with the
                // applied nodal correction by 1/dx on non-square grids.
                let cell_volume = dx * dy;

                let ap_u_val = NumericElement::abs(ap_u.at(i, j));
                d_u.set(
                    i,
                    j,
                    if ap_u_val > min_ap {
                        cell_volume / ap_u_val
                    } else {
                        scalar::zero::<T>()
                    },
                );

                let ap_v_val = NumericElement::abs(ap_v.at(i, j));
                d_v.set(
                    i,
                    j,
                    if ap_v_val > min_ap {
                        cell_volume / ap_v_val
                    } else {
                        scalar::zero::<T>()
                    },
                );
            }
        }
    }
}
