use super::algorithm::SimpleAlgorithm;
use super::STAGNANT_CELL_AP_THRESHOLD;
use crate::fields::SimulationFields;
use crate::grid::StructuredGrid2D;
use crate::physics::momentum::{MomentumComponent, MomentumSolver};
use cfd_core::error::Result;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::Debug> SimpleAlgorithm<T> {
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
        let min_ap = T::from_f64(STAGNANT_CELL_AP_THRESHOLD).unwrap();

        let d_u = self.d_u.as_mut().unwrap();
        let d_v = self.d_v.as_mut().unwrap();

        for j in 0..grid.ny {
            for i in 0..grid.nx {
                let ap_u_val = ap_u.at(i, j).abs();
                d_u.set(
                    i,
                    j,
                    if ap_u_val > min_ap {
                        dy / ap_u_val
                    } else {
                        T::zero()
                    },
                );

                let ap_v_val = ap_v.at(i, j).abs();
                d_v.set(
                    i,
                    j,
                    if ap_v_val > min_ap {
                        dx / ap_v_val
                    } else {
                        T::zero()
                    },
                );
            }
        }
    }
}
