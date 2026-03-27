use super::algorithm::SimpleAlgorithm;
use crate::fields::SimulationFields;
use crate::grid::StructuredGrid2D;
use cfd_core::error::Result;
use cfd_math::linear_solver::preconditioners::IdentityPreconditioner;
use cfd_math::linear_solver::{BiCGSTAB, IterativeLinearSolver, IterativeSolverConfig};
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::RealField;
use num_traits::{FromPrimitive, ToPrimitive};

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::Debug> SimpleAlgorithm<T> {
    pub(crate) fn assemble_pressure_correction(
        &mut self,
        fields: &SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> Result<T> {
        let nx = grid.nx;
        let ny = grid.ny;
        let dx = grid.dx;
        let dy = grid.dy;
        let half = T::from_f64(0.5).unwrap();
        let two = T::from_f64(2.0).unwrap();

        let matrix_builder = self.matrix_builder.as_mut().unwrap();
        matrix_builder.clear();
        let rhs = self.rhs.as_mut().unwrap();
        rhs.fill(T::zero());
        let p_prime = self.p_prime.as_mut().unwrap();
        p_prime.fill(T::zero());

        let d_u = self.d_u.as_ref().unwrap();
        let d_v = self.d_v.as_ref().unwrap();

        let mut max_residual = T::zero();

        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;

                let u_p = fields.u.at(i, j);
                let u_e = fields.u.at(i + 1, j);
                let u_w = fields.u.at(i - 1, j);
                let v_p = fields.v.at(i, j);
                let v_n = fields.v.at(i, j + 1);
                let v_s = fields.v.at(i, j - 1);

                let p_p = fields.p.at(i, j);
                let p_e = fields.p.at(i + 1, j);
                let p_w = fields.p.at(i - 1, j);
                let p_n = fields.p.at(i, j + 1);
                let p_s = fields.p.at(i, j - 1);

                let p_ee = if i + 2 < nx { fields.p.at(i + 2, j) } else { two * p_e - p_p };
                let p_ww = if i >= 2 { fields.p.at(i - 2, j) } else { two * p_w - p_p };
                let p_nn = if j + 2 < ny { fields.p.at(i, j + 2) } else { two * p_n - p_p };
                let p_ss = if j >= 2 { fields.p.at(i, j - 2) } else { two * p_s - p_p };

                let d_face_e = (d_u.at(i, j) + d_u.at(i + 1, j)) * half;
                let d_face_w = (d_u.at(i - 1, j) + d_u.at(i, j)) * half;
                let d_face_n = (d_v.at(i, j) + d_v.at(i, j + 1)) * half;
                let d_face_s = (d_v.at(i, j - 1) + d_v.at(i, j)) * half;

                let grad_p_p_x = (p_e - p_w) / (two * dx);
                let avg_grad_pe = (grad_p_p_x + (p_ee - p_p) / (two * dx)) * half;
                let u_face_e = (u_p + u_e) * half + d_face_e * (avg_grad_pe - (p_e - p_p) / dx);

                let avg_grad_pw = ((p_p - p_ww) / (two * dx) + grad_p_p_x) * half;
                let u_face_w = (u_w + u_p) * half + d_face_w * (avg_grad_pw - (p_p - p_w) / dx);

                let grad_p_p_y = (p_n - p_s) / (two * dy);
                let avg_grad_pn = (grad_p_p_y + (p_nn - p_p) / (two * dy)) * half;
                let v_face_n = (v_p + v_n) * half + d_face_n * (avg_grad_pn - (p_n - p_p) / dy);

                let avg_grad_ps = ((p_p - p_ss) / (two * dy) + grad_p_p_y) * half;
                let v_face_s = (v_s + v_p) * half + d_face_s * (avg_grad_ps - (p_p - p_s) / dy);

                let rho = fields.density.at(i, j);
                let flux_e = rho * u_face_e * dy;
                let flux_w = rho * u_face_w * dy;
                let flux_n = rho * v_face_n * dx;
                let flux_s = rho * v_face_s * dx;
                let mass_imbalance = flux_e - flux_w + flux_n - flux_s;

                let a_e = rho * d_face_e * dy / dx;
                let a_w = rho * d_face_w * dy / dx;
                let a_n = rho * d_face_n * dx / dy;
                let a_s = rho * d_face_s * dx / dy;
                let a_p = -(a_e + a_w + a_n + a_s);

                matrix_builder.add_entry(idx, idx, a_p)?;
                matrix_builder.add_entry(idx, idx + 1, a_e)?;
                matrix_builder.add_entry(idx, idx - 1, a_w)?;
                matrix_builder.add_entry(idx, idx + nx, a_n)?;
                matrix_builder.add_entry(idx, idx - nx, a_s)?;

                rhs[idx] = mass_imbalance;

                if mass_imbalance.abs() > max_residual {
                    max_residual = mass_imbalance.abs();
                }
            }
        }

        for j in 0..ny {
            for i in 0..nx {
                if i > 0 && i < nx - 1 && j > 0 && j < ny - 1 {
                    continue;
                }
                let idx = j * nx + i;
                if i == 0 || i == nx - 1 {
                    matrix_builder.add_entry(idx, idx, T::one())?;
                    rhs[idx] = T::zero();
                } else {
                    let neighbor = if j == 0 { idx + nx } else { idx - nx };
                    matrix_builder.add_entry(idx, idx, T::one())?;
                    matrix_builder.add_entry(idx, neighbor, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
        }
        Ok(max_residual)
    }

    pub(crate) fn solve_pressure_correction(&mut self) -> Result<()> {
        let builder_owned = self.matrix_builder.take().unwrap();
        let n_rows = builder_owned.num_rows();
        let matrix = builder_owned.build()?;
        self.matrix_builder = Some(SparseMatrixBuilder::new(n_rows, n_rows));
        let solver_config = IterativeSolverConfig {
            tolerance: self.tolerance,
            max_iterations: 2000,
            ..Default::default()
        };
        let linear_solver = BiCGSTAB::new(solver_config);
        linear_solver.solve(
            &matrix,
            self.rhs.as_ref().unwrap(),
            self.p_prime.as_mut().unwrap(),
            None::<&IdentityPreconditioner>,
        )?;

        let pp = self.p_prime.as_mut().unwrap();
        if pp.iter().any(|v| !v.is_finite()) {
            pp.fill(T::zero());
        }
        Ok(())
    }

    pub(crate) fn apply_corrections(
        &mut self,
        fields: &mut SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) {
        let nx = grid.nx;
        let ny = grid.ny;
        let dx = grid.dx;
        let dy = grid.dy;
        let two = T::from_f64(2.0).unwrap();

        let p_prime = self.p_prime.as_ref().unwrap();
        let d_u = self.d_u.as_ref().unwrap();
        let d_v = self.d_v.as_ref().unwrap();

        for j in 1..ny - 1 {
            for i in 1..nx - 1 {
                let idx = j * nx + i;
                let pp = p_prime[idx];
                
                if let Some(p) = fields.p.at_mut(i, j) {
                    *p += self.pressure_relaxation * pp;
                }

                let dp_dx = (p_prime[idx + 1] - p_prime[idx - 1]) / (two * dx);
                let dp_dy = (p_prime[idx + nx] - p_prime[idx - nx]) / (two * dy);

                if let Some(u) = fields.u.at_mut(i, j) {
                    *u -= d_u.at(i, j) * dp_dx * dx;
                }
                if let Some(v) = fields.v.at_mut(i, j) {
                    *v -= d_v.at(i, j) * dp_dy * dy;
                }
            }
        }
    }
}
