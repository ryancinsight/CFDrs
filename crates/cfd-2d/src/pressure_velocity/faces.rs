//! Face-based pressure correction equation assembly and solution.
//!
//! Separated from correction.rs to adhere to file length limits.

use super::pressure::PressureCorrectionSolver;
use crate::grid::array2d::Array2D;
use crate::physics::momentum::validate_boundary_consistency;
use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_math::sparse::SparseMatrixBuilder;
use eunomia::{FloatElement, NumericElement};
use leto::Array1;
use leto_ops::norm_l2;
use std::fmt::Debug;

impl<T: Cfd2dScalar + Copy + Debug + FloatElement> PressureCorrectionSolver<T> {
    /// Solve pressure correction equation using face velocities (Rhie-Chow)
    pub fn solve_pressure_correction_from_faces(
        &self,
        u_face: &Array2D<T>,
        v_face: &Array2D<T>,
        d_x: &Array2D<T>,
        d_y: &Array2D<T>,
        _rho: T,
        fields: &crate::fields::SimulationFields<T>,
        boundary_conditions: &std::collections::HashMap<
            String,
            cfd_core::physics::boundary::BoundaryCondition<T>,
        >,
        _rebuild_matrix: bool,
        output_correction: &mut Array2D<T>,
    ) -> cfd_core::error::Result<()> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        let dx = self.grid.dx;
        let dy = self.grid.dy;

        validate_boundary_consistency(boundary_conditions, &self.grid)
            .map_err(|error| cfd_core::error::Error::InvalidConfiguration(error.to_string()))?;

        let n = (nx - 2) * (ny - 2);
        if n <= 1 {
            for i in 0..nx {
                for j in 0..ny {
                    output_correction[(i, j)] = scalar::zero::<T>();
                }
            }
            return Ok(());
        }

        let is_dirichlet = |side: &str| -> bool {
            if let Some(bc) = boundary_conditions.get(side) {
                matches!(
                    bc,
                    cfd_core::physics::boundary::BoundaryCondition::PressureOutlet { .. }
                        | cfd_core::physics::boundary::BoundaryCondition::PressureInlet { .. }
                        | cfd_core::physics::boundary::BoundaryCondition::CharacteristicOutlet { .. }
                )
            } else {
                false
            }
        };

        let has_dirichlet = is_dirichlet("west")
            || is_dirichlet("east")
            || is_dirichlet("south")
            || is_dirichlet("north");

        let (system_size, reference_idx) = if has_dirichlet {
            (n, None)
        } else {
            let ref_idx = (1..nx - 1)
                .flat_map(|i| (1..ny - 1).map(move |j| (i, j)))
                .enumerate()
                .find(|(_, (i, j))| fields.mask.at(*i, *j))
                .map_or(0usize, |(idx, _)| idx);
            (n - 1, Some(ref_idx))
        };

        let map_index = |idx: usize| -> Option<usize> {
            if let Some(ref_idx) = reference_idx {
                match idx.cmp(&ref_idx) {
                    std::cmp::Ordering::Equal => None,
                    std::cmp::Ordering::Less => Some(idx),
                    std::cmp::Ordering::Greater => Some(idx - 1),
                }
            } else {
                Some(idx)
            }
        };

        let mut builder = self.take_matrix_builder(system_size, system_size);
        let mut rhs = self
            ._rhs_cache
            .borrow_mut()
            .take()
            .filter(|vector| vector.shape()[0] == system_size)
            .unwrap_or_else(|| Array1::from_elem([system_size], scalar::zero::<T>()));
        rhs.fill(scalar::zero::<T>());

        let dx2_inv = scalar::one::<T>() / (dx * dx);
        let dy2_inv = scalar::one::<T>() / (dy * dy);

        let mut max_residual = scalar::zero::<T>();

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                if Some(idx) == reference_idx {
                    continue;
                }
                let row_idx = map_index(idx).expect("row index must exist");

                if !fields.mask.at(i, j) {
                    builder.add_entry(row_idx, row_idx, scalar::one::<T>())?;
                    rhs[row_idx] = scalar::zero::<T>();
                    continue;
                }

                let mut ap = scalar::zero::<T>();

                // West neighbour
                let aw = d_x[(i - 1, j)] * dx2_inv;
                if i > 1 && fields.mask.at(i - 1, j) {
                    ap += aw;
                    if let Some(ci) = map_index(idx - (ny - 2)) {
                        builder.add_entry(row_idx, ci, -aw)?;
                    }
                } else if i == 1 && is_dirichlet("west") {
                    ap += aw;
                }

                // East neighbour
                let ae = d_x[(i, j)] * dx2_inv;
                if i < nx - 2 && fields.mask.at(i + 1, j) {
                    ap += ae;
                    if let Some(ci) = map_index(idx + (ny - 2)) {
                        builder.add_entry(row_idx, ci, -ae)?;
                    }
                } else if i == nx - 2 && is_dirichlet("east") {
                    ap += ae;
                }

                // South neighbour
                let as_ = d_y[(i, j - 1)] * dy2_inv;
                if j > 1 && fields.mask.at(i, j - 1) {
                    ap += as_;
                    if let Some(ci) = map_index(idx - 1) {
                        builder.add_entry(row_idx, ci, -as_)?;
                    }
                } else if j == 1 && is_dirichlet("south") {
                    ap += as_;
                }

                // North neighbour
                let an = d_y[(i, j)] * dy2_inv;
                if j < ny - 2 && fields.mask.at(i, j + 1) {
                    ap += an;
                    if let Some(ci) = map_index(idx + 1) {
                        builder.add_entry(row_idx, ci, -an)?;
                    }
                } else if j == ny - 2 && is_dirichlet("north") {
                    ap += an;
                }

                builder.add_entry(row_idx, row_idx, ap)?;

                let div_u = (u_face[(i, j)] - u_face[(i - 1, j)]) / dx
                    + (v_face[(i, j)] - v_face[(i, j - 1)]) / dy;
                rhs[row_idx] = -div_u;

                let abs_rhs = NumericElement::abs(rhs[row_idx]);
                if abs_rhs > max_residual {
                    max_residual = abs_rhs;
                }
            }
        }

        tracing::debug!(
            "Pressure Solve (faces): n={n}, rhs_norm={:?}, max_residual={max_residual:?}",
            norm_l2(&rhs.view())
                .expect("invariant: face-pressure RHS Leto vector has a valid layout")
        );

        let matrix = builder.build()?;

        self.reset_matrix_builder_cache(system_size, system_size);
        let mut p_correction_vec = self
            ._solution_cache
            .borrow_mut()
            .take()
            .filter(|vector| vector.shape()[0] == matrix.nrows())
            .unwrap_or_else(|| Array1::from_elem([matrix.nrows()], scalar::zero::<T>()));
        p_correction_vec.fill(scalar::zero::<T>());
        self.dispatch_solve(&matrix, &rhs, &mut p_correction_vec)?;

        self.scatter_correction(
            &p_correction_vec,
            reference_idx,
            &map_index,
            Some(boundary_conditions),
            output_correction,
        )?;

        *self._rhs_cache.borrow_mut() = Some(rhs);
        *self._solution_cache.borrow_mut() = Some(p_correction_vec);

        Ok(())
    }

    /// Scatter solution vector back to 2D grid with appropriate boundary conditions
    pub(super) fn scatter_correction(
        &self,
        solution: &Array1<T>,
        reference_idx: Option<usize>,
        map_index: &dyn Fn(usize) -> Option<usize>,
        boundary_conditions: Option<
            &std::collections::HashMap<String, cfd_core::physics::boundary::BoundaryCondition<T>>,
        >,
        output_correction: &mut Array2D<T>,
    ) -> cfd_core::error::Result<()> {
        let nx = self.grid.nx;
        let ny = self.grid.ny;
        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let idx = (i - 1) * (ny - 2) + (j - 1);
                let value = if Some(idx) == reference_idx {
                    scalar::zero::<T>()
                } else if let Some(col_idx) = map_index(idx) {
                    solution[col_idx]
                } else {
                    scalar::zero::<T>()
                };
                output_correction[(i, j)] = value;
            }
        }

        let is_dirichlet = |side: &str| -> bool {
            if let Some(bcs) = boundary_conditions {
                if let Some(bc) = bcs.get(side) {
                    return matches!(
                        bc,
                        cfd_core::physics::boundary::BoundaryCondition::PressureOutlet { .. }
                            | cfd_core::physics::boundary::BoundaryCondition::PressureInlet { .. }
                            | cfd_core::physics::boundary::BoundaryCondition::CharacteristicOutlet { .. }
                    );
                }
            }
            false
        };

        // South boundary
        if is_dirichlet("south") {
            for i in 0..nx {
                output_correction[(i, 0)] = scalar::zero::<T>();
            }
        } else {
            for i in 0..nx {
                output_correction[(i, 0)] = output_correction[(i, 1)];
            }
        }

        // North boundary
        if is_dirichlet("north") {
            for i in 0..nx {
                output_correction[(i, ny - 1)] = scalar::zero::<T>();
            }
        } else {
            for i in 0..nx {
                output_correction[(i, ny - 1)] = output_correction[(i, ny - 2)];
            }
        }

        // West boundary
        if is_dirichlet("west") {
            for j in 0..ny {
                output_correction[(0, j)] = scalar::zero::<T>();
            }
        } else {
            for j in 0..ny {
                output_correction[(0, j)] = output_correction[(1, j)];
            }
        }

        // East boundary
        if is_dirichlet("east") {
            for j in 0..ny {
                output_correction[(nx - 1, j)] = scalar::zero::<T>();
            }
        } else {
            for j in 0..ny {
                output_correction[(nx - 1, j)] = output_correction[(nx - 2, j)];
            }
        }

        Ok(())
    }

    pub(super) fn take_matrix_builder(&self, rows: usize, cols: usize) -> SparseMatrixBuilder<T> {
        self._matrix_builder_cache
            .borrow_mut()
            .take()
            .filter(|builder| builder.num_rows() == rows && builder.num_cols() == cols)
            .unwrap_or_else(|| SparseMatrixBuilder::new(rows, cols))
    }

    pub(super) fn reset_matrix_builder_cache(&self, rows: usize, cols: usize) {
        *self._matrix_builder_cache.borrow_mut() = Some(SparseMatrixBuilder::new(rows, cols));
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::SimulationFields;
    use crate::grid::array2d::Array2D;
    use crate::grid::StructuredGrid2D;
    use crate::pressure_velocity::config::PressureLinearSolver;
    use cfd_core::physics::boundary::BoundaryCondition;
    use std::collections::HashMap;

    #[test]
    fn face_based_pressure_correction_with_walls_is_well_posed() {
        let grid = StructuredGrid2D::new(6, 6, 0.0, 1.0, 0.0, 1.0).unwrap();
        let solver = PressureCorrectionSolver::new(
            grid.clone(),
            PressureLinearSolver::GMRES { restart_dim: 10 },
        )
        .unwrap();
        let fields: SimulationFields<f64> = SimulationFields::new(6, 6);

        let boundary_conditions = HashMap::from([
            ("west".to_string(), BoundaryCondition::wall_no_slip()),
            ("east".to_string(), BoundaryCondition::wall_no_slip()),
            ("north".to_string(), BoundaryCondition::wall_no_slip()),
            ("south".to_string(), BoundaryCondition::wall_no_slip()),
        ]);

        let u_face = Array2D::new(grid.nx - 1, grid.ny, 0.0);
        let v_face = Array2D::new(grid.nx, grid.ny - 1, 0.0);
        let d_x = Array2D::new(grid.nx - 1, grid.ny, 1.0);
        let d_y = Array2D::new(grid.nx, grid.ny - 1, 1.0);
        let mut p_corr = Array2D::new(grid.nx, grid.ny, 1.0);

        solver
            .solve_pressure_correction_from_faces(
                &u_face,
                &v_face,
                &d_x,
                &d_y,
                1.0,
                &fields,
                &boundary_conditions,
                true,
                &mut p_corr,
            )
            .unwrap();

        for i in 0..grid.nx {
            for j in 0..grid.ny {
                let value = p_corr[(i, j)];
                assert!(value.is_finite(), "pressure correction must be finite");
                assert!(value.abs() < 1e-10, "expected zero correction, got {value}");
            }
        }
    }
}
