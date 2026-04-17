//! Rhie-Chow velocity interpolation and boundary condition enforcement
//!
//! Provides consistent face-velocity interpolation using the Rhie-Chow formulation
//! and boundary condition application on the face-velocity field.
//!
//! # Theorem (Rhie-Chow Consistency — Rhie & Chow 1983)
//!
//! For collocated grid arrangements, the pressure correction equation must use face
//! velocities interpolated consistently with the momentum equation discretization to
//! prevent checkerboard pressure oscillations.
//!
//! The Rhie-Chow face velocity is:
//!
//! ```text
//! u_f = ū_f + d_f · [(∇p)_P - (∇p)_f] + transient correction
//! ```
//!
//! where `ū_f` is the linear interpolation of cell-centred velocities, `d_f = (V/a_P)_f`
//! is the face pressure-gradient coefficient, `(∇p)_P` is the cell-centred pressure
//! gradient, and `(∇p)_f` is the face-normal pressure gradient.
//!
//! **Proof sketch**: The momentum equation produces cell-centred velocities whose linear
//! interpolation to faces does not reproduce the discrete pressure gradient stencil
//! (3-point vs 2-point). The Rhie-Chow correction adds the difference between the
//! cell-centred and face-centred pressure gradients, scaled by `d_f`, which cancels the
//! even-odd decoupling mode and enforces a compact 3-point Laplacian for pressure.
//! This is equivalent to adding a fourth-order artificial dissipation that vanishes
//! at grid convergence.
//!
//! ## References
//! - Rhie, C.M. & Chow, W.L. (1983). *AIAA Journal*, 21(11), 1525–1532.

use super::solver::SimplecPimpleSolver;
use crate::fields::SimulationFields;
use crate::grid::array2d::Array2D;
use crate::pressure_velocity::RhieChowInterpolation;
use nalgebra::{RealField, Vector2};
use num_traits::{FromPrimitive, ToPrimitive};

impl<T: RealField + Copy + FromPrimitive + ToPrimitive + std::fmt::LowerExp>
    SimplecPimpleSolver<T>
{
    /// Interpolate consistent face velocities using Rhie-Chow interpolation
    pub(super) fn interpolate_consistent_velocity(
        &self,
        rhie_chow: &RhieChowInterpolation<T>,
        fields: &SimulationFields<T>,
        dt: Option<T>,
        velocity_field: &mut crate::fields::Field2D<nalgebra::Vector2<T>>,
        consistent_velocity: &mut Array2D<Vector2<T>>,
    ) {
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                if i < self.grid.nx - 1 {
                    let u_face = rhie_chow.face_velocity_x(
                        velocity_field,
                        &fields.p,
                        self.grid.dx,
                        self.grid.dy,
                        dt,
                        i,
                        j,
                    );
                    consistent_velocity[(i, j)].x = u_face;
                }

                if j < self.grid.ny - 1 {
                    let v_face = rhie_chow.face_velocity_y(
                        velocity_field,
                        &fields.p,
                        self.grid.dx,
                        self.grid.dy,
                        dt,
                        i,
                        j,
                    );
                    consistent_velocity[(i, j)].y = v_face;
                }
            }
        }

        self.apply_velocity_boundary_conditions(consistent_velocity, fields);
    }

    /// Apply boundary conditions to face velocities
    pub(super) fn apply_velocity_boundary_conditions(
        &self,
        face_velocity: &mut Array2D<Vector2<T>>,
        _fields: &SimulationFields<T>,
    ) {
        use cfd_core::physics::boundary::{BoundaryCondition, WallType};

        let bcs = self.momentum_solver.boundary_conditions();

        // North boundary (j = ny-2)
        if let Some(BoundaryCondition::Wall { wall_type }) = bcs.get("north") {
            for i in 0..self.grid.nx {
                match wall_type {
                    WallType::NoSlip => {
                        face_velocity[(i, self.grid.ny - 2)] = Vector2::zeros();
                    }
                    WallType::Moving { velocity } => {
                        face_velocity[(i, self.grid.ny - 2)] =
                            Vector2::new(velocity[0], velocity[1]);
                    }
                    _ => {}
                }
            }
        }

        // South boundary (j = 0)
        if let Some(BoundaryCondition::Wall { wall_type }) = bcs.get("south") {
            for i in 0..self.grid.nx {
                match wall_type {
                    WallType::NoSlip => {
                        face_velocity[(i, 0)] = Vector2::zeros();
                    }
                    WallType::Moving { velocity } => {
                        face_velocity[(i, 0)] = Vector2::new(velocity[0], velocity[1]);
                    }
                    _ => {}
                }
            }
        }

        // West boundary (i = 0)
        if let Some(BoundaryCondition::Wall { wall_type }) = bcs.get("west") {
            // Leave the corner nodes to the horizontal walls so the lid owns
            // the top corners and the lower wall owns the bottom corners.
            for j in 1..self.grid.ny {
                if j == self.grid.ny - 2 {
                    continue;
                }
                match wall_type {
                    WallType::NoSlip => {
                        face_velocity[(0, j)] = Vector2::zeros();
                    }
                    WallType::Moving { velocity } => {
                        face_velocity[(0, j)] = Vector2::new(velocity[0], velocity[1]);
                    }
                    _ => {}
                }
            }
        }

        // East boundary (i = nx-2)
        if let Some(BoundaryCondition::Wall { wall_type }) = bcs.get("east") {
            // Leave the corner nodes to the horizontal walls so the lid owns
            // the top corners and the lower wall owns the bottom corners.
            for j in 1..self.grid.ny {
                if j == self.grid.ny - 2 {
                    continue;
                }
                match wall_type {
                    WallType::NoSlip => {
                        face_velocity[(self.grid.nx - 2, j)] = Vector2::zeros();
                    }
                    WallType::Moving { velocity } => {
                        face_velocity[(self.grid.nx - 2, j)] =
                            Vector2::new(velocity[0], velocity[1]);
                    }
                    _ => {}
                }
            }
        }
    }

    /// Solve pressure correction equation with Rhie-Chow awareness
    pub(super) fn solve_pressure_correction(
        &self,
        fields: &mut SimulationFields<T>,
        dt: T,
        rho: T,
        rebuild_matrix: bool,
        output_correction: &mut Array2D<T>,
    ) -> cfd_core::error::Result<()> {
        if let Some(ref rhie_chow) = self.rhie_chow {
            let mut vfc = self._vel_field_cache.borrow_mut();
            if vfc.as_ref().is_none_or(|v| {
                let (nx, ny) = v.dimensions();
                nx != self.grid.nx || ny != self.grid.ny
            }) {
                *vfc = Some(crate::fields::Field2D::new(
                    self.grid.nx,
                    self.grid.ny,
                    Vector2::zeros(),
                ));
            }
            let velocity_field = vfc.as_mut().unwrap();

            let mut cvc = self._cons_vel_cache.borrow_mut();
            if cvc
                .as_ref()
                .is_none_or(|v| v.rows() != self.grid.nx || v.cols() != self.grid.ny)
            {
                *cvc = Some(Array2D::new(self.grid.nx, self.grid.ny, Vector2::zeros()));
            }
            let consistent_velocity = cvc.as_mut().unwrap();

            self.interpolate_consistent_velocity(
                rhie_chow,
                fields,
                Some(dt),
                velocity_field,
                consistent_velocity,
            );

            let mut ufc = self._u_face_cache.borrow_mut();
            if ufc
                .as_ref()
                .is_none_or(|v| v.rows() != self.grid.nx - 1 || v.cols() != self.grid.ny)
            {
                *ufc = Some(Array2D::new(self.grid.nx - 1, self.grid.ny, T::zero()));
            }
            let u_face = ufc.as_mut().unwrap();

            let mut vfc_out = self._v_face_cache.borrow_mut();
            if vfc_out
                .as_ref()
                .is_none_or(|v| v.rows() != self.grid.nx || v.cols() != self.grid.ny - 1)
            {
                *vfc_out = Some(Array2D::new(self.grid.nx, self.grid.ny - 1, T::zero()));
            }
            let v_face = vfc_out.as_mut().unwrap();

            let mut dxc = self._d_x_cache.borrow_mut();
            if dxc
                .as_ref()
                .is_none_or(|v| v.rows() != self.grid.nx - 1 || v.cols() != self.grid.ny)
            {
                *dxc = Some(Array2D::new(self.grid.nx - 1, self.grid.ny, T::zero()));
            }
            let d_x = dxc.as_mut().unwrap();

            let mut dyc = self._d_y_cache.borrow_mut();
            if dyc
                .as_ref()
                .is_none_or(|v| v.rows() != self.grid.nx || v.cols() != self.grid.ny - 1)
            {
                *dyc = Some(Array2D::new(self.grid.nx, self.grid.ny - 1, T::zero()));
            }
            let d_y = dyc.as_mut().unwrap();

            for i in 0..self.grid.nx - 1 {
                for j in 0..self.grid.ny {
                    u_face[(i, j)] = consistent_velocity[(i, j)].x;
                    d_x[(i, j)] = rhie_chow.d_face_x(i, j, self.grid.dx, self.grid.dy);
                }
            }
            for i in 0..self.grid.nx {
                for j in 0..self.grid.ny - 1 {
                    v_face[(i, j)] = consistent_velocity[(i, j)].y;
                    d_y[(i, j)] = rhie_chow.d_face_y(i, j, self.grid.dx, self.grid.dy);
                }
            }

            let boundary_conditions = self.momentum_solver.boundary_conditions();
            self.pressure_solver.solve_pressure_correction_from_faces(
                u_face,
                v_face,
                d_x,
                d_y,
                rho,
                fields,
                boundary_conditions,
                rebuild_matrix,
                output_correction,
            )
        } else {
            self.pressure_solver.solve_pressure_correction(
                fields,
                dt,
                rho,
                rebuild_matrix,
                output_correction,
            )
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::fields::SimulationFields;
    use crate::grid::array2d::Array2D;
    use crate::grid::StructuredGrid2D;
    use cfd_core::physics::boundary::WallType;
    use nalgebra::Vector2;

    #[test]
    fn north_wall_owns_the_top_corners() {
        let grid = StructuredGrid2D::new(4, 4, 0.0_f64, 1.0_f64, 0.0_f64, 1.0_f64)
            .expect("grid creation failed");
        let config = crate::simplec_pimple::config::SimplecPimpleConfig::simplec();
        let mut solver =
            SimplecPimpleSolver::new(grid.clone(), config).expect("solver creation failed");
        solver.set_boundary(
            "north".to_string(),
            cfd_core::physics::boundary::BoundaryCondition::Wall {
                wall_type: WallType::Moving {
                    velocity: nalgebra::Vector3::new(1.0, 0.0, 0.0),
                },
            },
        );
        solver.set_boundary(
            "south".to_string(),
            cfd_core::physics::boundary::BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        );
        solver.set_boundary(
            "west".to_string(),
            cfd_core::physics::boundary::BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        );
        solver.set_boundary(
            "east".to_string(),
            cfd_core::physics::boundary::BoundaryCondition::Wall {
                wall_type: WallType::NoSlip,
            },
        );

        let mut face_velocity = Array2D::new(grid.nx, grid.ny, Vector2::new(-1.0, -1.0));
        let fields = SimulationFields::new(grid.nx, grid.ny);

        solver.apply_velocity_boundary_conditions(&mut face_velocity, &fields);

        assert_eq!(face_velocity[(0, grid.ny - 2)], Vector2::new(1.0, 0.0));
        assert_eq!(
            face_velocity[(grid.nx - 2, grid.ny - 2)],
            Vector2::new(1.0, 0.0)
        );
        assert_eq!(face_velocity[(0, 1)], Vector2::zeros());
        assert_eq!(face_velocity[(grid.nx - 2, 1)], Vector2::zeros());
        assert_eq!(face_velocity[(1, 0)], Vector2::zeros());
    }
}
