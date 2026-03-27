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
    ) -> Array2D<Vector2<T>> {
        let mut velocity_field =
            crate::fields::Field2D::new(self.grid.nx, self.grid.ny, Vector2::zeros());
        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                velocity_field.set(i, j, Vector2::new(fields.u.at(i, j), fields.v.at(i, j)));
            }
        }

        let mut consistent_velocity = Array2D::new(self.grid.nx, self.grid.ny, Vector2::zeros());

        for i in 0..self.grid.nx {
            for j in 0..self.grid.ny {
                if i < self.grid.nx - 1 {
                    let u_face = rhie_chow.face_velocity_x(
                        &velocity_field,
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
                        &velocity_field,
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

        self.apply_velocity_boundary_conditions(&mut consistent_velocity, fields);
        consistent_velocity
    }

    /// Apply boundary conditions to face velocities
    pub(super) fn apply_velocity_boundary_conditions(
        &self,
        face_velocity: &mut Array2D<Vector2<T>>,
        _fields: &SimulationFields<T>,
    ) {
        use cfd_core::physics::boundary::{BoundaryCondition, WallType};

        let bcs = self.momentum_solver.boundary_conditions();

        // North boundary (j = ny-1)
        if let Some(BoundaryCondition::Wall { wall_type }) = bcs.get("north") {
            for i in 0..self.grid.nx {
                match wall_type {
                    WallType::NoSlip => {
                        face_velocity[(i, self.grid.ny - 1)] = Vector2::zeros();
                    }
                    WallType::Moving { velocity } => {
                        face_velocity[(i, self.grid.ny - 1)] = Vector2::new(velocity[0], velocity[1]);
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
            for j in 0..self.grid.ny {
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

        // East boundary (i = nx-1)
        if let Some(BoundaryCondition::Wall { wall_type }) = bcs.get("east") {
            for j in 0..self.grid.ny {
                match wall_type {
                    WallType::NoSlip => {
                        face_velocity[(self.grid.nx - 1, j)] = Vector2::zeros();
                    }
                    WallType::Moving { velocity } => {
                        face_velocity[(self.grid.nx - 1, j)] = Vector2::new(velocity[0], velocity[1]);
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
    ) -> cfd_core::error::Result<Array2D<T>> {
        if let Some(ref rhie_chow) = self.rhie_chow {
            let consistent_velocity =
                self.interpolate_consistent_velocity(rhie_chow, fields, Some(dt));

            let mut u_face = Array2D::new(self.grid.nx - 1, self.grid.ny, T::zero());
            let mut v_face = Array2D::new(self.grid.nx, self.grid.ny - 1, T::zero());
            let mut d_x = Array2D::new(self.grid.nx - 1, self.grid.ny, T::zero());
            let mut d_y = Array2D::new(self.grid.nx, self.grid.ny - 1, T::zero());

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

            self.pressure_solver
                .solve_pressure_correction_from_faces(&u_face, &v_face, &d_x, &d_y, rho, fields)
        } else {
            self.pressure_solver
                .solve_pressure_correction(fields, dt, rho)
        }
    }
}
