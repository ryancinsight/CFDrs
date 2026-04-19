//! Boundary condition handling for LBM D2Q9.
//!
//! All methods accept the flat distribution buffer `f[j*nx*9 + i*9 + q]`,
//! flat density `rho[j*nx + i]`, and flat velocity `vel[(j*nx+i)*2 + d]`.
//!
//! # Theorem — Bounce-Back No-Slip (Ladd 1994)
//!
//! **Statement**: Full bounce-back at node $\mathbf{x}_w$ yields
//! $\mathbf{u}(\mathbf{x}_w) = 0$ to first-order accuracy (half-way scheme:
//! second order).
//!
//! **Proof**: Bounce-back sets $f_{\bar{q}}(\mathbf{x}_w) \leftarrow f_q(\mathbf{x}_w)$
//! for each antipodal pair $(q, \bar{q})$. The macroscopic velocity
//! $\mathbf{u} = \sum_q \mathbf{e}_q f_q / \rho$ sums terms $\mathbf{e}_q f_q
//! + \mathbf{e}_{\bar{q}} f_{\bar{q}} = \mathbf{e}_q f_q - \mathbf{e}_q f_q = 0$
//!   for each pair (because $\mathbf{e}_{\bar{q}} = -\mathbf{e}_q$), so $\mathbf{u} = 0$. □
//!
//! # Theorem — Zou-He Boundary Reconstruction (Zou & He 1997)
//!
//! **Statement**: On a rectangular D2Q9 boundary, the missing populations are
//! recovered by solving the discrete mass and momentum equations for the face
//! normal velocity and the two diagonal populations adjacent to that face.
//! For a west inlet, for example,
//!
//! $$
//! \rho = \frac{f_0 + f_2 + f_4 + 2(f_3 + f_6 + f_7)}{1 - u_x}
//! $$
//!
//! and the remaining unknown populations follow the standard Zou-He closure.
//!
//! **Proof**: The D2Q9 population set gives one mass equation and two momentum
//! equations. On a boundary face, three populations are unknown, so the system
//! is closed by the standard tangential non-equilibrium relation from Zou & He
//! (1997, §3). Solving that local linear system yields the face-specific update.
//! □

use crate::solvers::lbm::lattice::D2Q9;
use cfd_core::physics::boundary::BoundaryCondition;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use std::collections::HashMap;

/// Types of LBM boundary conditions
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryType {
    /// Bounce-back (no-slip wall)
    BounceBack,
    /// Velocity boundary (Dirichlet)
    Velocity,
    /// Pressure boundary (Dirichlet)
    Pressure,
    /// Open boundary (Neumann)
    Open,
    /// Periodic boundary
    Periodic,
}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
enum BoundaryFace {
    West,
    East,
    South,
    North,
}

/// Boundary handler for applying boundary conditions
pub struct BoundaryHandler<T: RealField + Copy> {
    /// Boundary type for each edge
    boundary_types: HashMap<String, BoundaryType>,
    /// Phantom data for type parameter
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive> Default for BoundaryHandler<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive> BoundaryHandler<T> {
    /// Create new boundary handler
    #[must_use]
    pub fn new() -> Self {
        Self {
            boundary_types: HashMap::new(),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Set boundary type for an edge
    pub fn set_boundary(&mut self, edge: String, boundary_type: BoundaryType) {
        self.boundary_types.insert(edge, boundary_type);
    }

    #[inline]
    fn boundary_face(
        i: usize,
        j: usize,
        nx: usize,
        ny: usize,
        velocity_hint: [T; 2],
    ) -> Option<BoundaryFace> {
        let on_west = i == 0;
        let on_east = i + 1 == nx;
        let on_south = j == 0;
        let on_north = j + 1 == ny;

        match (on_west || on_east, on_south || on_north) {
            (true, false) => Some(if on_west {
                BoundaryFace::West
            } else {
                BoundaryFace::East
            }),
            (false, true) => Some(if on_south {
                BoundaryFace::South
            } else {
                BoundaryFace::North
            }),
            (true, true) => {
                let velocity_x = velocity_hint[0].abs();
                let velocity_y = velocity_hint[1].abs();
                if velocity_x >= velocity_y {
                    Some(if on_west {
                        BoundaryFace::West
                    } else {
                        BoundaryFace::East
                    })
                } else {
                    Some(if on_south {
                        BoundaryFace::South
                    } else {
                        BoundaryFace::North
                    })
                }
            }
            (false, false) => None,
        }
    }

    #[inline]
    fn load_cell_populations(f: &[T], i: usize, j: usize, nx: usize) -> [T; 9] {
        use crate::solvers::lbm::streaming::f_idx;

        let mut values = [T::zero(); 9];
        for q in 0..9 {
            values[q] = f[f_idx(j, i, q, nx)];
        }
        values
    }

    #[inline]
    fn store_cell_populations(f: &mut [T], i: usize, j: usize, nx: usize, values: [T; 9]) {
        use crate::solvers::lbm::streaming::f_idx;

        for q in 0..9 {
            f[f_idx(j, i, q, nx)] = values[q];
        }
    }

    #[inline]
    fn reconstruct_zou_he_face(
        face: BoundaryFace,
        values: &mut [T; 9],
        density: T,
        velocity: [T; 2],
    ) {
        let half = T::from_f64(0.5).expect("0.5 is an exact f64 constant");
        let one_sixth = T::from_f64(1.0 / 6.0).expect("1/6 is an exact f64 constant");
        let two_thirds = T::from_f64(2.0 / 3.0).expect("2/3 is an exact f64 constant");
        let u_x = velocity[0];
        let u_y = velocity[1];

        match face {
            BoundaryFace::West => {
                values[1] = values[3] + two_thirds * density * u_x;
                values[5] = values[7]
                    + half * (values[4] - values[2])
                    + one_sixth * density * u_x
                    + half * density * u_y;
                values[8] = values[6] + half * (values[2] - values[4]) + one_sixth * density * u_x
                    - half * density * u_y;
            }
            BoundaryFace::East => {
                values[3] = values[1] - two_thirds * density * u_x;
                values[6] = values[8] + half * (values[4] - values[2]) - one_sixth * density * u_x
                    + half * density * u_y;
                values[7] = values[5] + half * (values[2] - values[4])
                    - one_sixth * density * u_x
                    - half * density * u_y;
            }
            BoundaryFace::South => {
                values[2] = values[4] + two_thirds * density * u_y;
                values[5] = values[7]
                    + half * (values[3] - values[1])
                    + one_sixth * density * u_y
                    + half * density * u_x;
                values[6] = values[8] + half * (values[1] - values[3]) + one_sixth * density * u_y
                    - half * density * u_x;
            }
            BoundaryFace::North => {
                values[4] = values[2] - two_thirds * density * u_y;
                values[7] = values[5] + half * (values[1] - values[3])
                    - one_sixth * density * u_y
                    - half * density * u_x;
                values[8] = values[6] + half * (values[3] - values[1]) - one_sixth * density * u_y
                    + half * density * u_x;
            }
        }
    }

    #[inline]
    fn boundary_density_from_velocity(face: BoundaryFace, values: &[T; 9], velocity: [T; 2]) -> T {
        match face {
            BoundaryFace::West => {
                let known = values[0]
                    + values[2]
                    + values[4]
                    + T::from_f64(2.0).expect("2 is exact") * (values[3] + values[6] + values[7]);
                known / (T::one() - velocity[0])
            }
            BoundaryFace::East => {
                let known = values[0]
                    + values[2]
                    + values[4]
                    + T::from_f64(2.0).expect("2 is exact") * (values[1] + values[5] + values[8]);
                known / (T::one() + velocity[0])
            }
            BoundaryFace::South => {
                let known = values[0]
                    + values[1]
                    + values[3]
                    + T::from_f64(2.0).expect("2 is exact") * (values[4] + values[7] + values[8]);
                known / (T::one() - velocity[1])
            }
            BoundaryFace::North => {
                let known = values[0]
                    + values[1]
                    + values[3]
                    + T::from_f64(2.0).expect("2 is exact") * (values[2] + values[5] + values[6]);
                known / (T::one() + velocity[1])
            }
        }
    }

    #[inline]
    fn boundary_velocity_from_density(face: BoundaryFace, values: &[T; 9], density: T) -> [T; 2] {
        match face {
            BoundaryFace::West => {
                let known = values[0]
                    + values[2]
                    + values[4]
                    + T::from_f64(2.0).expect("2 is exact") * (values[3] + values[6] + values[7]);
                [T::one() - known / density, T::zero()]
            }
            BoundaryFace::East => {
                let known = values[0]
                    + values[2]
                    + values[4]
                    + T::from_f64(2.0).expect("2 is exact") * (values[1] + values[5] + values[8]);
                [known / density - T::one(), T::zero()]
            }
            BoundaryFace::South => {
                let known = values[0]
                    + values[1]
                    + values[3]
                    + T::from_f64(2.0).expect("2 is exact") * (values[4] + values[7] + values[8]);
                [T::zero(), T::one() - known / density]
            }
            BoundaryFace::North => {
                let known = values[0]
                    + values[1]
                    + values[3]
                    + T::from_f64(2.0).expect("2 is exact") * (values[2] + values[5] + values[6]);
                [T::zero(), known / density - T::one()]
            }
        }
    }

    /// Apply bounce-back (no-slip wall) at node (i, j).
    ///
    /// Inverts each distribution with its antipodal: f_q ← f_{q̄}.
    /// Yields u=0 at the node (Theorem — Bounce-Back).
    pub fn apply_bounce_back(f: &mut [T], i: usize, j: usize, nx: usize) {
        use crate::solvers::lbm::streaming::f_idx;
        // Snapshot the 9 values before overwriting
        let mut f_buf = [T::zero(); 9];
        for q in 0..9 {
            f_buf[q] = f[f_idx(j, i, q, nx)];
        }
        for q in 0..9 {
            let q_opp = D2Q9::OPPOSITE[q];
            f[f_idx(j, i, q, nx)] = f_buf[q_opp];
        }
    }

    /// Apply velocity (Zou-He) boundary condition at node (i, j).
    ///
    /// Sets the prescribed boundary velocity, computes the matching density
    /// from the known populations on that face, and reconstructs the missing
    /// outgoing populations with the exact Zou-He closure.
    pub fn apply_velocity_boundary(
        f: &mut [T],
        density: &mut [T],
        velocity: &mut [T],
        i: usize,
        j: usize,
        nx: usize,
        ny: usize,
        u_boundary: [T; 2],
    ) {
        let cell = j * nx + i;
        let Some(face) = Self::boundary_face(i, j, nx, ny, u_boundary) else {
            debug_assert!(false, "Velocity boundary must lie on the domain boundary");
            return;
        };

        let mut cell_values = Self::load_cell_populations(f, i, j, nx);
        let rho = Self::boundary_density_from_velocity(face, &cell_values, u_boundary);
        Self::reconstruct_zou_he_face(face, &mut cell_values, rho, u_boundary);

        density[cell] = rho;
        velocity[cell * 2] = u_boundary[0];
        velocity[cell * 2 + 1] = u_boundary[1];
        Self::store_cell_populations(f, i, j, nx, cell_values);
    }

    /// Apply pressure boundary condition at node (i, j).
    ///
    /// Converts pressure to density via $\rho = p / c_s^2$, derives the
    /// face-normal velocity from the density constraint, preserves the
    /// tangential component from the adjacent interior face, and reconstructs
    /// the missing populations with the Zou-He closure.
    pub fn apply_pressure_boundary(
        f: &mut [T],
        density: &mut [T],
        velocity: &mut [T],
        i: usize,
        j: usize,
        nx: usize,
        ny: usize,
        p_boundary: T,
    ) {
        let cell = j * nx + i;

        let cs2 = T::from_f64(crate::constants::physics::LATTICE_SOUND_SPEED_SQUARED)
            .expect("cs² is an exact f64 constant");
        let rho = p_boundary / cs2;
        let mut cell_values = Self::load_cell_populations(f, i, j, nx);
        let velocity_hint = Self::extrapolate_velocity_flat(velocity, nx, ny, i, j);
        let Some(face) = Self::boundary_face(i, j, nx, ny, velocity_hint) else {
            debug_assert!(false, "Pressure boundary must lie on the domain boundary");
            return;
        };

        let mut u = Self::extrapolate_velocity_for_face(velocity, nx, ny, i, j, face);
        let normal_velocity = Self::boundary_velocity_from_density(face, &cell_values, rho);
        match face {
            BoundaryFace::West | BoundaryFace::East => {
                u[0] = normal_velocity[0];
            }
            BoundaryFace::South | BoundaryFace::North => {
                u[1] = normal_velocity[1];
            }
        }

        Self::reconstruct_zou_he_face(face, &mut cell_values, rho, u);

        density[cell] = rho;
        velocity[cell * 2] = u[0];
        velocity[cell * 2 + 1] = u[1];
        Self::store_cell_populations(f, i, j, nx, cell_values);
    }

    /// Apply all boundary conditions from the boundary map.
    pub fn apply_boundaries(
        &self,
        f: &mut [T],
        density: &mut [T],
        velocity: &mut [T],
        boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>,
        nx: usize,
        ny: usize,
    ) {
        for ((i, j), bc) in boundaries {
            match bc {
                BoundaryCondition::Wall { .. } => {
                    Self::apply_bounce_back(f, *i, *j, nx);
                }
                BoundaryCondition::VelocityInlet { velocity: vel } => {
                    let u = [vel[0], vel[1]];
                    Self::apply_velocity_boundary(f, density, velocity, *i, *j, nx, ny, u);
                }
                BoundaryCondition::PressureInlet { pressure, .. }
                | BoundaryCondition::PressureOutlet { pressure } => {
                    Self::apply_pressure_boundary(f, density, velocity, *i, *j, nx, ny, *pressure);
                }
                _ => {}
            }
        }
    }

    /// Compute density sum at node from flat buffer (used internally).
    fn _compute_boundary_density(f: &[T], i: usize, j: usize, nx: usize) -> T {
        use crate::solvers::lbm::streaming::f_idx;
        (0..9).fold(T::zero(), |acc, q| acc + f[f_idx(j, i, q, nx)])
    }

    /// First-order velocity extrapolation from interior at boundary node (i, j).
    fn extrapolate_velocity_flat(
        velocity: &[T],
        nx: usize,
        ny: usize,
        i: usize,
        j: usize,
    ) -> [T; 2] {
        let cell = |jj: usize, ii: usize| {
            [
                (velocity[(jj * nx + ii) * 2]),
                (velocity[(jj * nx + ii) * 2 + 1]),
            ]
        };
        if i == 0 && i + 1 < nx {
            cell(j, i + 1)
        } else if i + 1 == nx && i > 0 {
            cell(j, i - 1)
        } else if j == 0 && j + 1 < ny {
            cell(j + 1, i)
        } else if j + 1 == ny && j > 0 {
            cell(j - 1, i)
        } else {
            [T::zero(), T::zero()]
        }
    }

    fn extrapolate_velocity_for_face(
        velocity: &[T],
        nx: usize,
        ny: usize,
        i: usize,
        j: usize,
        face: BoundaryFace,
    ) -> [T; 2] {
        let cell = |jj: usize, ii: usize| {
            [
                velocity[(jj * nx + ii) * 2],
                velocity[(jj * nx + ii) * 2 + 1],
            ]
        };

        match face {
            BoundaryFace::West if i + 1 < nx => cell(j, i + 1),
            BoundaryFace::East if i > 0 => cell(j, i - 1),
            BoundaryFace::South if j + 1 < ny => cell(j + 1, i),
            BoundaryFace::North if j > 0 => cell(j - 1, i),
            _ => [T::zero(), T::zero()],
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    fn write_cell(f: &mut [f64], i: usize, j: usize, nx: usize, values: [f64; 9]) {
        use crate::solvers::lbm::streaming::f_idx;

        for q in 0..9 {
            f[f_idx(j, i, q, nx)] = values[q];
        }
    }

    #[test]
    fn test_bounce_back() {
        // Theorem — Bounce-Back: after apply_bounce_back, f[q] = f_original[q̄]
        let nx = 10_usize;
        let ny = 10_usize;
        let mut f = vec![0.1_f64; nx * ny * 9];
        use crate::solvers::lbm::streaming::f_idx;

        // Set specific values at node (5,5)
        for q in 0..9 {
            f[f_idx(5, 5, q, nx)] = q as f64;
        }

        BoundaryHandler::<f64>::apply_bounce_back(&mut f, 5, 5, nx);

        for q in 0..9 {
            let q_opp = D2Q9::OPPOSITE[q];
            assert_eq!(
                f[f_idx(5, 5, q, nx)],
                q_opp as f64,
                "Bounce-back: f[{q}] should equal original f[{q_opp}]"
            );
        }
    }

    #[test]
    fn velocity_boundary_uses_zou_he_reconstruction() {
        let nx = 4_usize;
        let ny = 4_usize;
        let mut f = vec![0.0_f64; nx * ny * 9];
        let mut density = vec![0.0_f64; nx * ny];
        let mut velocity = vec![0.0_f64; nx * ny * 2];
        let cell = 1 * nx;
        let values = [0.31, 0.12, 0.23, 0.47, 0.19, 0.15, 0.17, 0.21, 0.29];
        write_cell(&mut f, 0, 1, nx, values);

        let u = [0.08_f64, -0.02_f64];
        BoundaryHandler::<f64>::apply_velocity_boundary(
            &mut f,
            &mut density,
            &mut velocity,
            0,
            1,
            nx,
            ny,
            u,
        );

        let rho_expected =
            (values[0] + values[2] + values[4] + 2.0 * (values[3] + values[6] + values[7]))
                / (1.0 - u[0]);
        let f1_expected = values[3] + (2.0 / 3.0) * rho_expected * u[0];
        let f5_expected = values[7]
            + 0.5 * (values[4] - values[2])
            + (1.0 / 6.0) * rho_expected * u[0]
            + 0.5 * rho_expected * u[1];
        let f8_expected =
            values[6] + 0.5 * (values[2] - values[4]) + (1.0 / 6.0) * rho_expected * u[0]
                - 0.5 * rho_expected * u[1];

        assert_relative_eq!(density[cell], rho_expected, epsilon = 1e-12);
        assert_relative_eq!(velocity[cell * 2], u[0], epsilon = 1e-12);
        assert_relative_eq!(velocity[cell * 2 + 1], u[1], epsilon = 1e-12);
        let west_base = (1 * nx + 0) * 9;
        assert_relative_eq!(f[west_base + 1], f1_expected, epsilon = 1e-12);
        assert_relative_eq!(f[west_base + 5], f5_expected, epsilon = 1e-12);
        assert_relative_eq!(f[west_base + 8], f8_expected, epsilon = 1e-12);
        assert_relative_eq!(f[west_base + 0], values[0], epsilon = 1e-12);
        assert_relative_eq!(f[west_base + 2], values[2], epsilon = 1e-12);
        assert_relative_eq!(f[west_base + 3], values[3], epsilon = 1e-12);
        assert_relative_eq!(f[west_base + 4], values[4], epsilon = 1e-12);
        assert_relative_eq!(f[west_base + 6], values[6], epsilon = 1e-12);
        assert_relative_eq!(f[west_base + 7], values[7], epsilon = 1e-12);
    }

    #[test]
    fn pressure_boundary_reconstructs_missing_face_populations() {
        let nx = 4_usize;
        let ny = 4_usize;
        let mut f = vec![0.0_f64; nx * ny * 9];
        let mut density = vec![0.0_f64; nx * ny];
        let mut velocity = vec![0.0_f64; nx * ny * 2];
        let values = [0.22, 0.34, 0.18, 0.29, 0.27, 0.11, 0.16, 0.25, 0.19];
        write_cell(&mut f, 3, 1, nx, values);
        let left_neighbor = 1 * nx + 2;
        velocity[left_neighbor * 2] = 0.07;
        velocity[left_neighbor * 2 + 1] = -0.04;

        let cs2 = crate::constants::physics::LATTICE_SOUND_SPEED_SQUARED;
        let rho = 1.24_f64;
        let pressure = rho * cs2;
        BoundaryHandler::<f64>::apply_pressure_boundary(
            &mut f,
            &mut density,
            &mut velocity,
            3,
            1,
            nx,
            ny,
            pressure,
        );

        let velocity_hint = [0.07_f64, -0.04_f64];
        let known = values[0] + values[2] + values[4] + 2.0 * (values[1] + values[5] + values[8]);
        let u_x = known / rho - 1.0;
        let expected_velocity = [u_x, velocity_hint[1]];
        let f3_expected = values[1] - (2.0 / 3.0) * rho * u_x;
        let f6_expected = values[8] + 0.5 * (values[4] - values[2]) - (1.0 / 6.0) * rho * u_x
            + 0.5 * rho * velocity_hint[1];
        let f7_expected = values[5] + 0.5 * (values[2] - values[4])
            - (1.0 / 6.0) * rho * u_x
            - 0.5 * rho * velocity_hint[1];

        assert_relative_eq!(density[1 * nx + 3], rho, epsilon = 1e-12);
        assert_relative_eq!(
            velocity[(1 * nx + 3) * 2],
            expected_velocity[0],
            epsilon = 1e-12
        );
        assert_relative_eq!(
            velocity[(1 * nx + 3) * 2 + 1],
            expected_velocity[1],
            epsilon = 1e-12
        );
        let east_base = (1 * nx + 3) * 9;
        assert_relative_eq!(f[east_base + 3], f3_expected, epsilon = 1e-12);
        assert_relative_eq!(f[east_base + 6], f6_expected, epsilon = 1e-12);
        assert_relative_eq!(f[east_base + 7], f7_expected, epsilon = 1e-12);
    }
}
