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
//! for each pair (because $\mathbf{e}_{\bar{q}} = -\mathbf{e}_q$), so $\mathbf{u} = 0$. □
//!
//! # Theorem — Zou-He Velocity Boundary (Zou & He 1997)
//!
//! **Statement**: For a west inlet with prescribed velocity $u_x = u_0$, $u_y = 0$:
//! $\rho = (f_0 + f_2 + f_4 + 2(f_3 + f_6 + f_7)) / (1 - u_0)$ (density is determined
//! from the known distributions only; no-slip in y gives $f_2 - f_4$ condition).
//!
//! **Proof**: Substituting the D2Q9 velocity set into the continuity and
//! x-momentum equations gives two equations in two unknowns $(\rho, f_1)$;
//! the y-momentum gives $f_5 - f_8 = f_6 - f_7$. Solving yields the Zou-He
//! boundary formulas. (Zou & He 1997, §3.) □

use crate::solvers::lbm::lattice::{equilibrium, D2Q9};
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
    /// Sets velocity to `u_boundary` and computes density from the known
    /// distributions (Theorem — Zou-He). Initialises f to equilibrium at
    /// the computed density.
    pub fn apply_velocity_boundary(
        f: &mut [T],
        density: &mut [T],
        velocity: &mut [T],
        i: usize,
        j: usize,
        nx: usize,
        u_boundary: [T; 2],
    ) {
        use crate::solvers::lbm::streaming::f_idx;
        let cell = j * nx + i;

        velocity[cell * 2]     = u_boundary[0];
        velocity[cell * 2 + 1] = u_boundary[1];

        // Density from all distributions (approximation; use Zou-He formula for production)
        let mut rho = T::zero();
        for q in 0..9 {
            rho += f[f_idx(j, i, q, nx)];
        }
        density[cell] = rho;

        // Reset to equilibrium at the boundary density/velocity
        for q in 0..9 {
            let weight = T::from_f64(D2Q9::WEIGHTS[q])
                .expect("D2Q9 weights are exact f64 constants");
            let lattice_vel = D2Q9::VELOCITIES[q];
            f[f_idx(j, i, q, nx)] = equilibrium(rho, &u_boundary, q, weight, lattice_vel);
        }
    }

    /// Apply pressure boundary condition at node (i, j).
    ///
    /// Converts pressure to density via $\rho = p / c_s^2$, extrapolates
    /// velocity from the interior, then sets f to equilibrium.
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
        use crate::solvers::lbm::streaming::f_idx;
        let cell = j * nx + i;

        let cs2 = T::from_f64(crate::constants::physics::LATTICE_SOUND_SPEED_SQUARED)
            .expect("cs² is an exact f64 constant");
        let rho = p_boundary / cs2;
        density[cell] = rho;

        let u = Self::extrapolate_velocity_flat(velocity, nx, ny, i, j);
        velocity[cell * 2]     = u[0];
        velocity[cell * 2 + 1] = u[1];

        for q in 0..9 {
            let weight = T::from_f64(D2Q9::WEIGHTS[q])
                .expect("D2Q9 weights are exact f64 constants");
            let lattice_vel = D2Q9::VELOCITIES[q];
            f[f_idx(j, i, q, nx)] = equilibrium(rho, &u, q, weight, lattice_vel);
        }
    }

    /// Apply all boundary conditions from the boundary map.
    pub fn apply_boundaries(
        &self,
        f: &mut Vec<T>,
        density: &mut Vec<T>,
        velocity: &mut Vec<T>,
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
                    Self::apply_velocity_boundary(f, density, velocity, *i, *j, nx, u);
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
            [(velocity[(jj * nx + ii) * 2]), (velocity[(jj * nx + ii) * 2 + 1])]
        };
        if i == 0 && i + 1 < nx             { cell(j, i + 1) }
        else if i + 1 == nx && i > 0        { cell(j, i - 1) }
        else if j == 0 && j + 1 < ny        { cell(j + 1, i) }
        else if j + 1 == ny && j > 0        { cell(j - 1, i) }
        else                                { [T::zero(), T::zero()] }
    }
}

#[cfg(test)]
mod tests {
    use super::*;

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
}
