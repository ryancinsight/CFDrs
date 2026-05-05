//! Main LBM solver — integrates collision, streaming, boundary, and macroscopic.
//!
//! # Theorem — Chapman-Enskog Expansion (2nd order)
//!
//! **Statement**: The BGK-LBM with D2Q9 lattice, $c_s^2 = 1/3$, and relaxation
//! time $\tau$ recovers the incompressible 2D Navier-Stokes equations:
//!
//! $$\nabla \cdot \mathbf{u} = 0, \qquad
//!   \frac{\partial \mathbf{u}}{\partial t} + (\mathbf{u} \cdot \nabla)\mathbf{u}
//!   = -\nabla p / \rho + \nu \nabla^2 \mathbf{u}$$
//!
//! with $\nu = c_s^2 (\tau - \tfrac{1}{2}) \Delta t$ in the low Mach number limit.
//!
//! **Proof**:
//!
//! 1. *Expansion*: Write $f_i = f_i^{(0)} + \epsilon f_i^{(1)} + \epsilon^2 f_i^{(2)}$
//!    where $\epsilon = Kn$ (Knudsen number), and expand time and space derivatives
//!    as $\partial_t = \epsilon \partial_{t_1} + \epsilon^2 \partial_{t_2}$,
//!    $\nabla = \epsilon \nabla_1$.
//!
//! 2. *O(ε⁰)*: $f_i^{(0)} = f_i^{eq}$. Moments give the continuity equation
//!    $\partial_{t_1} \rho + \nabla_1 \cdot (\rho \mathbf{u}) = 0$.
//!
//! 3. *O(ε¹)*: The deviation $f_i^{(1)} = -\tau(\partial_{t_1} + \mathbf{e}_i \cdot \nabla_1) f_i^{eq}$.
//!    Taking the second moment yields the viscous stress tensor with
//!    $\mu = \rho c_s^2 (\tau - \tfrac{1}{2}) \Delta t$.
//!
//! 4. *Momentum equation*: Adding the $O(\epsilon^2)$ time derivative gives the full
//!    Navier-Stokes equation with $\nu = \mu/\rho = c_s^2 (\tau - \tfrac{1}{2}) \Delta t$. □
//!
//! **Reference**: He & Luo (1997), *Phys. Rev. E* 56, 6811; Succi (2001), §4.3.

use crate::grid::{Grid2D, StructuredGrid2D};
use crate::solvers::lbm::{
    boundary::BoundaryHandler,
    collision::{BgkCollision, CollisionOperator},
    lattice::{equilibrium, D2Q9},
    macroscopic::{compute_pressure, MacroscopicQuantities},
    streaming::{f_idx, StreamingOperator},
};
use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use std::collections::HashMap;

/// Configuration for the LBM solver.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct LbmConfig<T: RealField + Copy> {
    /// Relaxation time τ ∈ (0.5, ∞).  Stability requires τ > 0.5.
    pub tau: T,
    /// Maximum number of time steps.
    pub max_steps: usize,
    /// Convergence tolerance for ‖Δu‖_∞.
    pub tolerance: T,
    /// Steps between convergence checks and (optional) output.
    pub output_frequency: usize,
    /// Enable verbose stdout progress.
    pub verbose: bool,
}

impl<T: RealField + Copy + FromPrimitive> Default for LbmConfig<T> {
    fn default() -> Self {
        Self {
            tau: T::from_f64(1.0).expect("T must represent f64; τ = 1.0"),
            max_steps: 10_000,
            tolerance: T::from_f64(1e-6).expect("T must represent f64; tol = 1e-6"),
            output_frequency: 100,
            verbose: false,
        }
    }
}

/// Lattice Boltzmann Method solver for 2D incompressible flows.
///
/// # Memory Layout
///
/// All distribution functions are stored in a single flat `Vec<T>`:
/// ```text
/// f[j * nx * 9 + i * 9 + q]
/// ```
/// This achieves stride-1 access over the 9 velocity directions at each node —
/// the innermost loop in collision and macroscopic computations (Theorem:
/// reduces expected cache misses by factor O(N_y)).
///
/// Macroscopic density and velocity are stored in separate flat `Vec<T>` buffers:
/// ```text
/// density[j * nx + i]
/// velocity[(j * nx + i) * 2 + d]   // d=0→x, d=1→y
/// ```
pub struct LbmSolver<T: RealField + Copy> {
    config: LbmConfig<T>,
    /// Flat distribution buffer (layout: j*nx*9 + i*9 + q)
    f: Vec<T>,
    /// Streaming double-buffer (pre-allocated, swapped each step)
    f_buffer: Vec<T>,
    /// Macroscopic quantities
    macroscopic: MacroscopicQuantities<T>,
    /// Collision operator (boxed for polymorphism)
    collision: Box<dyn CollisionOperator<T>>,
    /// Boundary condition handler
    boundary_handler: BoundaryHandler<T>,
    nx: usize,
    ny: usize,
    dx: T,
    dy: T,
    step_count: usize,
    /// Previous-step velocity buffer for convergence checking (pre-allocated)
    previous_velocity: Vec<T>,
    /// Optional flat distribution buffer for passive scalar nuclei (layout: j*nx*9 + i*9 + q)
    g: Option<Vec<T>>,
    /// Optional streaming double-buffer for nuclei
    g_buffer: Option<Vec<T>>,
    /// Optional nuclei transport physics evaluator
    nuclei_transport: Option<cfd_core::physics::cavitation::nuclei_transport::NucleiTransport<T>>,
    /// Relaxation time for the scalar lattice $\tau_g$
    tau_g: Option<T>,
}

impl<T: RealField + Copy + FromPrimitive> LbmSolver<T>
where
    T: Send + Sync + std::fmt::LowerExp,
{
    /// Construct from config and grid.
    #[must_use]
    pub fn new(config: LbmConfig<T>, grid: &StructuredGrid2D<T>) -> Self {
        let nx = grid.nx();
        let ny = grid.ny();
        let n = nx * ny;

        // Single flat allocation for all distribution functions
        let f = vec![T::zero(); n * 9];
        let f_buffer = vec![T::zero(); n * 9];

        let macroscopic = MacroscopicQuantities::new(nx, ny);
        let collision = Box::new(BgkCollision::new(config.tau));
        let boundary_handler = BoundaryHandler::new();
        let previous_velocity = vec![T::zero(); n * 2];

        Self {
            config,
            f,
            f_buffer,
            macroscopic,
            collision,
            boundary_handler,
            nx,
            ny,
            dx: grid.dx,
            dy: grid.dy,
            step_count: 0,
            previous_velocity,
            g: None,
            g_buffer: None,
            nuclei_transport: None,
            tau_g: None,
        }
    }

    #[inline]
    fn lattice_cavitation_source(density: T) -> T {
        let lattice_pressure = compute_pressure(density);
        let lattice_vapor_threshold = compute_pressure(T::one());
        if lattice_pressure < lattice_vapor_threshold {
            lattice_vapor_threshold - lattice_pressure
        } else {
            T::zero()
        }
    }

    /// Enable passive scalar advection-diffusion for nuclei transport.
    #[must_use]
    pub fn with_nuclei_transport(
        mut self,
        config: cfd_core::physics::cavitation::nuclei_transport::NucleiTransportConfig<T>,
    ) -> Self {
        let n = self.nx * self.ny;
        self.g = Some(vec![T::zero(); n * 9]);
        self.g_buffer = Some(vec![T::zero(); n * 9]);
        self.macroscopic = self.macroscopic.with_nuclei();
        self.nuclei_transport =
            Some(cfd_core::physics::cavitation::nuclei_transport::NucleiTransport::new(config));

        // Assuming dt = 1 in lattice units for this solver's inner loop if not specified
        self.tau_g = Some(self.config.tau);
        self
    }

    /// Compute the equilibrium distribution at given density and velocity.
    pub fn equilibrium_distribution(&self, density: T, velocity: Vector2<T>) -> Vec<T> {
        let u = [velocity.x, velocity.y];
        (0..9)
            .map(|q| {
                let weight =
                    T::from_f64(D2Q9::WEIGHTS[q]).expect("D2Q9 weights are exact f64 constants");
                equilibrium(density, &u, q, weight, D2Q9::VELOCITIES[q])
            })
            .collect()
    }

    #[inline]
    fn validate_low_mach_velocity(velocity: Vector2<T>) -> Result<()> {
        let cs = T::from_f64(crate::constants::physics::LATTICE_SOUND_SPEED_SQUARED)
            .expect("cs² is an exact f64 constant")
            .sqrt();
        let mach_limit = T::from_f64(0.1).expect("0.1 is an exact f64 constant");
        let speed = (velocity.x * velocity.x + velocity.y * velocity.y).sqrt();

        if speed / cs > mach_limit {
            return Err(Error::InvalidConfiguration(
                "LBM initialization violates Ma <= 0.1 low-Mach incompressible limit".to_string(),
            ));
        }
        Ok(())
    }

    /// Get macroscopic density and velocity at node (i, j).
    pub fn compute_macroscopic(&self, i: usize, j: usize) -> (T, Vector2<T>) {
        let rho = self.macroscopic.density_at(i, j);
        let [ux, uy] = self.macroscopic.velocity_at(i, j);
        (rho, Vector2::new(ux, uy))
    }

    /// Initialise distribution functions using user-provided density and velocity fields.
    ///
    /// Sets all f_q(i,j) = f_q^eq(ρ(x,y), **u**(x,y)).
    pub fn initialize<F1, F2>(&mut self, density_fn: F1, velocity_fn: F2) -> Result<()>
    where
        F1: Fn(T, T) -> T,
        F2: Fn(T, T) -> Vector2<T>,
    {
        let nx = self.nx;
        let ny = self.ny;

        for j in 0..ny {
            for i in 0..nx {
                let x = T::from_usize(i).expect("grid index fits in T") * self.dx;
                let y = T::from_usize(j).expect("grid index fits in T") * self.dy;

                let rho = density_fn(x, y);
                let vel = velocity_fn(x, y);
                Self::validate_low_mach_velocity(vel)?;
                let u = [vel.x, vel.y];

                let cell = j * nx + i;
                self.macroscopic.density[cell] = rho;
                self.macroscopic.velocity[cell * 2] = u[0];
                self.macroscopic.velocity[cell * 2 + 1] = u[1];

                for q in 0..9 {
                    let weight = T::from_f64(D2Q9::WEIGHTS[q])
                        .expect("D2Q9 weights are exact f64 constants");
                    self.f[f_idx(j, i, q, nx)] =
                        equilibrium(rho, &u, q, weight, D2Q9::VELOCITIES[q]);
                }
            }
        }

        self.step_count = 0;
        Ok(())
    }

    /// Perform one time step: macroscopic → collision → streaming → boundary.
    pub fn step(
        &mut self,
        boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>,
    ) -> Result<()> {
        let nx = self.nx;
        let ny = self.ny;

        // 1. Update macroscopic fields from current f (and g)
        self.macroscopic
            .update_from_distributions(&self.f, self.g.as_deref());

        // 2. Collision step (operates in-place on self.f)
        self.collision.collide(
            &mut self.f,
            &self.macroscopic.density,
            &self.macroscopic.velocity,
            nx,
            ny,
        );

        // 3. Streaming step (pull scheme, double-buffer swap — zero allocation)
        StreamingOperator::stream(&self.f, &mut self.f_buffer, nx, ny);
        std::mem::swap(&mut self.f, &mut self.f_buffer);

        // 3b. Collision and streaming for passive scalar `g`
        if let (Some(g), Some(g_buf), Some(nuclei_fraction), Some(transport), Some(tau_g)) = (
            &mut self.g,
            &mut self.g_buffer,
            &self.macroscopic.nuclei_fraction,
            &self.nuclei_transport,
            self.tau_g,
        ) {
            let omega_g = T::one() / tau_g;
            for j in 0..ny {
                for i in 0..nx {
                    let cell = j * nx + i;
                    let u = [
                        self.macroscopic.velocity[cell * 2],
                        self.macroscopic.velocity[cell * 2 + 1],
                    ];
                    let phi = nuclei_fraction[cell];

                    // Lattice-space cavitation source: once the local pressure drops below the
                    // unit-density equilibrium threshold, seed nuclei generation proportionally
                    // to the pressure deficit.
                    let macroscopic_source =
                        Self::lattice_cavitation_source(self.macroscopic.density[cell]);

                    let s_net = transport.calculate_net_reaction_rate(phi, macroscopic_source);

                    for q in 0..9 {
                        let weight = T::from_f64(D2Q9::WEIGHTS[q])
                            .expect("D2Q9 weights are exact f64 constants");
                        let g_eq = equilibrium(phi, &u, q, weight, D2Q9::VELOCITIES[q]);
                        let idx = f_idx(j, i, q, nx);
                        // Advection-diffusion collision + source/sink
                        g[idx] = g[idx] - omega_g * (g[idx] - g_eq) + weight * s_net;
                    }
                }
            }
            StreamingOperator::stream(g, g_buf, nx, ny);
            std::mem::swap(g, g_buf);
        }

        // 4. Apply boundary conditions
        self.boundary_handler.apply_boundaries(
            &mut self.f,
            &mut self.macroscopic.density,
            &mut self.macroscopic.velocity,
            boundaries,
            nx,
            ny,
        )?;

        if let Some(g) = &mut self.g {
            super::scalar_boundary::apply_scalar_boundaries(g, boundaries, nx, ny);
        }

        self.step_count += 1;
        Ok(())
    }
    /// Run the solver until convergence (‖Δu‖_∞ < tol) or max_steps.
    ///
    /// Non-convergence does not return an error — the caller inspects the
    /// returned velocity field and decides. This is consistent with the
    /// LBM design contract: steady-state is a user concern.
    pub fn solve(
        &mut self,
        boundaries: HashMap<(usize, usize), BoundaryCondition<T>>,
        initial_density: T,
        initial_velocity: Vector2<T>,
    ) -> Result<()> {
        let v = initial_velocity;
        self.initialize(|_, _| initial_density, |_, _| v)?;

        let mut converged = false;
        self.copy_velocity_to_buffer();

        for step in 0..self.config.max_steps {
            self.step(&boundaries)?;

            if step % self.config.output_frequency == 0 {
                let max_change = self.compute_max_velocity_change();

                if self.config.verbose {
                    tracing::debug!("Step {step}: max velocity change = {max_change:e}");
                }

                if max_change < self.config.tolerance {
                    converged = true;
                    if self.config.verbose {
                        tracing::debug!("Converged after {step} steps");
                    }
                    break;
                }

                // Update reference snapshot for next comparison (reuses pre-allocated buffer)
                self.copy_velocity_to_buffer();
            }
        }

        if !converged && self.config.verbose {
            tracing::debug!(
                "Warning: Did not converge after {} steps",
                self.config.max_steps
            );
        }

        Ok(())
    }

    /// Copy velocity buffer for convergence checking (zero additional allocation).
    fn copy_velocity_to_buffer(&mut self) {
        self.previous_velocity
            .copy_from_slice(&self.macroscopic.velocity);
    }

    /// Compute ‖u^{n+1} − u^n‖_∞ for convergence check.
    fn compute_max_velocity_change(&self) -> T {
        self.macroscopic
            .velocity
            .iter()
            .zip(self.previous_velocity.iter())
            .fold(T::zero(), |acc, (&cur, &prev)| acc.max((cur - prev).abs()))
    }

    /// Get the flat velocity field slice.
    pub fn velocity_field(&self) -> &[T] {
        &self.macroscopic.velocity
    }

    /// Get the flat density field slice.
    pub fn density_field(&self) -> &[T] {
        &self.macroscopic.density
    }

    /// Get velocity at node (i, j).
    pub fn velocity_at(&self, i: usize, j: usize) -> Option<[T; 2]> {
        if i < self.nx && j < self.ny {
            Some(self.macroscopic.velocity_at(i, j))
        } else {
            None
        }
    }

    /// Get density at node (i, j).
    pub fn density_at(&self, i: usize, j: usize) -> Option<T> {
        if i < self.nx && j < self.ny {
            Some(self.macroscopic.density_at(i, j))
        } else {
            None
        }
    }

    /// Get macroscopic field slices (velocity, density).
    pub fn get_macroscopic(&self) -> (&[T], &[T]) {
        (&self.macroscopic.velocity, &self.macroscopic.density)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::grid::StructuredGrid2D;
    use approx::assert_relative_eq;

    #[test]
    fn test_equilibrium_distribution() -> Result<()> {
        // Theorem — Moment Consistency: ∑_q f_q^eq = ρ exactly.
        let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)?;
        let config = LbmConfig::<f64>::default();
        let solver = LbmSolver::new(config, &grid);

        let rho = 1.0;
        let u = Vector2::new(0.1, 0.0);
        let feq = solver.equilibrium_distribution(rho, u);

        let sum: f64 = feq.iter().sum();
        assert_relative_eq!(sum, rho, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn test_initialization() -> Result<()> {
        let grid = StructuredGrid2D::<f64>::new(10, 10, 0.0, 1.0, 0.0, 1.0)?;
        let config = LbmConfig::<f64>::default();
        let mut solver = LbmSolver::new(config, &grid);

        solver.initialize(|_, _| 1.0, |_, _| Vector2::new(0.0, 0.0))?;

        let (rho, u) = solver.compute_macroscopic(5, 5);
        assert_relative_eq!(rho, 1.0, epsilon = 1e-10);
        assert_relative_eq!(u.x, 0.0, epsilon = 1e-10);
        assert_relative_eq!(u.y, 0.0, epsilon = 1e-10);
        Ok(())
    }

    #[test]
    fn initialization_rejects_high_mach_velocity() -> Result<()> {
        let grid = StructuredGrid2D::<f64>::new(4, 4, 0.0, 1.0, 0.0, 1.0)?;
        let config = LbmConfig::<f64>::default();
        let mut solver = LbmSolver::new(config, &grid);

        let result = solver.initialize(|_, _| 1.0, |_, _| Vector2::new(0.2, 0.0));

        assert!(result.is_err(), "D2Q9 initialization must reject Ma > 0.1");
        Ok(())
    }

    #[test]
    fn test_total_mass_conservation_one_step() -> Result<()> {
        // Chapman-Enskog: mass is conserved ∑_{i,j} ρ(i,j) = const.
        let nx = 8_usize;
        let ny = 8_usize;
        let grid = StructuredGrid2D::<f64>::new(nx, ny, 0.0, 1.0, 0.0, 1.0)?;
        let config = LbmConfig::<f64>::default();
        let mut solver = LbmSolver::new(config, &grid);

        solver.initialize(|_, _| 1.0, |_, _| Vector2::new(0.05, 0.0))?;

        // Update macroscopic once to populate density field
        solver
            .macroscopic
            .update_from_distributions(&solver.f, None);
        let mass_before: f64 = solver.macroscopic.density.iter().sum();

        solver.step(&HashMap::new())?;
        solver
            .macroscopic
            .update_from_distributions(&solver.f, None);
        let mass_after: f64 = solver.macroscopic.density.iter().sum();

        assert_relative_eq!(mass_before, mass_after, epsilon = 1e-8);
        Ok(())
    }

    #[test]
    fn test_lattice_cavitation_source_tracks_pressure_deficit() -> Result<()> {
        let grid = StructuredGrid2D::<f64>::new(4, 4, 0.0, 1.0, 0.0, 1.0)?;
        let config = LbmConfig::<f64>::default();
        let _solver = LbmSolver::new(config, &grid);

        let below_threshold = LbmSolver::lattice_cavitation_source(0.92);
        let at_threshold = LbmSolver::lattice_cavitation_source(1.0);
        let above_threshold = LbmSolver::lattice_cavitation_source(1.08);

        assert!(below_threshold > 0.0);
        assert_relative_eq!(at_threshold, 0.0, epsilon = 1e-15);
        assert_relative_eq!(above_threshold, 0.0, epsilon = 1e-15);
        Ok(())
    }
}
