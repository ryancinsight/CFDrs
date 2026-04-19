//! # Level Set Method for 3D Interface Tracking
//!
//! This module implements the Level Set method for accurate tracking of evolving
//! interfaces in multiphase flows using signed distance functions.
//!
//! ## Mathematical Foundation
//!
//! ### Signed Distance Function Theorem
//!
//! **Statement**: A signed distance function φ(x,t) represents the interface Γ(t)
//! as its zero level set, with φ > 0 outside the interface and φ < 0 inside.
//!
//! **Mathematical Properties**:
//!
//! ```math
//! φ(x,t) = ±d(x, Γ(t))    where d is Euclidean distance
//! |∇φ| = 1                 (unit gradient magnitude — Eikonal equation)
//! φ(Γ(t), t) = 0           (zero level set condition)
//! ```
//!
//! **Reference**: Osher, S., Fedkiw, R. (2003). *Level Set Methods and Dynamic
//! Implicit Surfaces*. Cambridge University Press.
//!
//! ### Level Set Evolution Theorem
//!
//! **Statement**: The evolution of interfaces in a velocity field u(x,t) satisfies
//! the Hamilton-Jacobi advection equation:
//!
//! ```math
//! ∂φ/∂t + u·∇φ = 0
//! ```
//!
//! **Reference**: Sethian, J.A. (1999). *Level Set Methods and Fast Marching Methods*.
//! Cambridge University Press.
//!
//! ### WENO5-Z Stability Theorem (Borges et al. 2008)
//!
//! **Statement**: The 5th-order WENO-Z scheme applied to
//!
//! ```math
//! φ_t + H(∇φ) = 0
//! ```
//!
//! is (i) non-oscillatory near discontinuities, (ii) 5th-order accurate in smooth
//! regions, and (iii) stable under the combined CFL condition:
//!
//! ```math
//! CFL = max_i (|u_x|/Δx + |u_y|/Δy + |u_z|/Δz) Δt ≤ 0.5
//! ```
//!
//! when paired with at least a third-order Total Variation Diminishing (TVD)
//! Runge-Kutta time integrator.
//!
//! **Reference**: Borges, R., Carmona, M., Costa, B. & Don, W. (2008).
//! "An improved weighted essentially non-oscillatory scheme for hyperbolic
//! conservation laws." J. Comput. Phys. 227:3191-3211.
//!
//! ### Reinitialization Theorem (Sussman et al. 1994)
//!
//! **Statement**: The signed distance property |∇φ| = 1 is restored by iterating
//! the reinitialization equation until convergence:
//!
//! ```math
//! ∂φ/∂τ + S(φ₀)(|∇φ| − 1) = 0
//! ```
//!
//! where the smoothed sign function `S(φ₀) = φ₀ / √(φ₀² + Δx²)` ensures subcell
//! resolution of the zero level set position.
//!
//! **Convergence**: The iteration converges exponentially when `∂τ ≤ 0.5 min(Δx,Δy,Δz)`,
//! and the zero level set position is preserved to `O(Δx²)`.
//!
//! **Reference**: Sussman, M., Smereka, P. & Osher, S. (1994). "A level set approach
//! for computing solutions to incompressible two-phase flow."
//! J. Comput. Phys. 114:146-159.
//!
//! ### Algorithm
//!
//! 1. **Initialization**: Construct signed distance function from geometry.
//! 2. **Evolution**: Solve ∂φ/∂t + u·∇φ = 0 using WENO5-Z (one Euler step per advance).
//! 3. **Reinitialization**: Restore |∇φ| = 1 every `config.reinitialization_interval` steps
//!    using convergence-based iteration.
//! 4. **Narrow Band**: Update index set of cells within `band_width` grid spacings
//!    of the interface.

use super::{config::LevelSetConfig, weno::weno5_derivative};
use cfd_core::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::FromPrimitive;

/// Level Set solver for interface tracking
pub struct LevelSetSolver<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> {
    config: LevelSetConfig,
    /// Grid dimensions
    nx: usize,
    ny: usize,
    nz: usize,
    /// Grid spacing
    dx: T,
    dy: T,
    dz: T,
    /// Level set function (signed distance)
    phi: Vec<T>,
    /// Previous timestep / working buffer for double-buffering
    phi_previous: Vec<T>,
    /// Reinitialization reference buffer (stores φ₀ for S(φ₀) sign function)
    phi_reinit: Vec<T>,
    /// Velocity field
    velocity: Vec<Vector3<T>>,
    /// Narrow band indices (if using narrow band)
    narrow_band: Vec<usize>,
    /// Time step counter
    time_step: usize,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy> LevelSetSolver<T> {
    /// Create a new Level Set solver
    pub fn new(
        config: LevelSetConfig,
        nx: usize,
        ny: usize,
        nz: usize,
        dx: T,
        dy: T,
        dz: T,
    ) -> Self {
        let grid_size = nx * ny * nz;
        Self {
            config,
            nx,
            ny,
            nz,
            dx,
            dy,
            dz,
            phi: vec![T::zero(); grid_size],
            phi_previous: vec![T::zero(); grid_size],
            phi_reinit: vec![T::zero(); grid_size],
            velocity: vec![Vector3::zeros(); grid_size],
            narrow_band: Vec::new(),
            time_step: 0,
        }
    }

    /// Get grid index from (i, j, k) coordinates
    #[inline]
    pub fn index(&self, i: usize, j: usize, k: usize) -> usize {
        k * self.ny * self.nx + j * self.nx + i
    }

    /// Grid x-dimension
    pub fn nx(&self) -> usize {
        self.nx
    }
    /// Grid y-dimension
    pub fn ny(&self) -> usize {
        self.ny
    }
    /// Grid z-dimension
    pub fn nz(&self) -> usize {
        self.nz
    }

    /// Get the level set function
    pub fn phi(&self) -> &[T] {
        &self.phi
    }

    /// Get mutable reference to level set function
    pub fn phi_mut(&mut self) -> &mut [T] {
        &mut self.phi
    }

    /// Set velocity field
    pub fn set_velocity(&mut self, velocity: Vec<Vector3<T>>) {
        self.velocity = velocity;
    }

    /// Get narrow band indices
    pub fn narrow_band(&self) -> &[usize] {
        &self.narrow_band
    }

    /// Update narrow band indices based on current level set
    pub fn update_narrow_band(&mut self) {
        self.narrow_band.clear();
        let band_width = <T as FromPrimitive>::from_f64(self.config.band_width)
            .expect("band_width config is an IEEE 754 representable f64");

        for idx in 0..self.phi.len() {
            if num_traits::Float::abs(self.phi[idx])
                <= band_width
                    * num_traits::Float::min(num_traits::Float::min(self.dx, self.dy), self.dz)
            {
                self.narrow_band.push(idx);
            }
        }
    }

    /// Advance level set by one time step using WENO5-Z advection.
    pub fn advance(&mut self, dt: T) -> Result<()> {
        std::mem::swap(&mut self.phi, &mut self.phi_previous);
        self.advect_weno5(dt)?;

        self.time_step += 1;
        if self
            .time_step
            .is_multiple_of(self.config.reinitialization_interval)
        {
            self.reinitialize()?;
        }

        if self.config.use_narrow_band {
            self.update_narrow_band();
        }

        Ok(())
    }

    // ────────────────────────────────────────────────────────────────────────────
    // WENO5-Z Advection
    // ────────────────────────────────────────────────────────────────────────────

    /// Copy a halo from the previous field into the current field.
    ///
    /// The WENO advection step requires a 3-cell halo, while the
    /// reinitialization step only requires a 1-cell boundary layer. Copying
    /// just the requested shell preserves the current state without paying a
    /// full-domain memory copy on every time step.
    fn copy_halo_from_previous(&mut self, halo: usize) {
        if halo == 0 || self.nx == 0 || self.ny == 0 || self.nz == 0 {
            return;
        }

        let halo_x = halo.min(self.nx);
        let halo_y = halo.min(self.ny);
        let halo_z = halo.min(self.nz);
        let x_max = self.nx.saturating_sub(halo_x);
        let y_max = self.ny.saturating_sub(halo_y);
        let z_max = self.nz.saturating_sub(halo_z);

        for k in 0..self.nz {
            let k_shell = k < halo_z || k >= z_max;
            for j in 0..self.ny {
                let j_shell = j < halo_y || j >= y_max;
                for i in 0..self.nx {
                    if k_shell || j_shell || i < halo_x || i >= x_max {
                        let idx = self.index(i, j, k);
                        self.phi[idx] = self.phi_previous[idx];
                    }
                }
            }
        }
    }

    /// Advect the level set using the 5th-order WENO-Z scheme.
    ///
    /// For each direction we compute the one-sided WENO5 reconstruction of the
    /// spatial derivative:
    ///
    /// ```text
    /// φ_x⁻  (left-biased, used when u > 0)
    /// φ_x⁺  (right-biased, used when u < 0)
    /// ```
    ///
    /// Then the Godunov upwind Hamiltonian is:
    ///
    /// ```math
    /// H(φ_x⁻, φ_x⁺, u) = u φ_x   where φ_x = φ_x⁻ if u > 0, φ_x⁺ if u < 0
    /// ```
    ///
    /// A single Forward-Euler time step is used here.  For production use couple
    /// with TVD-RK3 (Shu & Osher 1988).
    fn advect_weno5(&mut self, dt: T) -> Result<()> {
        let nx = self.nx;
        let ny = self.ny;
        let nz = self.nz;

        // WENO5 flux-differencing requires 3 ghost cells per side.
        for k in 3..nz.saturating_sub(3) {
            for j in 3..ny.saturating_sub(3) {
                for i in 3..nx.saturating_sub(3) {
                    let idx = self.index(i, j, k);
                    let vel = self.velocity[idx];

                    // ── x-direction: 7-point stencil {φ_{i-3},...,φ_{i+3}} ──
                    let vx = [
                        self.phi_previous[self.index(i - 3, j, k)],
                        self.phi_previous[self.index(i - 2, j, k)],
                        self.phi_previous[self.index(i - 1, j, k)],
                        self.phi_previous[idx],
                        self.phi_previous[self.index(i + 1, j, k)],
                        self.phi_previous[self.index(i + 2, j, k)],
                        self.phi_previous[self.index(i + 3, j, k)],
                    ];
                    let dphi_dx = weno5_derivative(vx, self.dx, vel.x);

                    // ── y-direction ──
                    let vy = [
                        self.phi_previous[self.index(i, j - 3, k)],
                        self.phi_previous[self.index(i, j - 2, k)],
                        self.phi_previous[self.index(i, j - 1, k)],
                        self.phi_previous[idx],
                        self.phi_previous[self.index(i, j + 1, k)],
                        self.phi_previous[self.index(i, j + 2, k)],
                        self.phi_previous[self.index(i, j + 3, k)],
                    ];
                    let dphi_dy = weno5_derivative(vy, self.dy, vel.y);

                    // ── z-direction ──
                    let vz = [
                        self.phi_previous[self.index(i, j, k - 3)],
                        self.phi_previous[self.index(i, j, k - 2)],
                        self.phi_previous[self.index(i, j, k - 1)],
                        self.phi_previous[idx],
                        self.phi_previous[self.index(i, j, k + 1)],
                        self.phi_previous[self.index(i, j, k + 2)],
                        self.phi_previous[self.index(i, j, k + 3)],
                    ];
                    let dphi_dz = weno5_derivative(vz, self.dz, vel.z);

                    self.phi[idx] = self.phi_previous[idx]
                        - dt * (vel.x * dphi_dx + vel.y * dphi_dy + vel.z * dphi_dz);
                }
            }
        }

        self.copy_halo_from_previous(3);

        Ok(())
    }

    // ────────────────────────────────────────────────────────────────────────────
    // Reinitialization
    // ────────────────────────────────────────────────────────────────────────────

    /// Reinitialize the level set to a signed distance function.
    ///
    /// Solves the reinitialization PDE (Sussman, Smereka & Osher 1994):
    ///
    /// ```math
    /// ∂φ/∂τ + S(φ₀)(|∇φ| − 1) = 0
    /// ```
    ///
    /// using the Godunov upwind Hamiltonian with first-order one-sided differences.
    ///
    /// # Theorem (Sussman reinitialization)
    ///
    /// If the initial zero level set is not modified and the PDE is marched
    /// in pseudo-time τ to steady state, then the solution converges to the
    /// signed distance function: |∇φ| = 1 everywhere, with the same zero
    /// level set as φ₀.
    ///
    /// **Proof sketch**: The smoothed sign function
    /// `S(φ₀) = φ₀ / √(φ₀² + Δx²)` vanishes at the interface, pinning it
    /// in place. Away from the interface, the Godunov Hamiltonian propagates
    /// information outward (S > 0) or inward (S < 0) until |∇φ| = 1.
    fn reinitialize(&mut self) -> Result<()> {
        use num_traits::Float;

        let dx_min = Float::min(Float::min(self.dx, self.dy), self.dz);
        // CFL-stable pseudo-time step for 3D: dtau ≤ dx/(2*dim) to stay below CFL=1.
        let dtau = dx_min
            / <T as FromPrimitive>::from_f64(6.0)
                .expect("6.0 is representable in all IEEE 754 types");
        let tol = <T as FromPrimitive>::from_f64(1e-6)
            .expect("1e-6 is an IEEE 754 representable f64 constant");

        // Store φ₀ for the sign function (must not be overwritten during iteration).
        self.phi_reinit.copy_from_slice(&self.phi);

        let max_iters = 100usize;

        for _ in 0..max_iters {
            std::mem::swap(&mut self.phi, &mut self.phi_previous);

            let mut max_err = T::zero();

            for k in 1..self.nz.saturating_sub(1) {
                for j in 1..self.ny.saturating_sub(1) {
                    for i in 1..self.nx.saturating_sub(1) {
                        let idx = self.index(i, j, k);

                        // S(φ₀) = φ₀ / √(φ₀² + Δx²)
                        let phi0_val = self.phi_reinit[idx];
                        let sign_phi =
                            phi0_val / Float::sqrt(phi0_val * phi0_val + dx_min * dx_min);

                        let grad_mag = self.godunov_gradient_magnitude(i, j, k, sign_phi);

                        self.phi[idx] =
                            self.phi_previous[idx] - dtau * sign_phi * (grad_mag - T::one());

                        let err = Float::abs(grad_mag - T::one());
                        if err > max_err {
                            max_err = err;
                        }
                    }
                }
            }

            self.copy_halo_from_previous(1);

            if max_err < tol {
                break;
            }
        }

        Ok(())
    }

    /// First-order Godunov upwind gradient magnitude for the reinitialization PDE.
    ///
    /// Uses one-sided first-order differences (Sussman et al. 1994):
    ///
    /// ```text
    /// D⁻x φ = (φ_i − φ_{i−1}) / Δx   (backward difference)
    /// D⁺x φ = (φ_{i+1} − φ_i) / Δx   (forward difference)
    /// ```
    ///
    /// Godunov Hamiltonian selects:
    /// - If S > 0: |∇φ|² = Σ_d [ max(max(D⁻,0)², max(−D⁺,0)²) ]
    /// - If S < 0: |∇φ|² = Σ_d [ max(max(D⁺,0)², max(−D⁻,0)²) ]
    fn godunov_gradient_magnitude(&self, i: usize, j: usize, k: usize, sign_phi: T) -> T {
        use num_traits::Float;

        let zero = T::zero();
        let idx = self.index(i, j, k);

        // First-order one-sided differences from phi_previous
        let dm_x = (self.phi_previous[idx]
            - self.phi_previous[self.index(i.saturating_sub(1), j, k)])
            / self.dx;
        let dp_x = (self.phi_previous[self.index((i + 1).min(self.nx - 1), j, k)]
            - self.phi_previous[idx])
            / self.dx;

        let dm_y = (self.phi_previous[idx]
            - self.phi_previous[self.index(i, j.saturating_sub(1), k)])
            / self.dy;
        let dp_y = (self.phi_previous[self.index(i, (j + 1).min(self.ny - 1), k)]
            - self.phi_previous[idx])
            / self.dy;

        let dm_z = (self.phi_previous[idx]
            - self.phi_previous[self.index(i, j, k.saturating_sub(1))])
            / self.dz;
        let dp_z = (self.phi_previous[self.index(i, j, (k + 1).min(self.nz - 1))]
            - self.phi_previous[idx])
            / self.dz;

        let mut grad_sq = T::zero();

        if sign_phi > zero {
            // Propagate outward: use max(D⁻,0)² and max(-D⁺,0)²
            let ax = Float::max(dm_x, zero);
            let bx = Float::min(dp_x, zero);
            grad_sq += Float::max(ax * ax, bx * bx);

            let ay = Float::max(dm_y, zero);
            let by = Float::min(dp_y, zero);
            grad_sq += Float::max(ay * ay, by * by);

            let az = Float::max(dm_z, zero);
            let bz = Float::min(dp_z, zero);
            grad_sq += Float::max(az * az, bz * bz);
        } else {
            // Propagate inward: use max(D⁺,0)² and max(-D⁻,0)²
            let ax = Float::min(dm_x, zero);
            let bx = Float::max(dp_x, zero);
            grad_sq += Float::max(ax * ax, bx * bx);

            let ay = Float::min(dm_y, zero);
            let by = Float::max(dp_y, zero);
            grad_sq += Float::max(ay * ay, by * by);

            let az = Float::min(dm_z, zero);
            let bz = Float::max(dp_z, zero);
            grad_sq += Float::max(az * az, bz * bz);
        }

        Float::sqrt(grad_sq)
    }
}
