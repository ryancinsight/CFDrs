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
//! ### WENO5-JS Stability Theorem (Jiang & Shu 1996)
//!
//! **Statement**: The 5th-order WENO scheme of Jiang & Shu (1996) applied to
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
//! **Reference**: Jiang, G.-S. & Shu, C.-W. (1996). "Efficient implementation of
//! Weighted ENO schemes." J. Comput. Phys. 126:202-228.
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
//! 2. **Evolution**: Solve ∂φ/∂t + u·∇φ = 0 using WENO5-JS (one Euler step per advance).
//! 3. **Reinitialization**: Restore |∇φ| = 1 every `config.reinitialization_interval` steps
//!    using convergence-based iteration.
//! 4. **Narrow Band**: Update index set of cells within `band_width` grid spacings
//!    of the interface.

use super::config::LevelSetConfig;
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
        let band_width =
            <T as FromPrimitive>::from_f64(self.config.band_width)
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

    /// Advance level set by one time step using WENO5-JS advection.
    pub fn advance(&mut self, dt: T) -> Result<()> {
        self.phi_previous.copy_from_slice(&self.phi);
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
    // WENO5-JS Advection
    // ────────────────────────────────────────────────────────────────────────────

    /// Advect the level set using the 5th-order WENO-JS scheme (Jiang & Shu 1996).
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

        std::mem::swap(&mut self.phi, &mut self.phi_previous);

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
        let dtau = dx_min / <T as FromPrimitive>::from_f64(6.0)
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

// ────────────────────────────────────────────────────────────────────────────
// WENO5-JS Free Functions
// ────────────────────────────────────────────────────────────────────────────

/// Compute the upwind WENO5-JS spatial derivative `dφ/dx` at point `i`
/// via flux differencing (Jiang & Shu 1996, §2).
///
/// Takes a 7-point stencil `v = [φ_{i-3}, …, φ_{i+3}]`, cell spacing `h`,
/// and the local velocity component `u`.
///
/// # Algorithm
///
/// The derivative is assembled from left-biased (`D⁻`) and right-biased (`D⁺`)
/// flux-difference reconstructions:
///
/// ```text
/// D⁻φ_i = (φ̂⁻_{i+½} − φ̂⁻_{i−½}) / h   (u > 0: upwind from left)
/// D⁺φ_i = (φ̂⁺_{i+½} − φ̂⁺_{i−½}) / h   (u < 0: upwind from right)
/// ```
///
/// where each `φ̂` is a WENO5-JS reconstruction at a cell face from 5 point values.
fn weno5_derivative<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    v: [T; 7],
    h: T,
    u: T,
) -> T {
    let half = <T as FromPrimitive>::from_f64(0.5)
        .expect("0.5 is exactly representable in IEEE 754");

    // Left-biased derivative D⁻φ_i = (φ̂⁻_{i+½} - φ̂⁻_{i-½}) / h
    // φ̂⁻_{i+½} from {v[1]..v[5]} = {φ_{i-2},..,φ_{i+2}}
    // φ̂⁻_{i-½} from {v[0]..v[4]} = {φ_{i-3},..,φ_{i+1}}
    let fl_right = weno5_reconstruct_left([v[1], v[2], v[3], v[4], v[5]]);
    let fl_left = weno5_reconstruct_left([v[0], v[1], v[2], v[3], v[4]]);
    let dm = (fl_right - fl_left) / h;

    // Right-biased derivative D⁺φ_i = (φ̂⁺_{i+½} - φ̂⁺_{i-½}) / h
    // φ̂⁺_{i+½} from {v[2]..v[6]} = {φ_{i-1},..,φ_{i+3}} (mirrored)
    // φ̂⁺_{i-½} from {v[1]..v[5]} = {φ_{i-2},..,φ_{i+2}} (mirrored)
    let fr_right = weno5_reconstruct_right([v[2], v[3], v[4], v[5], v[6]]);
    let fr_left = weno5_reconstruct_right([v[1], v[2], v[3], v[4], v[5]]);
    let dp = (fr_right - fr_left) / h;

    if u > T::zero() {
        dm
    } else if u < T::zero() {
        dp
    } else {
        (dm + dp) * half
    }
}

/// Left-biased WENO5-JS reconstruction of `φ` at the right face of the
/// middle cell (Jiang & Shu 1996, Eq. 2.10).
///
/// Given `[φ_{i-2}, φ_{i-1}, φ_i, φ_{i+1}, φ_{i+2}]`, reconstructs
/// `φ̂⁻_{i+½}` using three overlapping 3rd-order sub-stencils:
///
/// ```text
/// q₀ = (2 φ_{i-2} − 7 φ_{i-1} + 11 φ_i) / 6         stencil S₀ = {i-2, i-1, i}
/// q₁ = (−  φ_{i-1} + 5 φ_i     +  2 φ_{i+1}) / 6     stencil S₁ = {i-1, i, i+1}
/// q₂ = (2 φ_i     + 5 φ_{i+1} −    φ_{i+2}) / 6      stencil S₂ = {i, i+1, i+2}
/// ```
///
/// Ideal weights for 5th-order: `d = (1/10, 6/10, 3/10)`.
fn weno5_reconstruct_left<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    v: [T; 5],
) -> T {
    let eps = <T as FromPrimitive>::from_f64(1e-36)
        .expect("1e-36 is an IEEE 754 representable f64 constant");
    let one = T::one();
    let two = one + one;
    let three = two + one;
    let five = three + two;
    let six = three + three;
    let seven = three + three + one;
    let eleven = five + six;

    // Candidate reconstructions
    let q0 = (two * v[0] - seven * v[1] + eleven * v[2]) / six;
    let q1 = (-v[1] + five * v[2] + two * v[3]) / six;
    let q2 = (two * v[2] + five * v[3] - v[4]) / six;

    // Smoothness indicators (Jiang & Shu 1996, Eq. 3.1)
    let b0 = smoothness_indicator(v[0], v[1], v[2]);
    let b1 = smoothness_indicator(v[1], v[2], v[3]);
    let b2 = smoothness_indicator(v[2], v[3], v[4]);

    // Nonlinear weights
    let (w0, w1, w2) = nonlinear_weights(
        b0,
        b1,
        b2,
        <T as FromPrimitive>::from_f64(0.1)
            .expect("0.1 is an IEEE 754 representable f64 constant"),
        <T as FromPrimitive>::from_f64(0.6)
            .expect("0.6 is an IEEE 754 representable f64 constant"),
        <T as FromPrimitive>::from_f64(0.3)
            .expect("0.3 is an IEEE 754 representable f64 constant"),
        eps,
    );

    w0 * q0 + w1 * q1 + w2 * q2
}

/// Right-biased WENO5-JS reconstruction of `φ` at the left face of the
/// middle cell, obtained by mirroring the left-biased stencil.
///
/// Given `[φ_{i-2}, φ_{i-1}, φ_i, φ_{i+1}, φ_{i+2}]`, reconstructs
/// `φ̂⁺_{i−½}` using:
///
/// ```text
/// q₀ = (2 φ_{i+2} − 7 φ_{i+1} + 11 φ_i) / 6         stencil S₀ = {i+2, i+1, i}
/// q₁ = (−  φ_{i+1} + 5 φ_i     +  2 φ_{i-1}) / 6     stencil S₁ = {i+1, i, i-1}
/// q₂ = (2 φ_i     + 5 φ_{i-1} −    φ_{i-2}) / 6      stencil S₂ = {i, i-1, i-2}
/// ```
///
/// Ideal weights: `d = (1/10, 6/10, 3/10)`.
fn weno5_reconstruct_right<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    v: [T; 5],
) -> T {
    // Mirror: weno5_left([v[4], v[3], v[2], v[1], v[0]])
    weno5_reconstruct_left([v[4], v[3], v[2], v[1], v[0]])
}

/// WENO5-JS smoothness indicator for a 3-point sub-stencil (Jiang & Shu 1996,
/// Eq. 3.1).
///
/// ```text
/// β = (13/12)(v₀ − 2v₁ + v₂)² + (1/4)(v₀ − v₂)²
/// ```
///
/// The first term measures the curvature (second derivative), the second
/// measures the slope difference.  Near discontinuities β is O(1), in smooth
/// regions β is O(h⁴), which drives the nonlinear weights toward the ideal values.
#[inline]
fn smoothness_indicator<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    v0: T,
    v1: T,
    v2: T,
) -> T {
    let thirteen_over_twelve = <T as FromPrimitive>::from_f64(13.0 / 12.0)
        .expect("13/12 is an IEEE 754 representable f64 constant");
    let quarter = <T as FromPrimitive>::from_f64(0.25)
        .expect("0.25 is exactly representable in IEEE 754");
    let two = T::one() + T::one();

    let diff1 = v0 - two * v1 + v2; // ≈ h² φ''
    let diff2 = v0 - v2; // ≈ 2h φ'
    thirteen_over_twelve * diff1 * diff1 + quarter * diff2 * diff2
}

/// Compute normalized WENO5 nonlinear weights from smoothness indicators.
///
/// `ω_k = α_k / Σ α_l` where `α_k = d_k / (ε + β_k)²`.
#[inline]
fn nonlinear_weights<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + Copy>(
    b0: T,
    b1: T,
    b2: T,
    d0: T,
    d1: T,
    d2: T,
    eps: T,
) -> (T, T, T) {
    let a0 = d0 / ((eps + b0) * (eps + b0));
    let a1 = d1 / ((eps + b1) * (eps + b1));
    let a2 = d2 / ((eps + b2) * (eps + b2));
    let sum = a0 + a1 + a2;
    if sum < eps {
        return (
            d0 / (d0 + d1 + d2),
            d1 / (d0 + d1 + d2),
            d2 / (d0 + d1 + d2),
        );
    }
    (a0 / sum, a1 / sum, a2 / sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use proptest::prelude::*;

    proptest! {
        #[test]
        fn test_smoothness_indicator_non_negative(
            v0 in -10.0..10.0f64,
            v1 in -10.0..10.0f64,
            v2 in -10.0..10.0f64,
        ) {
            let beta = smoothness_indicator(v0, v1, v2);
            assert!(beta >= 0.0);
        }

        #[test]
        fn test_nonlinear_weights_sum_to_one(
            b0 in 0.0..10.0f64,
            b1 in 0.0..10.0f64,
            b2 in 0.0..10.0f64,
        ) {
            let d0 = 0.1;
            let d1 = 0.6;
            let d2 = 0.3;
            let eps = 1e-36;

            let (w0, w1, w2) = nonlinear_weights(b0, b1, b2, d0, d1, d2, eps);
            let sum = w0 + w1 + w2;

            // Sum of weights must be 1.0 within machine precision
            assert!((sum - 1.0).abs() < 1e-14);

            // Weights must be non-negative
            assert!(w0 >= 0.0);
            assert!(w1 >= 0.0);
            assert!(w2 >= 0.0);
        }
    }
}
