//! # Immersed Boundary Method (IBM) Solver for 3D Complex Geometries
//!
//! This module implements the Immersed Boundary Method for simulating fluid flow
//! around complex geometries without body-fitted meshes.
//!
//! ## Mathematical Foundation
//!
//! ### Direct Forcing Theorem
//!
//! **Statement**: The direct forcing method enforces no-slip boundary conditions by
//! adding momentum source terms to the Navier-Stokes equations at immersed boundary points.
//!
//! **Mathematical Formulation**:
//!
//! ```math
//! ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + f
//! ```
//!
//! where the forcing term f is defined as:
//!
//! ```math
//! f(x,t) = ∫_Γ F(s,t) δ(x - X(s,t)) ds
//! ```
//!
//! **Discrete Implementation**: The forcing is applied as:
//!
//! ```math
//! f_i = (U_desired - u_i) / Δt    for points near boundary
//! ```
//!
//! **Convergence Properties**:
//! - **First-order accuracy** in boundary condition enforcement
//! - **Stability**: CFL condition must be satisfied for explicit forcing
//! - **Conservation**: Momentum conservation maintained through force balance
//!
//! **Assumptions**:
//! 1. **Linear forcing response** within each time step
//! 2. **Localized forcing** effects decay away from boundary
//! 3. **Small boundary displacement** compared to grid spacing
//!
//! **Literature**: Fadlun, E.A., Verzicco, R., Orlandi, P., Mohd-Yusof, J. (2000).
//! "Combined immersed-boundary finite-difference methods for three-dimensional
//! complex flow simulations". Journal of Computational Physics, 161(1), 35-60.
//!
//! ### Feedback Forcing Theorem
//!
//! **Statement**: The feedback forcing method uses a proportional-integral-derivative
//! (PID) controller to enforce boundary conditions through corrective forcing.
//!
//! **Mathematical Formulation**: The forcing term is computed as:
//!
//! ```math
//! f^{n+1} = K_p (U_desired - u^n) + K_i ∫(U_desired - u) dt + K_d d(U_desired - u)/dt
//! ```
//!
//! **PID Controller Parameters**:
//! - **Proportional Gain (K_p)**: Immediate correction based on current error
//! - **Integral Gain (K_i)**: Correction based on accumulated error over time
//! - **Derivative Gain (K_d)**: Anticipatory correction based on error trend
//!
//! **Stability Analysis**: The controller is stable if:
//!
//! ```math
//! |K_p| < 2/Δt,    K_i > 0,    K_d < Δt/2
//! ```
//!
//! **Convergence Properties**:
//! - **Exponential convergence** to desired boundary conditions
//! - **Robustness**: Less sensitive to grid resolution than direct forcing
//! - **Accuracy**: Higher-order boundary condition enforcement possible
//!
//! **Literature**: Goldstein, D., Handler, R., Sirovich, L. (1993).
//! "Modeling a no-slip flow boundary with an external force field".
//! Journal of Computational Physics, 105(2), 354-366.
//!
//! ### Interpolation Theory (Discrete Delta Functions)
//!
//! **Statement**: The discrete delta function provides a smooth interpolation kernel
//! for transferring quantities between Eulerian grid and Lagrangian boundary points.
//!
//! **Mathematical Formulation**: The discrete delta function D_h(x) satisfies:
//!
//! ```math
//! ∫_Ω D_h(x) dx = 1,    D_h(x) ≥ 0,    supp(D_h) ⊂ [-2h, 2h]^d
//! ```
//!
//! **Roma-Peskin Kernel**: The 4-point Roma kernel is defined as:
//!
//! ```math
//! D_h(r) = (1/h^d) ∏_{i=1}^d φ(r_i/h)
//! ```
//!
//! where φ is the 1D kernel:
//!
//! ```math
//! φ(r) = (1/8)(3 - 2|r| + √(1 + 4|r| - 4r²))    for |r| ≤ 1
//! φ(r) = (1/8)(5 - 2|r| - √(-7 + 12|r| - 4r²))  for 1 < |r| ≤ 2
//! φ(r) = 0                                       for |r| > 2
//! ```
//!
//! **Interpolation Properties**:
//! - **Conservation**: ∫ f δ_h(x - X) dx = f(X) (interpolation)
//! - **Adjoint Property**: ∫ f δ_h(x - X) dx = f(X) (spreading)
//! - **Accuracy**: Second-order accurate for smooth functions
//!
//! **Literature**: Roma, A.M., Peskin, C.S., Berger, M.J. (1999).
//! "An adaptive version of the immersed boundary method".
//! Journal of Computational Physics, 153(2), 509-534.
//!
//! Peskin, C.S. (2002). "The immersed boundary method".
//! Acta Numerica, 11, 479-517.
//!
//! ### Algorithm Implementation
//!
//! 1. **Boundary Representation**: Lagrangian points define immersed geometry
//! 2. **Force Computation**: Boundary forces computed from desired conditions
//! 3. **Force Spreading**: Forces spread to Eulerian grid using δ_h function
//! 4. **Flow Solution**: Navier-Stokes equations solved with forcing terms
//! 5. **Velocity Interpolation**: Velocities interpolated back to boundary points
//! 6. **Condition Enforcement**: Boundary conditions enforced via forcing update
//!
//! ### Convergence and Stability
//!
//! **Theorem (IBM Convergence)**: For sufficiently resolved boundaries (h << geometric features),
//! the IBM converges to the body-fitted solution with order O(h²) accuracy.
//!
//! **Stability Condition**: The forcing time scale must satisfy τ_force << Δt_flow
//! to maintain numerical stability.
//!
//! **Literature**: Mittal, R., Iaccarino, G. (2005). "Immersed boundary methods".
//! Annual Review of Fluid Mechanics, 37, 239-261.

use super::{
    config::IbmConfig,
    forcing::{DirectForcing, FeedbackForcing, ForcingMethod},
    interpolation::{DeltaFunction, InterpolationKernel},
    lagrangian::LagrangianPoint,
};
use cfd_core::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::{FromPrimitive, ToPrimitive};

// Feedback control constants
const DEFAULT_PROPORTIONAL_GAIN: f64 = 10.0;
const DEFAULT_INTEGRAL_GAIN: f64 = 1.0;

/// IBM solver for 3D flow around immersed boundaries
pub struct IbmSolver<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + ToPrimitive + Copy> {
    /// Configuration
    #[allow(dead_code)]
    config: IbmConfig,
    /// Lagrangian points representing the immersed boundary
    lagrangian_points: Vec<LagrangianPoint<T>>,
    /// Interpolation kernel
    kernel: InterpolationKernel<T>,
    /// Forcing method
    forcing: Box<dyn ForcingMethod<T>>,
    /// Grid spacing
    dx: Vector3<T>,
    /// Grid dimensions
    grid_size: (usize, usize, usize),
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + ToPrimitive + Copy> IbmSolver<T> {
    /// Create a new IBM solver
    pub fn new(config: IbmConfig, dx: Vector3<T>, grid_size: (usize, usize, usize)) -> Self {
        let kernel = InterpolationKernel::new(
            DeltaFunction::RomaPeskin4,
            <T as FromPrimitive>::from_f64(config.smoothing_width).unwrap_or_else(T::one),
        );

        let forcing: Box<dyn ForcingMethod<T>> = if config.use_direct_forcing {
            Box::new(DirectForcing::new(
                <T as FromPrimitive>::from_f64(config.force_scale).unwrap_or_else(T::one),
            ))
        } else {
            Box::new(FeedbackForcing::new(
                <T as FromPrimitive>::from_f64(DEFAULT_PROPORTIONAL_GAIN).unwrap_or_else(T::one),
                <T as FromPrimitive>::from_f64(DEFAULT_INTEGRAL_GAIN).unwrap_or_else(T::one),
            ))
        };

        Self {
            config,
            lagrangian_points: Vec::new(),
            kernel,
            forcing,
            dx,
            grid_size,
        }
    }

    /// Add a Lagrangian point
    pub fn add_lagrangian_point(&mut self, point: LagrangianPoint<T>) {
        self.lagrangian_points.push(point);
    }

    /// Interpolate velocity from Eulerian to Lagrangian grid
    pub fn interpolate_velocity(&mut self, eulerian_velocity: &[Vector3<T>]) -> Result<()> {
        // Collect velocities first to avoid borrow conflict
        let velocities: Vec<_> = self
            .lagrangian_points
            .iter()
            .map(|point| self.interpolate_at_point(&point.position, eulerian_velocity))
            .collect::<Result<Vec<_>>>()?;

        // Then update the points
        for (point, velocity) in self.lagrangian_points.iter_mut().zip(velocities) {
            point.velocity = velocity;
        }
        Ok(())
    }

    /// Calculate forcing terms
    pub fn calculate_forcing(&mut self, desired_velocity: &[Vector3<T>], dt: T) -> Result<()> {
        for (i, point) in self.lagrangian_points.iter_mut().enumerate() {
            point.force = self
                .forcing
                .calculate_force(&desired_velocity[i], &point.velocity, dt);
        }
        Ok(())
    }

    /// Spread forces from Lagrangian to Eulerian grid
    pub fn spread_forces(&self) -> Result<Vec<Vector3<T>>> {
        let total_points = self.grid_size.0 * self.grid_size.1 * self.grid_size.2;
        let mut eulerian_forces = vec![Vector3::zeros(); total_points];

        for point in &self.lagrangian_points {
            self.spread_from_point(point, &mut eulerian_forces)?;
        }

        Ok(eulerian_forces)
    }

    /// Interpolate at a specific point
    fn interpolate_at_point(
        &self,
        position: &Vector3<T>,
        field: &[Vector3<T>],
    ) -> Result<Vector3<T>> {
        let mut result = Vector3::zeros();
        let stencil = self.kernel.stencil_size();

        // Find grid indices and convert to integers
        let i_int = (num_traits::Float::floor(position.x / self.dx.x)).to_isize().unwrap_or(0);
        let j_int = (num_traits::Float::floor(position.y / self.dx.y)).to_isize().unwrap_or(0);
        let k_int = (num_traits::Float::floor(position.z / self.dx.z)).to_isize().unwrap_or(0);

        let i_start = i_int - (stencil as isize / 2);
        let j_start = j_int - (stencil as isize / 2);
        let k_start = k_int - (stencil as isize / 2);

        for di in 0..stencil {
            for dj in 0..stencil {
                for dk in 0..stencil {
                    let ii = (i_start + di as isize) as usize;
                    let jj = (j_start + dj as isize) as usize;
                    let kk = (k_start + dk as isize) as usize;

                    if ii < self.grid_size.0 && jj < self.grid_size.1 && kk < self.grid_size.2 {
                        let idx =
                            ii + jj * self.grid_size.0 + kk * self.grid_size.0 * self.grid_size.1;

                        let rx = position.x / self.dx.x
                            - T::from_usize(ii).unwrap_or_else(T::zero);
                        let ry = position.y / self.dx.y
                            - T::from_usize(jj).unwrap_or_else(T::zero);
                        let rz = position.z / self.dx.z
                            - T::from_usize(kk).unwrap_or_else(T::zero);

                        let weight =
                            self.kernel.delta(rx) * self.kernel.delta(ry) * self.kernel.delta(rz);
                        result += field[idx] * weight;
                    }
                }
            }
        }

        Ok(result)
    }

    /// Spread from a Lagrangian point to Eulerian grid
    fn spread_from_point(
        &self,
        point: &LagrangianPoint<T>,
        field: &mut [Vector3<T>],
    ) -> Result<()> {
        let stencil = self.kernel.stencil_size();

        // Find grid indices and convert to integers
        let i_int = (num_traits::Float::floor(point.position.x / self.dx.x))
            .to_isize()
            .unwrap_or(0);
        let j_int = (num_traits::Float::floor(point.position.y / self.dx.y))
            .to_isize()
            .unwrap_or(0);
        let k_int = (num_traits::Float::floor(point.position.z / self.dx.z))
            .to_isize()
            .unwrap_or(0);

        let i_start = i_int - (stencil as isize / 2);
        let j_start = j_int - (stencil as isize / 2);
        let k_start = k_int - (stencil as isize / 2);

        for di in 0..stencil {
            for dj in 0..stencil {
                for dk in 0..stencil {
                    let ii = (i_start + di as isize) as usize;
                    let jj = (j_start + dj as isize) as usize;
                    let kk = (k_start + dk as isize) as usize;

                    if ii < self.grid_size.0 && jj < self.grid_size.1 && kk < self.grid_size.2 {
                        let idx =
                            ii + jj * self.grid_size.0 + kk * self.grid_size.0 * self.grid_size.1;

                        let rx = point.position.x / self.dx.x
                            - T::from_usize(ii).unwrap_or_else(T::zero);
                        let ry = point.position.y / self.dx.y
                            - T::from_usize(jj).unwrap_or_else(T::zero);
                        let rz = point.position.z / self.dx.z
                            - T::from_usize(kk).unwrap_or_else(T::zero);

                        let weight =
                            self.kernel.delta(rx) * self.kernel.delta(ry) * self.kernel.delta(rz);
                        field[idx] += point.force * weight * point.weight;
                    }
                }
            }
        }

        Ok(())
    }

    /// Update Lagrangian point positions
    pub fn update_positions(&mut self, dt: T) {
        for point in &mut self.lagrangian_points {
            point.update_position(dt);
        }
    }

    /// Get number of Lagrangian points
    pub fn num_points(&self) -> usize {
        self.lagrangian_points.len()
    }

    /// Get reference to Lagrangian points
    pub fn points(&self) -> &[LagrangianPoint<T>] {
        &self.lagrangian_points
    }

    /// Get mutable reference to Lagrangian points
    pub fn points_mut(&mut self) -> &mut [LagrangianPoint<T>] {
        &mut self.lagrangian_points
    }
}
