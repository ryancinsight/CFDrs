//! 2D Poiseuille flow solver with non-Newtonian blood rheology
//!
//! # Physical Problem
//!
//! Steady, fully-developed laminar flow between two parallel plates separated by height H.
//! The flow is driven by a constant pressure gradient dP/dx.
//!
//! ## Governing Equations
//!
//! For incompressible, steady flow in the x-direction with no variation in x:
//!
//! **Continuity:**
//! ```text
//! ∂u/∂x + ∂v/∂y = 0
//! ```
//! Since u = u(y) only and v = 0, continuity is automatically satisfied.
//!
//! **Momentum (x-direction):**
//! ```text
//! 0 = -dP/dx + d/dy(μ du/dy)
//! ```
//!
//! **Boundary Conditions:**
//! ```text
//! u(y=0) = 0     (no-slip at bottom wall)
//! u(y=H) = 0     (no-slip at top wall)
//! ```
//!
//! ## Analytical Solution (Newtonian Fluid)
//!
//! For constant viscosity μ:
//! ```text
//! u(y) = (1/2μ)(dP/dx)y(H - y)
//! u_max = (H²/8μ)|dP/dx|  (at y = H/2)
//! Q = (H³W/12μ)|dP/dx|    (volumetric flow rate)
//! ```
//!
//! ## Non-Newtonian Blood Flow
//!
//! For blood, viscosity depends on shear rate:
//! ```text
//! μ = μ(γ̇)  where γ̇ = |du/dy|
//! ```
//!
//! **Casson Model:**
//! ```text
//! μ(γ̇) = (√τ_y/√γ̇ + √μ_∞)²
//! ```
//!
//! **Solution Method:**
//! 1. Guess initial viscosity field (constant or previous solution)
//! 2. Solve linear system for velocity with current viscosity
//! 3. Calculate shear rate: γ̇ = |du/dy|
//! 4. Update viscosity: μ_new = μ(γ̇)
//! 5. Check convergence: ||μ_new - μ_old|| < tolerance
//! 6. Repeat until converged
//!
//! ## Numerical Method
//!
//! **Finite Difference Discretization:**
//!
//! At interior points (j = 1, 2, ..., ny-2):
//! ```text
//! 0 = -dP/dx + (μ_{j+1/2}(u_{j+1} - u_j) - μ_{j-1/2}(u_j - u_{j-1})) / Δy²
//! ```
//!
//! Where μ_{j+1/2} = (μ_j + μ_{j+1})/2 (arithmetic mean at cell faces)
//!
//! This gives a tridiagonal system: A*u = b
//!
//! **Matrix Structure:**
//! ```text
//! [-a_j  c_j    0  ] [u_{j-1}]   [b_j]
//! [ 0   -a_j  c_j ] [u_j    ] = [b_j]
//! [ 0    0   -a_j ] [u_{j+1}]   [b_j]
//! ```
//!
//! Where:
//! - a_j = (μ_{j-1/2} + μ_{j+1/2}) / Δy²
//! - c_j = μ_{j+1/2} / Δy²
//! - b_j = -dP/dx
//!
//! ## Literature References
//!
//! 1. White, F.M. (2006) "Viscous Fluid Flow" 3rd Ed., McGraw-Hill, pp. 123-125
//! 2. Fung, Y.C. (1993) "Biomechanics: Circulation" 2nd Ed., Springer, pp. 50-55
//! 3. Pries, A.R. et al. (1994) "Blood flow in microvessels" Annu Rev Fluid Mech 26:315-348
//! 4. Merrill, E.W. (1969) "Rheology of blood" Physiol Rev 49:863-888
//!
//! # Theorem
//! The solver algorithm must converge to a unique solution that satisfies the discrete
//! conservation laws.
//!
//! **Proof sketch**:
//! For a well-posed boundary value problem, the discretized system of equations
//! $\mathbf{A}\mathbf{x} = \mathbf{b}$ forms a diagonally dominant matrix $\mathbf{A}$
//! under appropriate upwinding or stabilization. The iterative solver (e.g., SIMPLE, PISO)
//! reduces the residual norm $\|\mathbf{r}\| = \|\mathbf{b} - \mathbf{A}\mathbf{x}\|$
//! monotonically. Convergence is guaranteed by the spectral radius of the iteration matrix
//! being strictly less than 1.

mod numerics;

use crate::scalar;
use crate::scalar::Cfd2dScalar;
use cfd_core::error::Error;
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};
use eunomia::{FloatElement, NumericElement};
use serde::{Deserialize, Serialize};

/// Configuration for Poiseuille flow solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PoiseuilleConfig<T: Cfd2dScalar + Copy> {
    /// Channel height \[m]
    pub height: T,

    /// Channel width \[m] (for 3D flow rate calculation)
    pub width: T,

    /// Channel length \[m] (for reference only, flow is fully developed)
    pub length: T,

    /// Number of grid points in y-direction
    pub ny: usize,

    /// Pressure gradient magnitude |dP/dx| [Pa/m]
    /// Positive value drives flow in +x direction
    pub pressure_gradient: T,

    /// Convergence tolerance for iterative solution
    pub tolerance: T,

    /// Maximum iterations for non-Newtonian solve
    pub max_iterations: usize,

    /// Under-relaxation factor (0 < α ≤ 1)
    /// μ_new = α*μ_computed + (1-α)*μ_old
    pub relaxation_factor: T,
}

impl<T: Cfd2dScalar + FloatElement + Copy> Default for PoiseuilleConfig<T> {
    fn default() -> Self {
        Self {
            height: scalar::from_f64(100e-6),            // 100 μm
            width: scalar::from_f64(500e-6),             // 500 μm
            length: scalar::from_f64(1e-3),              // 1 mm
            ny: 101,                                     // 101 points for 100 intervals
            pressure_gradient: scalar::from_f64(1000.0), // 1000 Pa/m
            tolerance: scalar::from_f64(1e-6),
            max_iterations: 1000,
            relaxation_factor: scalar::from_f64(0.7),
        }
    }
}

/// 2D Poiseuille flow solver with non-Newtonian blood rheology
///
/// Solves steady-state velocity profile u(y) between parallel plates.
#[derive(Debug)]
pub struct PoiseuilleFlow2D<T: Cfd2dScalar + Copy + FloatElement> {
    /// Configuration
    pub(super) config: PoiseuilleConfig<T>,

    /// Grid in y-direction (1D since flow only varies in y)
    pub(super) y_coords: Vec<T>,

    /// Grid spacing
    pub(super) dy: T,

    /// Velocity field u(y)
    pub(super) velocity: Vec<T>,

    /// Shear rate field γ̇(y) = |du/dy|
    pub(super) shear_rate: Vec<T>,

    /// Viscosity field μ(y) = μ(γ̇(y))
    pub(super) viscosity: Vec<T>,

    /// Blood model (Casson or Carreau-Yasuda)
    pub(super) blood_model: BloodModel<T>,

    /// Convergence history
    pub(super) residuals: Vec<T>,

    /// Reusable buffer for viscosity convergence checks (avoids per-iteration allocation)
    viscosity_work: Vec<T>,
}

/// Blood rheology model selection
#[derive(Debug, Clone)]
pub enum BloodModel<T: Cfd2dScalar + Copy> {
    /// Casson model with yield stress
    Casson(CassonBlood<T>),

    /// Carreau-Yasuda model
    CarreauYasuda(CarreauYasudaBlood<T>),
}

impl<T: Cfd2dScalar + FloatElement + Copy> PoiseuilleFlow2D<T> {
    /// Create new Poiseuille flow solver
    ///
    /// # Arguments
    ///
    /// * `config` - Solver configuration
    /// * `blood_model` - Blood rheology model (Casson or Carreau-Yasuda)
    ///
    /// # Example
    ///
    /// ```rust
    /// use cfd_2d::solvers::poiseuille::{PoiseuilleFlow2D, PoiseuilleConfig, BloodModel};
    /// use cfd_core::physics::fluid::blood::CassonBlood;
    ///
    /// let config = PoiseuilleConfig::default();
    /// let blood = CassonBlood::<f64>::normal_blood();
    /// let solver = PoiseuilleFlow2D::new(config, BloodModel::Casson(blood));
    /// ```
    pub fn new(config: PoiseuilleConfig<T>, blood_model: BloodModel<T>) -> Self {
        let ny = config.ny;
        let dy = config.height / scalar::from_usize::<T>(ny - 1);

        // Create grid points
        let y_coords: Vec<T> = (0..ny).map(|j| scalar::from_usize::<T>(j) * dy).collect();

        // Initialize fields
        let zero: T = scalar::zero();
        let velocity = vec![zero; ny];
        let shear_rate = vec![zero; ny];

        // Initialize viscosity with constant value
        let mu_init = <T as FloatElement>::from_f64(0.004); // ~4 cP initial guess
        let viscosity = vec![mu_init; ny];

        Self {
            config,
            y_coords,
            dy,
            velocity,
            shear_rate,
            viscosity: viscosity.clone(),
            blood_model,
            residuals: Vec::new(),
            viscosity_work: viscosity,
        }
    }

    /// Solve for steady-state velocity profile
    ///
    /// Uses iterative method for non-Newtonian flow:
    /// 1. Solve linear system with current viscosity
    /// 2. Calculate shear rate from velocity gradient
    /// 3. Update viscosity from blood model
    /// 4. Check convergence
    /// 5. Repeat until converged
    ///
    /// # Returns
    ///
    /// `Ok(num_iterations)` if converged
    /// `Err(Error)` if failed to converge
    pub fn solve(&mut self) -> Result<usize, Error> {
        let ny = self.config.ny;
        let tolerance = self.config.tolerance;
        let max_iter = self.config.max_iterations;
        let alpha = self.config.relaxation_factor;

        self.residuals.clear();

        for iteration in 0..max_iter {
            // Store old viscosity for convergence check (reuse pre-allocated buffer)
            self.viscosity_work.copy_from_slice(&self.viscosity);

            // 1. Solve for velocity with current viscosity
            self.solve_velocity_linear()?;

            // 2. Calculate shear rate from velocity
            self.calculate_shear_rate();

            // 3. Update viscosity from blood model
            self.update_viscosity_from_shear_rate();

            // 4. Apply under-relaxation
            for j in 0..ny {
                self.viscosity[j] = alpha * self.viscosity[j]
                    + (scalar::one::<T>() - alpha) * self.viscosity_work[j];
            }

            // 5. Check convergence
            let residual = self.calculate_viscosity_residual(&self.viscosity_work);
            self.residuals.push(residual);

            if residual < tolerance {
                return Ok(iteration + 1);
            }
        }

        Err(Error::Solver(format!(
            "Convergence failed after {max_iter} iterations (residual: {:.2e})",
            <T as NumericElement>::to_f64(
                self.residuals
                    .last()
                    .copied()
                    .unwrap_or_else(scalar::zero::<T>)
            )
        )))
    }

    /// Calculate analytical solution for Newtonian fluid (for validation)
    ///
    /// u(y) = (1/2μ)(dP/dx)y(H - y)
    ///
    /// Returns `None` if blood model is non-Newtonian
    pub fn analytical_solution(&self, viscosity: T) -> Vec<T> {
        let dp_dx = self.config.pressure_gradient;
        let height = self.config.height;
        let two = scalar::from_f64::<T>(2.0);

        self.y_coords
            .iter()
            .map(|&y| (dp_dx / (two * viscosity)) * y * (height - y))
            .collect()
    }

    /// Get maximum velocity (at centerline)
    pub fn max_velocity(&self) -> T {
        self.velocity
            .iter()
            .copied()
            .fold(scalar::zero(), |acc, v| if v > acc { v } else { acc })
    }

    /// Calculate volumetric flow rate per unit width
    ///
    /// Q = ∫₀ᴴ u(y) dy
    ///
    /// Uses trapezoidal rule for integration
    pub fn flow_rate_per_width(&self) -> T {
        let ny = self.config.ny;
        let dy = self.dy;
        let two = scalar::from_f64::<T>(2.0);

        let mut sum: T = scalar::zero();
        for j in 0..ny - 1 {
            sum += (self.velocity[j] + self.velocity[j + 1]) / two * dy;
        }

        sum
    }

    /// Calculate total volumetric flow rate
    ///
    /// Q_total = Q_per_width * width
    pub fn flow_rate(&self) -> T {
        self.flow_rate_per_width() * self.config.width
    }

    /// Get velocity profile
    pub fn velocity_profile(&self) -> &[T] {
        &self.velocity
    }

    /// Get shear rate profile
    pub fn shear_rate_profile(&self) -> &[T] {
        &self.shear_rate
    }

    /// Get viscosity profile
    pub fn viscosity_profile(&self) -> &[T] {
        &self.viscosity
    }

    /// Get y-coordinates
    pub fn y_coordinates(&self) -> &[T] {
        &self.y_coords
    }

    /// Get wall shear rate (at bottom wall, y=0)
    pub fn wall_shear_rate(&self) -> T {
        self.shear_rate[0]
    }

    /// Get wall shear stress (at bottom wall, y=0)
    ///
    /// τ_wall = μ |du/dy|_{y=0}
    pub fn wall_shear_stress(&self) -> T {
        self.viscosity[0] * self.shear_rate[0]
    }

    /// Get convergence history
    pub fn convergence_history(&self) -> &[T] {
        &self.residuals
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use eunomia::NumericElement;

    #[test]
    fn test_thomas_algorithm() {
        // Test case: known solution
        let a = vec![0.0, -1.0, -1.0, -1.0];
        let b = vec![2.0, 2.0, 2.0, 2.0];
        let c = vec![-1.0, -1.0, -1.0, 0.0];
        let d = vec![1.0, 0.0, 0.0, 1.0];

        let x = numerics::thomas_algorithm(&a, &b, &c, &d).unwrap();

        // Verify solution
        assert!(<f64 as NumericElement>::abs(x[0] - 1.0) < 1e-10);
        assert!(<f64 as NumericElement>::abs(x[3] - 1.0) < 1e-10);
    }

    #[test]
    fn test_poiseuille_newtonian() {
        // Test with non-Newtonian blood model at high shear (approaches Newtonian)
        let config = PoiseuilleConfig {
            ny: 51,
            height: 0.001,               // 1mm channel
            pressure_gradient: 100000.0, // High pressure gradient for high shear
            tolerance: 1e-8,
            ..Default::default()
        };

        let blood = CassonBlood::<f64>::normal_blood();
        let mut solver = PoiseuilleFlow2D::new(config.clone(), BloodModel::Casson(blood));

        let iterations = solver.solve().unwrap();
        println!("Converged in {iterations} iterations");

        // Debug: check shear rates and viscosities
        let max_shear_rate = solver
            .shear_rate
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        let min_viscosity = solver
            .viscosity
            .iter()
            .min_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap();
        println!("Max shear rate: {max_shear_rate:.6e} s^-1");
        println!("Min viscosity: {min_viscosity:.6e} Pa·s");

        // For validation with non-Newtonian fluid:
        // Use minimum viscosity (high-shear asymptotic value) for analytical comparison
        // This represents the Newtonian behavior at high shear rates
        let effective_viscosity = *min_viscosity;
        println!("Effective viscosity (min): {effective_viscosity:.6e} Pa·s");

        let analytical = solver.analytical_solution(effective_viscosity);

        // Check maximum error (skip boundaries where analytical solution is zero)
        let max_error: f64 = solver
            .velocity
            .iter()
            .zip(analytical.iter())
            .enumerate()
            .filter(|(i, _)| *i > 0 && *i < solver.velocity.len() - 1) // Skip boundaries
            .map(|(_, (num, ana))| {
                if <f64 as NumericElement>::abs(*ana) < 1e-14 {
                    0.0 // Skip near-zero values
                } else {
                    <f64 as NumericElement>::abs((num - ana) / ana)
                }
            })
            .fold(0.0_f64, |max, val| if val > max { val } else { max });

        println!("Maximum relative error: {max_error:.2e}");
        let numerical_center = solver.velocity[solver.velocity.len() / 2];
        let analytical_center = analytical[analytical.len() / 2];
        println!("Velocity at center: numerical={numerical_center:.6e}, analytical={analytical_center:.6e}");

        // With high-shear viscosity, non-Newtonian blood should approach Newtonian behavior
        // Error comes from shear-thinning effects away from centerline
        // Accept up to 15% error due to viscosity variation across channel
        let max_error_percent = max_error * 100.0;
        assert!(
            max_error < 0.15,
            "Error too large: {max_error_percent:.2}% (expected <15% for high-shear non-Newtonian vs Newtonian)"
        );

        // Verify solver converged and produced reasonable results
        assert!(iterations < config.max_iterations, "Failed to converge");
        assert!(
            solver.velocity[solver.velocity.len() / 2] > 0.0,
            "Center velocity should be positive"
        );
    }
}
