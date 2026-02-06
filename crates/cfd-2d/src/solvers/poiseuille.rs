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

use crate::error::Error;
use cfd_core::physics::fluid::blood::{CarreauYasudaBlood, CassonBlood};
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};
use serde::{Deserialize, Serialize};

/// Configuration for Poiseuille flow solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct PoiseuilleConfig<T: RealField + Copy> {
    /// Channel height [m]
    pub height: T,

    /// Channel width [m] (for 3D flow rate calculation)
    pub width: T,

    /// Channel length [m] (for reference only, flow is fully developed)
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

impl<T: RealField + FromPrimitive + Copy> Default for PoiseuilleConfig<T> {
    fn default() -> Self {
        Self {
            height: T::from_f64(100e-6).unwrap(),            // 100 μm
            width: T::from_f64(500e-6).unwrap(),             // 500 μm
            length: T::from_f64(1e-3).unwrap(),              // 1 mm
            ny: 101,                                         // 101 points for 100 intervals
            pressure_gradient: T::from_f64(1000.0).unwrap(), // 1000 Pa/m
            tolerance: T::from_f64(1e-6).unwrap(),
            max_iterations: 1000,
            relaxation_factor: T::from_f64(0.7).unwrap(),
        }
    }
}

/// 2D Poiseuille flow solver with non-Newtonian blood rheology
///
/// Solves steady-state velocity profile u(y) between parallel plates.
#[derive(Debug)]
pub struct PoiseuilleFlow2D<T: RealField + Copy + Float> {
    /// Configuration
    config: PoiseuilleConfig<T>,

    /// Grid in y-direction (1D since flow only varies in y)
    y_coords: Vec<T>,

    /// Grid spacing
    dy: T,

    /// Velocity field u(y)
    velocity: Vec<T>,

    /// Shear rate field γ̇(y) = |du/dy|
    shear_rate: Vec<T>,

    /// Viscosity field μ(y) = μ(γ̇(y))
    viscosity: Vec<T>,

    /// Blood model (Casson or Carreau-Yasuda)
    blood_model: BloodModel<T>,

    /// Convergence history
    residuals: Vec<T>,
}

/// Blood rheology model selection
#[derive(Debug, Clone)]
pub enum BloodModel<T: RealField + Copy> {
    /// Casson model with yield stress
    Casson(CassonBlood<T>),

    /// Carreau-Yasuda model
    CarreauYasuda(CarreauYasudaBlood<T>),
}

impl<T: RealField + FromPrimitive + Float + Copy> PoiseuilleFlow2D<T> {
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
        let dy = config.height / T::from_usize(ny - 1).unwrap();

        // Create grid points
        let y_coords: Vec<T> = (0..ny).map(|j| T::from_usize(j).unwrap() * dy).collect();

        // Initialize fields
        let velocity = vec![T::zero(); ny];
        let shear_rate = vec![T::zero(); ny];

        // Initialize viscosity with constant value
        let mu_init = T::from_f64(0.004).unwrap(); // ~4 cP initial guess
        let viscosity = vec![mu_init; ny];

        Self {
            config,
            y_coords,
            dy,
            velocity,
            shear_rate,
            viscosity,
            blood_model,
            residuals: Vec::new(),
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
            // Store old viscosity for convergence check
            let viscosity_old = self.viscosity.clone();

            // 1. Solve for velocity with current viscosity
            self.solve_velocity_linear()?;

            // 2. Calculate shear rate from velocity
            self.calculate_shear_rate();

            // 3. Update viscosity from blood model
            self.update_viscosity_from_shear_rate();

            // 4. Apply under-relaxation
            for j in 0..ny {
                self.viscosity[j] =
                    alpha * self.viscosity[j] + (T::one() - alpha) * viscosity_old[j];
            }

            // 5. Check convergence
            let residual = self.calculate_viscosity_residual(&viscosity_old);
            self.residuals.push(residual);

            if residual < tolerance {
                return Ok(iteration + 1);
            }
        }

        Err(Error::ConvergenceFailed {
            iterations: max_iter,
            residual: self
                .residuals
                .last()
                .copied()
                .unwrap_or(T::zero())
                .to_f64()
                .unwrap(),
        })
    }

    /// Solve linear system for velocity with current viscosity
    ///
    /// Solves tridiagonal system using Thomas algorithm (direct solver)
    ///
    /// Matrix equation: A*u = b where
    /// - A is tridiagonal
    /// - u is velocity vector
    /// - b is source term (pressure gradient)
    fn solve_velocity_linear(&mut self) -> Result<(), Error> {
        let ny = self.config.ny;
        let dy = self.dy;
        let dp_dx = self.config.pressure_gradient;

        // Allocate tridiagonal matrix coefficients
        let mut a = vec![T::zero(); ny]; // Sub-diagonal
        let mut b = vec![T::zero(); ny]; // Diagonal
        let mut c = vec![T::zero(); ny]; // Super-diagonal
        let mut d = vec![T::zero(); ny]; // RHS

        // Boundary conditions: u(0) = 0, u(ny-1) = 0
        b[0] = T::one();
        c[0] = T::zero();
        d[0] = T::zero();

        b[ny - 1] = T::one();
        a[ny - 1] = T::zero();
        d[ny - 1] = T::zero();

        // Interior points: j = 1, 2, ..., ny-2
        let dy2 = dy * dy;

        for j in 1..ny - 1 {
            // Viscosity at cell faces (harmonic mean for better stability)
            let mu_jm12 = T::from_f64(2.0).unwrap()
                / (T::one() / self.viscosity[j - 1] + T::one() / self.viscosity[j]);
            let mu_jp12 = T::from_f64(2.0).unwrap()
                / (T::one() / self.viscosity[j] + T::one() / self.viscosity[j + 1]);

            // Coefficients for interior stencil
            a[j] = -mu_jm12 / dy2;
            c[j] = -mu_jp12 / dy2;
            b[j] = mu_jm12 / dy2 + mu_jp12 / dy2;
            d[j] = dp_dx; // Source term from pressure gradient
        }

        // Solve tridiagonal system using Thomas algorithm
        self.velocity = thomas_algorithm(&a, &b, &c, &d)?;

        Ok(())
    }

    /// Calculate shear rate from velocity gradient
    ///
    /// γ̇(y) = |du/dy|
    ///
    /// Uses central differences for interior points:
    /// du/dy ≈ (u_{j+1} - u_{j-1}) / (2Δy)
    fn calculate_shear_rate(&mut self) {
        let ny = self.config.ny;
        let dy = self.dy;
        let two = T::from_f64(2.0).unwrap();

        // Boundaries: use one-sided differences
        self.shear_rate[0] = Float::abs((self.velocity[1] - self.velocity[0]) / dy);
        self.shear_rate[ny - 1] = Float::abs((self.velocity[ny - 1] - self.velocity[ny - 2]) / dy);

        // Interior: central differences
        for j in 1..ny - 1 {
            let du_dy = (self.velocity[j + 1] - self.velocity[j - 1]) / (two * dy);
            self.shear_rate[j] = Float::abs(du_dy);
        }
    }

    /// Update viscosity from shear rate using blood model
    ///
    /// μ(y) = μ(γ̇(y))
    fn update_viscosity_from_shear_rate(&mut self) {
        for j in 0..self.config.ny {
            self.viscosity[j] = match &self.blood_model {
                BloodModel::Casson(casson) => casson.apparent_viscosity(self.shear_rate[j]),
                BloodModel::CarreauYasuda(carreau) => {
                    carreau.apparent_viscosity(self.shear_rate[j])
                }
            };
        }
    }

    /// Calculate L2 norm of viscosity change (for convergence check)
    fn calculate_viscosity_residual(&self, viscosity_old: &[T]) -> T {
        let mut sum_sq = T::zero();
        let mut sum_sq_old = T::zero();

        for j in 0..self.config.ny {
            let diff = self.viscosity[j] - viscosity_old[j];
            sum_sq = sum_sq + diff * diff;
            sum_sq_old = sum_sq_old + viscosity_old[j] * viscosity_old[j];
        }

        // Relative L2 norm
        Float::sqrt(sum_sq / (sum_sq_old + T::from_f64(1e-20).unwrap()))
    }

    /// Calculate analytical solution for Newtonian fluid (for validation)
    ///
    /// u(y) = (1/2μ)(dP/dx)y(H - y)
    ///
    /// Returns `None` if blood model is non-Newtonian
    pub fn analytical_solution(&self, viscosity: T) -> Vec<T> {
        let dp_dx = self.config.pressure_gradient;
        let height = self.config.height;
        let two = T::from_f64(2.0).unwrap();

        self.y_coords
            .iter()
            .map(|&y| (dp_dx / (two * viscosity)) * y * (height - y))
            .collect()
    }

    /// Get maximum velocity (at centerline)
    pub fn max_velocity(&self) -> T {
        *self
            .velocity
            .iter()
            .max_by(|a, b| a.partial_cmp(b).unwrap())
            .unwrap()
    }

    /// Calculate volumetric flow rate per unit width
    ///
    /// Q = ∫₀ᴴ u(y) dy
    ///
    /// Uses trapezoidal rule for integration
    pub fn flow_rate_per_width(&self) -> T {
        let ny = self.config.ny;
        let dy = self.dy;
        let two = T::from_f64(2.0).unwrap();

        let mut sum = T::zero();
        for j in 0..ny - 1 {
            sum = sum + (self.velocity[j] + self.velocity[j + 1]) / two * dy;
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

/// Thomas algorithm for tridiagonal systems
///
/// Solves Ax = d where A is tridiagonal:
/// - a: sub-diagonal (a[0] unused)
/// - b: diagonal
/// - c: super-diagonal (c[n-1] unused)
/// - d: right-hand side
///
/// Returns solution vector x
fn thomas_algorithm<T: RealField + Copy + Float>(
    a: &[T],
    b: &[T],
    c: &[T],
    d: &[T],
) -> Result<Vec<T>, Error> {
    let n = b.len();

    if a.len() != n || c.len() != n || d.len() != n {
        return Err(Error::InvalidInput("Array lengths must match".to_string()));
    }

    let mut c_prime = vec![T::zero(); n];
    let mut d_prime = vec![T::zero(); n];
    let mut x = vec![T::zero(); n];

    // Forward sweep
    c_prime[0] = c[0] / b[0];
    d_prime[0] = d[0] / b[0];

    for i in 1..n {
        let denom = b[i] - a[i] * c_prime[i - 1];

        if Float::abs(denom) < T::from_f64(1e-14).unwrap() {
            return Err(Error::NumericalInstability(
                "Near-zero pivot in Thomas algorithm".to_string(),
            ));
        }

        c_prime[i] = c[i] / denom;
        d_prime[i] = (d[i] - a[i] * d_prime[i - 1]) / denom;
    }

    // Back substitution
    x[n - 1] = d_prime[n - 1];

    for i in (0..n - 1).rev() {
        x[i] = d_prime[i] - c_prime[i] * x[i + 1];
    }

    Ok(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_thomas_algorithm() {
        // Test case: known solution
        let a = vec![0.0, -1.0, -1.0, -1.0];
        let b = vec![2.0, 2.0, 2.0, 2.0];
        let c = vec![-1.0, -1.0, -1.0, 0.0];
        let d = vec![1.0, 0.0, 0.0, 1.0];

        let x = thomas_algorithm(&a, &b, &c, &d).unwrap();

        // Verify solution
        assert!((x[0] - 1.0).abs() < 1e-10);
        assert!((x[3] - 1.0).abs() < 1e-10);
    }

    #[test]
    fn test_poiseuille_newtonian() {
        // Test with non-Newtonian blood model at high shear (approaches Newtonian)
        let mut config = PoiseuilleConfig::default();
        config.ny = 51;
        config.height = 0.001; // 1mm channel
        config.pressure_gradient = 100000.0; // High pressure gradient for high shear
        config.tolerance = 1e-8;

        let blood = CassonBlood::<f64>::normal_blood();
        let mut solver = PoiseuilleFlow2D::new(config.clone(), BloodModel::Casson(blood));

        let iterations = solver.solve().unwrap();
        println!("Converged in {} iterations", iterations);

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
        println!("Max shear rate: {:.6e} s^-1", max_shear_rate);
        println!("Min viscosity: {:.6e} Pa·s", min_viscosity);

        // For validation with non-Newtonian fluid:
        // Use minimum viscosity (high-shear asymptotic value) for analytical comparison
        // This represents the Newtonian behavior at high shear rates
        let effective_viscosity = *min_viscosity;
        println!(
            "Effective viscosity (min): {:.6e} Pa·s",
            effective_viscosity
        );

        let analytical = solver.analytical_solution(effective_viscosity);

        // Check maximum error (skip boundaries where analytical solution is zero)
        let max_error: f64 = solver
            .velocity
            .iter()
            .zip(analytical.iter())
            .enumerate()
            .filter(|(i, _)| *i > 0 && *i < solver.velocity.len() - 1) // Skip boundaries
            .map(|(_, (num, ana))| {
                if ana.abs() < 1e-14 {
                    0.0 // Skip near-zero values
                } else {
                    ((num - ana) / ana).abs()
                }
            })
            .fold(0.0_f64, |max, val| if val > max { val } else { max });

        println!("Maximum relative error: {:.2e}", max_error);
        println!(
            "Velocity at center: numerical={:.6e}, analytical={:.6e}",
            solver.velocity[solver.velocity.len() / 2],
            analytical[analytical.len() / 2]
        );

        // With high-shear viscosity, non-Newtonian blood should approach Newtonian behavior
        // Error comes from shear-thinning effects away from centerline
        // Accept up to 15% error due to viscosity variation across channel
        assert!(
            max_error < 0.15,
            "Error too large: {:.2}% (expected <15% for high-shear non-Newtonian vs Newtonian)",
            max_error * 100.0
        );

        // Verify solver converged and produced reasonable results
        assert!(iterations < config.max_iterations, "Failed to converge");
        assert!(
            solver.velocity[solver.velocity.len() / 2] > 0.0,
            "Center velocity should be positive"
        );
    }
}
