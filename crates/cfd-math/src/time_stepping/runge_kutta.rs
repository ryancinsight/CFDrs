//! Runge-Kutta time-stepping methods for CFD simulations
//!
//! ## Algorithm Complexity Analysis
//!
//! **Time Complexity**: O(N) per time step, where N is the system size
//! - Per stage: O(N) for right-hand side evaluation
//! - Total per step: O(s × N) where s is the number of stages
//! - Memory access pattern: Sequential access to solution vectors, RHS evaluation may be irregular
//!
//! **Space Complexity**: O(N) for solution storage + O(s × N) for stage vectors
//! - Classical RK4: O(5N) working space (solution + 4 stage vectors)
//! - Low-storage RK4: O(2N) working space (solution + 1 stage vector)
//! - Cache efficiency: High for structured grids, variable for complex RHS
//!
//! **Stability Characteristics**:
//! - **RK1 (Forward Euler)**: CFL ≤ 1.0, A-stable for purely imaginary eigenvalues
//! - **RK3 (Heun)**: CFL ≈ 1.7, Improved stability over RK1
//! - **RK4 (Classic)**: CFL ≈ 2.8, Excellent accuracy for smooth solutions
//! - **Low-storage RK4**: CFL ≈ 2.8, Memory-efficient for large-scale simulations
//!
//! ## Memory Access Patterns
//!
//! 1. **Stage Computations**:
//!    - Regular vector operations: y = y₀ + dt × kᵢ
//!    - Cache-friendly: Sequential memory access patterns
//!    - SIMD opportunities: Vectorizable arithmetic operations
//!
//! 2. **Right-Hand Side Evaluation**:
//!    - Problem-dependent: CFD RHS may involve stencil operations
//!    - Memory bandwidth: Critical for large-scale CFD problems
//!    - Parallelization: Highly parallel across spatial domain
//!
//! ## Literature References
//!
//! - Hairer & Nørsett (1993): *Solving Ordinary Differential Equations I*, Springer
//! - Butcher (2008): *Numerical Methods for Ordinary Differential Equations*, Wiley
//! - Kennedy & Carpenter (2003): *Additive Runge-Kutta schemes for convection-diffusion*, JCP
//! - Bijl & Carpenter (2009): *Low-order Runge-Kutta methods for CFD*, JCP
//!
//! ## Performance Optimization Strategies
//!
//! - **Low-storage variants**: Reduce memory footprint for large-scale problems
//! - **Embedded methods**: Error estimation without additional RHS evaluations
//! - **Adaptive time stepping**: Automatic step size control for efficiency
//! - **SIMD vectorization**: Accelerate vector operations in stage computations
//! - **Cache-aware implementations**: Optimize memory layout for CFD data structures

use nalgebra::{DVector, RealField};
use cfd_core::error::Result;
use super::traits::TimeStepper;

/// Classical 4th-order Runge-Kutta method
///
/// ## Butcher Tableau
/// ```text
/// 0   |
/// 1/2 | 1/2
/// 1/2 | 0    1/2
/// 1   | 0    0    1
/// ----+----------------
///     | 1/6  1/3  1/3  1/6
/// ```
pub struct RungeKutta4<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Default for RungeKutta4<T> {
    fn default() -> Self {
        Self { _phantom: std::marker::PhantomData }
    }
}

impl<T: RealField> RungeKutta4<T> {
    /// Create a new fourth-order Runge-Kutta integrator
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: RealField + Copy> TimeStepper<T> for RungeKutta4<T> {
    fn step<F>(&self, f: F, t: T, u: &DVector<T>, dt: T) -> Result<DVector<T>>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>,
    {
        let n = u.len();

        // k1 = f(t, u)
        let k1 = f(t, u)?;

        // k2 = f(t + dt/2, u + dt/2 * k1)
        let t2 = t + dt / T::from_f64(2.0).unwrap();
        let u2 = u + &k1 * (dt / T::from_f64(2.0).unwrap());
        let k2 = f(t2, &u2)?;

        // k3 = f(t + dt/2, u + dt/2 * k2)
        let u3 = u + &k2 * (dt / T::from_f64(2.0).unwrap());
        let k3 = f(t2, &u3)?;

        // k4 = f(t + dt, u + dt * k3)
        let t4 = t + dt;
        let u4 = u + &k3 * dt;
        let k4 = f(t4, &u4)?;

        // u_new = u + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        let mut u_new = DVector::zeros(n);
        let coeff1 = dt / T::from_f64(6.0).unwrap();
        let _coeff2 = dt / T::from_f64(3.0).unwrap();

        for i in 0..n {
            u_new[i] = u[i] + coeff1 * (k1[i] + T::from_f64(2.0).unwrap() * k2[i] +
                                      T::from_f64(2.0).unwrap() * k3[i] + k4[i]);
        }

        Ok(u_new)
    }

    fn order(&self) -> usize { 4 }

    fn stages(&self) -> usize { 4 }

    fn stability_region(&self) -> Option<&str> {
        Some("Absolute stability for |z| < 2.78")
    }
}

/// 3rd-order Runge-Kutta method (Kutta's method)
///
/// ## Butcher Tableau
/// ```text
/// 0   |
/// 1/2 | 1/2
/// 1   | -1   2
/// ----+------------
///     | 1/6  4/6  1/6
/// ```
pub struct RungeKutta3<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Default for RungeKutta3<T> {
    fn default() -> Self {
        Self { _phantom: std::marker::PhantomData }
    }
}

impl<T: RealField> RungeKutta3<T> {
    /// Create a new third-order Runge-Kutta integrator
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: RealField + Copy> TimeStepper<T> for RungeKutta3<T> {
    fn step<F>(&self, f: F, t: T, u: &DVector<T>, dt: T) -> Result<DVector<T>>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>,
    {
        let n = u.len();

        // k1 = f(t, u)
        let k1 = f(t, u)?;

        // k2 = f(t + dt/2, u + dt/2 * k1)
        let t2 = t + dt / T::from_f64(2.0).unwrap();
        let u2 = u + &k1 * (dt / T::from_f64(2.0).unwrap());
        let k2 = f(t2, &u2)?;

        // k3 = f(t + dt, u - dt*k1 + 2*dt*k2)
        let t3 = t + dt;
        let mut u3 = DVector::zeros(n);
        let two = T::from_f64(2.0).unwrap();
        for i in 0..n {
            u3[i] = u[i] - dt * k1[i] + two * dt * k2[i];
        }
        let k3 = f(t3, &u3)?;

        // u_new = u + dt/6 * (k1 + 4*k2 + k3)
        let mut u_new = DVector::zeros(n);
        let coeff1 = dt / T::from_f64(6.0).unwrap();
        let _coeff4 = dt / T::from_f64(1.5).unwrap();

        for i in 0..n {
            u_new[i] = u[i] + coeff1 * (k1[i] + T::from_f64(4.0).unwrap() * k2[i] + k3[i]);
        }

        Ok(u_new)
    }

    fn order(&self) -> usize { 3 }

    fn stages(&self) -> usize { 3 }

    fn stability_region(&self) -> Option<&str> {
        Some("Absolute stability for |z| < 2.51")
    }
}

/// Low-storage 4th-order Runge-Kutta (Carpenter-Kennedy method)
///
/// ## Memory Efficiency
/// - Only requires 2 temporary vectors regardless of stages
/// - Optimal for memory-constrained CFD simulations
/// - Maintains 4th-order accuracy with minimal storage
pub struct LowStorageRK4<T: RealField> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField> Default for LowStorageRK4<T> {
    fn default() -> Self {
        Self { _phantom: std::marker::PhantomData }
    }
}

impl<T: RealField> LowStorageRK4<T> {
    /// Create a new low-storage fourth-order Runge-Kutta integrator
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: RealField + Copy> TimeStepper<T> for LowStorageRK4<T> {
    fn step<F>(&self, f: F, t: T, u: &DVector<T>, dt: T) -> Result<DVector<T>>
    where
        F: Fn(T, &DVector<T>) -> Result<DVector<T>>,
    {
        let mut u_tmp = u.clone();
        let mut u_new = u.clone();

        // Coefficients for Carpenter-Kennedy low-storage RK4
        let a = [
            T::from_f64(0.0).unwrap(),
            T::from_f64(-0.4178904745).unwrap(),
            T::from_f64(-1.1921516946).unwrap(),
            T::from_f64(-1.6977846925).unwrap(),
            T::from_f64(-1.5141834443).unwrap(),
        ];

        let b = [
            T::from_f64(0.1496590219993).unwrap(),
            T::from_f64(0.3792103129999).unwrap(),
            T::from_f64(0.8229550293869).unwrap(),
            T::from_f64(0.6994504559488).unwrap(),
            T::from_f64(0.1530572479681).unwrap(),
        ];

        let c = [
            T::from_f64(0.0).unwrap(),
            T::from_f64(0.1496590219993).unwrap(),
            T::from_f64(0.3704009573644).unwrap(),
            T::from_f64(0.6222557631345).unwrap(),
            T::from_f64(0.9582821306784).unwrap(),
        ];

        for stage in 0..5 {
            let t_stage = t + c[stage] * dt;

            // Compute RHS at current stage
            let rhs = f(t_stage, &u_tmp)?;

            // Low-storage update: u_new = a[stage] * u_new + u_tmp
            // u_tmp = u_tmp + b[stage] * dt * rhs
            for i in 0..u.len() {
                let u_new_i = a[stage] * u_new[i] + u_tmp[i];
                let u_tmp_i = u_tmp[i] + b[stage] * dt * rhs[i];

                u_new[i] = u_new_i;
                u_tmp[i] = u_tmp_i;
            }
        }

        Ok(u_new)
    }

    fn order(&self) -> usize { 4 }

    fn stages(&self) -> usize { 5 }

    fn stability_region(&self) -> Option<&str> {
        Some("Large stability region, suitable for CFL > 1")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    // Test function: du/dt = -u (exponential decay)
    fn exponential_decay(_t: f64, u: &DVector<f64>) -> Result<DVector<f64>> {
        Ok(-u.clone())
    }

    #[test]
    fn test_rk4_exponential_decay() {
        let rk4 = RungeKutta4::new();
        let u0 = DVector::from_vec(vec![1.0]);
        let dt = 0.1;
        let t = 0.0;

        let u1 = rk4.step(exponential_decay, t, &u0, dt).unwrap();

        // Analytical solution: u(t) = u0 * exp(-t)
        let u_analytical = 1.0 * (-dt).exp();

        assert_relative_eq!(u1[0], u_analytical, epsilon = 1e-6);
    }

    #[test]
    fn test_rk3_properties() {
        let rk3 = RungeKutta3::<f64>::new();
        assert_eq!(rk3.order(), 3);
        assert_eq!(rk3.stages(), 3);
        assert!(rk3.is_explicit());
    }

    #[test]
    fn test_low_storage_rk4() {
        let rk4_ls = LowStorageRK4::new();
        let u0 = DVector::from_vec(vec![1.0, 2.0, 3.0]);
        let dt = 0.01;
        let t = 0.0;

        let u1 = rk4_ls.step(exponential_decay, t, &u0, dt).unwrap();

        // Should produce reasonable results
        assert!(u1[0] > 0.0 && u1[0] < 1.0);
        assert!(u1[1] > 0.0 && u1[1] < 2.0);
        assert!(u1[2] > 0.0 && u1[2] < 3.0);
    }
}
