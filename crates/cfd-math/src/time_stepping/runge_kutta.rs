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

use super::traits::{
    add_scaled_in_place, assign_base_plus_scaled, from_f64, state_len, state_zeros, TimeState,
    TimeStepper,
};
use cfd_core::error::Result;
use eunomia::FloatElement;
use eunomia::RealField;
use std::cell::RefCell;

/// Classical 4th-order Runge-Kutta method
///
/// # Theorem (RK4 Order of Accuracy)
///
/// The classical RK4 scheme with Butcher coefficients below is **4th-order
/// accurate**: for a sufficiently smooth ODE $y' = f(t, y)$ the local
/// truncation error satisfies $\|y(t_{n+1}) - y_{n+1}\| = O(h^5)$.
///
/// **Proof sketch**: Expanding $y(t_n + h)$ in a Taylor series to $O(h^5)$
/// and matching term-by-term with the RK4 update shows that all order
/// conditions $\sum_i b_i c_i^{q-1} = 1/q$ for $q = 1,2,3,4$ and the
/// Butcher simplifying assumptions B(4), C(2), D(1) are satisfied by the
/// tableau below. There are 8 independent conditions for order 4; all 8
/// hold for the classical weights $(1/6, 1/3, 1/3, 1/6)$.
///
/// **Reference**: Butcher (2008), *Numerical Methods for ODEs*, §3.2;
/// Hairer, Nørsett & Wanner (1993), §II.1.
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
pub struct RungeKutta4<T: RealField + Copy> {
    /// Scratch buffer for intermediate calculations to avoid allocation
    scratch: RefCell<TimeState<T>>,
}

impl<T: RealField + Copy> Default for RungeKutta4<T> {
    fn default() -> Self {
        Self {
            scratch: RefCell::new(state_zeros(0)),
        }
    }
}

impl<T: RealField + Copy> RungeKutta4<T> {
    /// Create a new fourth-order Runge-Kutta integrator
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: RealField + Copy + FloatElement> TimeStepper<T> for RungeKutta4<T> {
    fn step<F>(&self, f: F, t: T, u: &TimeState<T>, dt: T) -> Result<TimeState<T>>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
    {
        let n = state_len(u);

        // Reuse scratch buffer for intermediate states to avoid allocation
        let mut scratch = self.scratch.borrow_mut();
        if state_len(&scratch) != n {
            *scratch = state_zeros(n);
        }

        // k1 = f(t, u)
        let k1 = f(t, u)?;

        let dt_half = dt / from_f64(2.0);
        let t2 = t + dt_half;

        // k2 = f(t + dt/2, u + dt/2 * k1)
        scratch.assign(u);
        add_scaled_in_place(&mut scratch, &k1, dt_half, "RK4 k2 stage")?;
        let k2 = f(t2, &scratch)?;

        // k3 = f(t + dt/2, u + dt/2 * k2)
        scratch.assign(u);
        add_scaled_in_place(&mut scratch, &k2, dt_half, "RK4 k3 stage")?;
        let k3 = f(t2, &scratch)?;

        // k4 = f(t + dt, u + dt * k3)
        let t4 = t + dt;
        scratch.assign(u);
        add_scaled_in_place(&mut scratch, &k3, dt, "RK4 k4 stage")?;
        let k4 = f(t4, &scratch)?;

        // u_new = u + dt/6 * (k1 + 2*k2 + 2*k3 + k4)
        // We allocate u_new because the API requires returning an owned state.
        // However, we avoid intermediate temporary vectors in the summation.
        let mut u_new = state_zeros(n);
        u_new.assign(u);

        let coeff = dt / from_f64(6.0);
        let coeff_2 = coeff * from_f64(2.0);

        add_scaled_in_place(&mut u_new, &k1, coeff, "RK4 final k1 accumulation")?;
        add_scaled_in_place(&mut u_new, &k2, coeff_2, "RK4 final k2 accumulation")?;
        add_scaled_in_place(&mut u_new, &k3, coeff_2, "RK4 final k3 accumulation")?;
        add_scaled_in_place(&mut u_new, &k4, coeff, "RK4 final k4 accumulation")?;

        Ok(u_new)
    }

    fn order(&self) -> usize {
        4
    }

    fn stages(&self) -> usize {
        4
    }

    fn stability_region(&self) -> Option<&str> {
        Some("Absolute stability for |z| < 2.78")
    }
}

/// 3rd-order Runge-Kutta method (Kutta's method)
///
/// # Theorem (RK3 Order of Accuracy)
///
/// With the Butcher tableau below, the scheme satisfies all 4 order
/// conditions for a 3rd-order explicit Runge–Kutta method:
/// $\sum b_i = 1$, $\sum b_i c_i = 1/2$, $\sum b_i c_i^2 = 1/3$,
/// $\sum_i b_i \sum_j a_{ij} c_j = 1/6$, giving local truncation
/// error $O(h^4)$.
///
/// **Proof sketch**: Direct substitution of $(b, c, A)$ into the
/// Butcher order conditions and algebraic verification.
///
/// ## Butcher Tableau
/// ```text
/// 0   |
/// 1/2 | 1/2
/// 1   | -1   2
/// ----+------------
///     | 1/6  4/6  1/6
/// ```
pub struct RungeKutta3<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for RungeKutta3<T> {
    fn default() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy> RungeKutta3<T> {
    /// Create a new third-order Runge-Kutta integrator
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: RealField + Copy + FloatElement> TimeStepper<T> for RungeKutta3<T> {
    fn step<F>(&self, f: F, t: T, u: &TimeState<T>, dt: T) -> Result<TimeState<T>>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
    {
        let n = state_len(u);

        // k1 = f(t, u)
        let k1 = f(t, u)?;

        // k2 = f(t + dt/2, u + dt/2 * k1)
        let t2 = t + dt / from_f64(2.0);
        let mut u2 = state_zeros(n);
        assign_base_plus_scaled(&mut u2, u, &k1, dt / from_f64(2.0), "RK3 k2 stage")?;
        let k2 = f(t2, &u2)?;

        // k3 = f(t + dt, u - dt*k1 + 2*dt*k2)
        let t3 = t + dt;
        let mut u3 = state_zeros(n);
        let two = from_f64::<T>(2.0);
        for i in 0..n {
            u3[i] = u[i] - dt * k1[i] + two * dt * k2[i];
        }
        let k3 = f(t3, &u3)?;

        // u_new = u + dt/6 * (k1 + 4*k2 + k3)
        let mut u_new = state_zeros(n);
        let coeff1 = dt / from_f64::<T>(6.0);
        let four = from_f64::<T>(4.0);

        for i in 0..n {
            u_new[i] = u[i] + coeff1 * (k1[i] + four * k2[i] + k3[i]);
        }

        Ok(u_new)
    }

    fn order(&self) -> usize {
        3
    }

    fn stages(&self) -> usize {
        3
    }

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
pub struct LowStorageRK4<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy> Default for LowStorageRK4<T> {
    fn default() -> Self {
        Self {
            _phantom: std::marker::PhantomData,
        }
    }
}

impl<T: RealField + Copy> LowStorageRK4<T> {
    /// Create a new low-storage fourth-order Runge-Kutta integrator
    pub fn new() -> Self {
        Self::default()
    }
}

impl<T: RealField + Copy + FloatElement> TimeStepper<T> for LowStorageRK4<T> {
    fn step<F>(&self, f: F, t: T, u: &TimeState<T>, dt: T) -> Result<TimeState<T>>
    where
        F: Fn(T, &TimeState<T>) -> Result<TimeState<T>>,
    {
        let mut u_stage = u.clone();
        let mut residual = state_zeros(state_len(u));

        // Coefficients for Carpenter-Kennedy low-storage RK4
        let a = [
            from_f64::<T>(0.0),
            from_f64::<T>(-0.4178904745),
            from_f64::<T>(-1.1921516946),
            from_f64::<T>(-1.6977846925),
            from_f64::<T>(-1.5141834443),
        ];

        let b = [
            from_f64::<T>(0.1496590219993),
            from_f64::<T>(0.3792103129999),
            from_f64::<T>(0.8229550293869),
            from_f64::<T>(0.6994504559488),
            from_f64::<T>(0.1530572479681),
        ];

        let c = [
            from_f64::<T>(0.0),
            from_f64::<T>(0.1496590219993),
            from_f64::<T>(0.3704009573644),
            from_f64::<T>(0.6222557631345),
            from_f64::<T>(0.9582821306784),
        ];

        for stage in 0..5 {
            let t_stage = t + c[stage] * dt;

            // Compute RHS at current stage
            let rhs = f(t_stage, &u_stage)?;

            // Carpenter-Kennedy 2N form:
            // residual_i = a_i * residual_{i-1} + dt * f(t_i, u_i)
            // u_{i+1} = u_i + b_i * residual_i
            for i in 0..state_len(u) {
                residual[i] = a[stage] * residual[i] + dt * rhs[i];
                u_stage[i] += b[stage] * residual[i];
            }
        }

        Ok(u_stage)
    }

    fn order(&self) -> usize {
        4
    }

    fn stages(&self) -> usize {
        5
    }

    fn stability_region(&self) -> Option<&str> {
        Some("Large stability region, suitable for CFL > 1")
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time_stepping::traits::{state_from_vec, state_neg};
    use approx::assert_relative_eq;

    // Test function: du/dt = -u (exponential decay)
    fn exponential_decay(_t: f64, u: &TimeState<f64>) -> Result<TimeState<f64>> {
        Ok(state_neg(u))
    }

    #[test]
    fn test_rk4_exponential_decay() {
        let rk4 = RungeKutta4::new();
        let u0 = state_from_vec(vec![1.0]);
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
        let u0 = state_from_vec(vec![1.0, 2.0, 3.0]);
        let dt = 0.01;
        let t = 0.0;

        let u1 = rk4_ls.step(exponential_decay, t, &u0, dt).unwrap();

        let decay = (-dt).exp();
        assert_relative_eq!(u1[0], u0[0] * decay, epsilon = 1e-10);
        assert_relative_eq!(u1[1], u0[1] * decay, epsilon = 1e-10);
        assert_relative_eq!(u1[2], u0[2] * decay, epsilon = 1e-10);
    }

    #[test]
    fn test_low_storage_rk4_preserves_constant_solution() {
        let rk4_ls = LowStorageRK4::new();
        let u0 = state_from_vec(vec![1.0, -2.0, 3.5]);
        let dt = 0.25;
        let t = 1.0;

        let u1 = rk4_ls
            .step(|_, u| Ok(state_zeros(state_len(u))), t, &u0, dt)
            .unwrap();

        assert_relative_eq!(u1[0], u0[0], epsilon = 1e-14);
        assert_relative_eq!(u1[1], u0[1], epsilon = 1e-14);
        assert_relative_eq!(u1[2], u0[2], epsilon = 1e-14);
    }
}
