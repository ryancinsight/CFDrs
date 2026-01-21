//! Runge-Kutta-Chebyshev Methods for Stiff Problems
//!
//! Runge-Kutta-Chebyshev (RKC) methods are explicit schemes designed for solving
//! stiff diffusion-dominated problems. They use Chebyshev polynomials to extend
//! the stability region along the negative real axis.
//!
//! ## Mathematical Formulation (Second-Order RKC)
//!
//! Based on Sommeijer, B. P., Shampine, L. F., & Verwer, J. G. (1997).
//! "RKC: An explicit solver for parabolic PDEs". Journal of Computational and Applied Mathematics.
//!
//! The scheme uses the three-term recurrence relation of Chebyshev polynomials.
//!
//! ### Algorithm
//!
//! For $j = 2, \dots, s$:
//!
//! $Y_0 = y_n$
//! $Y_1 = Y_0 + \tilde{\mu}_1 \Delta t F(Y_0)$
//! $Y_j = (1 - \mu_j - \nu_j) Y_0 + \mu_j Y_{j-1} + \nu_j Y_{j-2} + \tilde{\mu}_j \Delta t F(Y_{j-1}) + \tilde{\gamma}_j \Delta t F(Y_0)$
//! $y_{n+1} = Y_s$
//!
//! ### Coefficients
//!
//! $\epsilon = 2/13$ (damping parameter)
//! $w_0 = 1 + \epsilon/s^2$
//!
//! $b_j = \frac{T''_j(w_0)}{(T'_j(w_0))^2}$
//!
//! $\mu_j = \frac{2 w_0 b_j}{b_{j-1}}$
//! $\nu_j = - \frac{b_j}{b_{j-2}}$
//! $\tilde{\mu}_j = \frac{2 w_1 b_j}{b_{j-1}}$
//! $\tilde{\gamma}_j = - (1 - b_{j-1} T_j(w_0)) \tilde{\mu}_j$
//!
//! Note: The exact recurrence coefficients for the first few stages and the relation to $T_j$
//! follow the properties of shifted Chebyshev polynomials.
//!
//! The stability limit is $\beta \approx (w_0 + 1) s^2 / 2 \approx 0.8 s^2$ for damped RKC.

use crate::error::Result;
use nalgebra::{DVector, RealField};
use num_traits::FromPrimitive;

/// RKC method configuration
#[derive(Debug, Clone, Copy)]
pub struct RkcConfig<T: RealField + Copy> {
    /// Number of stages (s)
    pub num_stages: usize,
    /// Damping parameter (epsilon), typically 2/13 ≈ 0.15
    pub damping: T,
    /// Absolute tolerance for error control
    pub atol: T,
    /// Relative tolerance for error control
    pub rtol: T,
    /// Safety factor for time step adaptation
    pub safety_factor: T,
}

impl<T: RealField + Copy + FromPrimitive> Default for RkcConfig<T> {
    fn default() -> Self {
        Self {
            num_stages: 10,
            damping: T::from_f64(2.0 / 13.0).unwrap(),
            atol: T::from_f64(1e-8).unwrap(),
            rtol: T::from_f64(1e-6).unwrap(),
            safety_factor: T::from_f64(0.9).unwrap(),
        }
    }
}

/// Runge-Kutta-Chebyshev method implementation
#[derive(Debug, Clone)]
pub struct RungeKuttaChebyshev<T: RealField + Copy> {
    config: RkcConfig<T>,
    /// Precomputed stage coefficients
    coeffs: Vec<RkcStageCoeffs<T>>,
}

#[derive(Debug, Clone, Copy)]
struct RkcStageCoeffs<T> {
    mu: T,
    nu: T,
    mu_tilde: T,
    gamma_tilde: T,
    c: T,
}

/// Right-hand side function trait
pub trait RhsFunction<T: RealField + Copy> {
    /// Evaluate the right-hand side f(t, y)
    fn evaluate(&self, t: T, y: &DVector<T>) -> Result<DVector<T>>;
}

impl<T: RealField + Copy + FromPrimitive> RungeKuttaChebyshev<T> {
    /// Create RKC method with default configuration
    pub fn new() -> Self {
        Self::with_config(RkcConfig::default())
    }

    /// Create RKC method with custom configuration
    pub fn with_config(config: RkcConfig<T>) -> Self {
        let coeffs = Self::compute_coefficients(config.num_stages, config.damping);
        Self { config, coeffs }
    }

    /// Compute Chebyshev polynomials T_k(x), T'_k(x), T''_k(x)
    fn chebyshev_eval(k: usize, x: T) -> (T, T, T) {
        if k == 0 {
            return (T::one(), T::zero(), T::zero());
        }
        if k == 1 {
            return (x, T::one(), T::zero());
        }

        let two = T::from_f64(2.0).unwrap();
        let mut t_prev = x;
        let mut t_prev2 = T::one();
        let mut tp_prev = T::one();
        let mut tp_prev2 = T::zero();
        let mut tpp_prev = T::zero();
        let mut tpp_prev2 = T::zero();

        let mut t_curr = T::zero();
        let mut tp_curr = T::zero();
        let mut tpp_curr = T::zero();

        for _ in 2..=k {
            // T_j(x) = 2x T_{j-1}(x) - T_{j-2}(x)
            t_curr = two * x * t_prev - t_prev2;
            // T'_j(x) = 2 T_{j-1}(x) + 2x T'_{j-1}(x) - T'_{j-2}(x)
            tp_curr = two * t_prev + two * x * tp_prev - tp_prev2;
            // T''_j(x) = 4 T'_{j-1}(x) + 2x T''_{j-1}(x) - T''_{j-2}(x)
            tpp_curr = T::from_f64(4.0).unwrap() * tp_prev + two * x * tpp_prev - tpp_prev2;

            t_prev2 = t_prev;
            t_prev = t_curr;
            tp_prev2 = tp_prev;
            tp_prev = tp_curr;
            tpp_prev2 = tpp_prev;
            tpp_prev = tpp_curr;
        }

        (t_curr, tp_curr, tpp_curr)
    }

    /// Precompute RKC coefficients for given number of stages and damping
    fn compute_coefficients(s: usize, epsilon: T) -> Vec<RkcStageCoeffs<T>> {
        let mut coeffs = Vec::with_capacity(s + 1);

        // Stage 0 and 1 are special or unused in the recurrence loop (j=2..s)
        // We define dummy values for index 0 and 1 to align with 1-based indexing logic
        coeffs.push(RkcStageCoeffs {
            mu: T::zero(),
            nu: T::zero(),
            mu_tilde: T::zero(),
            gamma_tilde: T::zero(),
            c: T::zero(),
        }); // Index 0

        // w0 = 1 + epsilon / s^2
        let s_sq = T::from_usize(s * s).unwrap();
        let w0 = T::one() + epsilon / s_sq;

        // Calculate b_j = T''_j(w0) / (T'_j(w0))^2
        let mut b = vec![T::zero(); s + 1];
        for j in 2..=s {
            let (_, tp, tpp) = Self::chebyshev_eval(j, w0);
            b[j] = tpp / (tp * tp);
        }
        // b_0 = b_2, b_1 = b_2 according to standard RKC
        b[0] = b[2];
        b[1] = b[2];

        // Calculate w1 = T'_s(w0) / T''_s(w0)
        let (_, tps, tpps) = Self::chebyshev_eval(s, w0);
        let w1 = tps / tpps;

        // First stage coefficient (tilde_mu_1)
        // Y_1 = Y_0 + tilde_mu_1 * dt * F_0
        // tilde_mu_1 = w1 * b_1
        let mu_tilde_1 = w1 * b[1];
        coeffs.push(RkcStageCoeffs {
            mu: T::zero(),
            nu: T::zero(),
            mu_tilde: mu_tilde_1,
            gamma_tilde: T::zero(),
            c: mu_tilde_1,
        }); // Index 1

        // Subsequent stages j = 2..s
        for j in 2..=s {
            let two = T::from_f64(2.0).unwrap();

            let mu_j = two * w0 * b[j] / b[j - 1];
            let nu_j = -b[j] / b[j - 2];
            let mu_tilde_j = two * w1 * b[j] / b[j - 1];

            let (t_j_minus1, _, _) = Self::chebyshev_eval(j - 1, w0);
            let gamma_tilde_j = -(T::one() - b[j - 1] * t_j_minus1) * mu_tilde_j;
            let c_prev = coeffs[j - 1].c;
            let c_prev2 = coeffs[j - 2].c;
            let c_j = mu_j * c_prev + nu_j * c_prev2 + mu_tilde_j + gamma_tilde_j;

            coeffs.push(RkcStageCoeffs {
                mu: mu_j,
                nu: nu_j,
                mu_tilde: mu_tilde_j,
                gamma_tilde: gamma_tilde_j,
                c: c_j,
            });
        }

        coeffs
    }

    /// Solve ODE system using RKC method
    pub fn step<F: RhsFunction<T>>(
        &self,
        rhs: &F,
        t0: T,
        y0: &DVector<T>,
        dt: T,
    ) -> Result<DVector<T>> {
        let n = y0.len();

        // Y_0 = y_n
        let mut y_prev2 = y0.clone(); // Y_{j-2} (initially Y_0)

        // Evaluate F_0 = F(t0, Y_0)
        let f0 = rhs.evaluate(t0, &y_prev2)?;

        // Y_1 = Y_0 + mu_tilde_1 * dt * F_0
        let mu_tilde_1 = self.coeffs[1].mu_tilde;
        let mut y_prev1 = &y_prev2 + &f0 * (mu_tilde_1 * dt); // Y_{j-1} (initially Y_1)

        // Stages 2 to s
        for j in 2..=self.config.num_stages {
            let coeff = &self.coeffs[j];

            let t_stage = t0 + self.coeffs[j - 1].c * dt;
            let f_prev = rhs.evaluate(t_stage, &y_prev1)?;

            // Y_j = (1 - mu - nu) Y_0 + mu Y_{j-1} + nu Y_{j-2} + mu_tilde dt F_{j-1} + gamma_tilde dt F_0
            let term_y0 = T::one() - coeff.mu - coeff.nu;

            let mut y_curr = DVector::zeros(n);
            for i in 0..n {
                y_curr[i] = term_y0 * y0[i]
                    + coeff.mu * y_prev1[i]
                    + coeff.nu * y_prev2[i]
                    + coeff.mu_tilde * dt * f_prev[i]
                    + coeff.gamma_tilde * dt * f0[i];
            }

            // Shift
            y_prev2 = y_prev1;
            y_prev1 = y_curr;
        }

        Ok(y_prev1)
    }

    /// Solve ODE system with adaptive time stepping
    pub fn solve_adaptive<F: RhsFunction<T>>(
        &self,
        rhs: &F,
        t0: T,
        y0: &DVector<T>,
        t_final: T,
        dt_initial: T,
    ) -> Result<(DVector<T>, T)> {
        let mut t = t0;
        let mut y = y0.clone();
        let mut dt = dt_initial;

        while t < t_final {
            if t + dt > t_final {
                dt = t_final - t;
            }

            let mut accepted = false;
            let mut attempts = 0;
            let max_attempts = 10;
            while !accepted {
                attempts += 1;
                if attempts > max_attempts {
                    return Err(cfd_core::error::Error::Convergence(
                        cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded { max: max_attempts },
                    ));
                }

                let y_full = self.step(rhs, t, &y, dt)?;
                let dt_half = dt / T::from_f64(2.0).unwrap();
                let y_half = self.step(rhs, t, &y, dt_half)?;
                let y_half_full = self.step(rhs, t + dt_half, &y_half, dt_half)?;

                let n = y.len();
                let mut err_sum = T::zero();
                for i in 0..n {
                    let scale = self.config.atol
                        + self.config.rtol * y[i].abs().max(y_half_full[i].abs());
                    let denom = if scale > T::zero() { scale } else { T::one() };
                    let diff = (y_half_full[i] - y_full[i]).abs()
                        / (T::from_f64(3.0).unwrap() * denom);
                    err_sum += diff * diff;
                }
                let n_t = T::from_usize(n).unwrap();
                let err_norm = (err_sum / n_t).sqrt();

                if err_norm <= T::one() {
                    y = y_half_full;
                    t += dt;
                    let factor = if err_norm == T::zero() {
                        T::from_f64(5.0).unwrap()
                    } else {
                        let exponent = T::from_f64(1.0 / 3.0).unwrap();
                        self.config.safety_factor * err_norm.powf(-exponent)
                    };
                    let factor = factor
                        .max(T::from_f64(0.1).unwrap())
                        .min(T::from_f64(5.0).unwrap());
                    dt *= factor;
                    accepted = true;
                } else {
                    let exponent = T::from_f64(1.0 / 3.0).unwrap();
                    let factor = self.config.safety_factor * err_norm.powf(-exponent);
                    let factor = factor
                        .max(T::from_f64(0.1).unwrap())
                        .min(T::from_f64(0.5).unwrap());
                    dt *= factor;
                    if t + dt > t_final {
                        dt = t_final - t;
                    }
                }
            }
        }

        Ok((y, dt))
    }

    /// Get the configuration of the RKC solver
    pub fn config(&self) -> &RkcConfig<T> {
        &self.config
    }
}

impl<T: RealField + Copy + FromPrimitive> Default for RungeKuttaChebyshev<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    struct ExponentialDecay<T: RealField + Copy> {
        lambda: T,
    }

    impl<T: RealField + Copy> ExponentialDecay<T> {
        fn new(lambda: T) -> Self {
            Self { lambda }
        }
    }

    impl<T: RealField + Copy + FromPrimitive> RhsFunction<T> for ExponentialDecay<T> {
        fn evaluate(&self, _t: T, y: &DVector<T>) -> Result<DVector<T>> {
            Ok(y.scale(-self.lambda))
        }
    }

    #[test]
    fn test_rkc_creation() {
        let rkc = RungeKuttaChebyshev::<f64>::new();
        assert_eq!(rkc.config().num_stages, 10);
    }

    #[test]
    fn test_exponential_decay() {
        let lambda: f64 = 1.0;
        let rhs = ExponentialDecay::new(lambda);
        let rkc = RungeKuttaChebyshev::<f64>::new();

        let y0 = DVector::from_vec(vec![1.0]);
        let dt = 0.1;
        let y_final = rkc.step(&rhs, 0.0, &y0, dt).unwrap();

        let analytical = (-lambda * dt).exp();
        assert_relative_eq!(y_final[0], analytical, epsilon = 1e-3);
    }

    #[test]
    fn test_stiff_problem() {
        let lambda: f64 = 100.0; // Stiff
        let rhs = ExponentialDecay::new(lambda);

        // Increase stages for stiffness
        let config = RkcConfig {
            num_stages: 20,
            ..RkcConfig::default()
        };
        let rkc = RungeKuttaChebyshev::with_config(config);

        let y0 = DVector::from_vec(vec![1.0]);
        let dt = 0.01; // lambda*dt = 1
        let y_final = rkc.step(&rhs, 0.0, &y0, dt).unwrap();

        let analytical = (-lambda * dt).exp();
        assert_relative_eq!(y_final[0], analytical, epsilon = 1e-3);
    }

    #[test]
    fn test_adaptive_solving() {
        let lambda: f64 = 1.0;
        let rhs = ExponentialDecay::new(lambda);
        let rkc = RungeKuttaChebyshev::<f64>::new();

        let y0 = DVector::from_vec(vec![1.0]);
        let (y_final, _) = rkc.solve_adaptive(&rhs, 0.0, &y0, 1.0, 0.1).unwrap();

        let analytical = (-lambda).exp();
        assert_relative_eq!(y_final[0], analytical, epsilon = 1e-3);
    }
}
