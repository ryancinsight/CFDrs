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
use crate::time_stepping::traits::{from_f64, one, state_len, state_zeros, zero, TimeState};
use eunomia::RealField;
use eunomia::{FloatElement, NumericElement};

fn from_usize<T: FloatElement>(value: usize) -> T {
    const MAX_EXACT_F64_INTEGER: u64 = 1_u64 << f64::MANTISSA_DIGITS;
    let value_u64 =
        u64::try_from(value).expect("invariant: usize value fits in u64 for RKC scalar conversion");
    assert!(
        value_u64 <= MAX_EXACT_F64_INTEGER,
        "RKC scalar conversion requires exact f64 representation of usize value"
    );
    from_f64(<u64 as NumericElement>::to_f64(value_u64))
}

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

impl<T: RealField + Copy + FloatElement> Default for RkcConfig<T> {
    fn default() -> Self {
        Self {
            num_stages: 10,
            damping: from_f64(2.0 / 13.0),
            atol: from_f64(1e-8),
            rtol: from_f64(1e-6),
            safety_factor: from_f64(0.9),
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
    fn evaluate(&self, t: T, y: &TimeState<T>) -> Result<TimeState<T>>;
}

impl<T: RealField + Copy + FloatElement> RungeKuttaChebyshev<T> {
    /// Create RKC method with default configuration
    pub fn new() -> Self {
        Self::with_config(RkcConfig::default())
    }

    /// Create RKC method with custom configuration
    pub fn with_config(config: RkcConfig<T>) -> Self {
        let coeffs = Self::compute_coefficients(config.num_stages, config.damping);
        Self { config, coeffs }
    }

    /// Precompute RKC coefficients for given number of stages and damping
    fn compute_coefficients(s: usize, epsilon: T) -> Vec<RkcStageCoeffs<T>> {
        let mut coeffs = Vec::with_capacity(s + 1);

        // Stage 0 and 1 are special or unused in the recurrence loop (j=2..s)
        // We define dummy values for index 0 and 1 to align with 1-based indexing logic
        coeffs.push(RkcStageCoeffs {
            mu: zero::<T>(),
            nu: zero::<T>(),
            mu_tilde: zero::<T>(),
            gamma_tilde: zero::<T>(),
            c: zero::<T>(),
        }); // Index 0

        // w0 = 1 + epsilon / s^2
        let two: T = from_f64(2.0);
        let four: T = from_f64(4.0);
        let s_t = from_usize(s);
        let s_sq = s_t * s_t;
        let w0 = one::<T>() + epsilon / s_sq;

        let temp1 = w0 * w0 - one::<T>();
        let temp2 = <T as NumericElement>::sqrt(temp1);
        let arg: T = s_t * <T as FloatElement>::ln(w0 + temp2);
        let sinh_arg = <T as FloatElement>::sinh(arg);
        let cosh_arg = <T as FloatElement>::cosh(arg);
        let w1 = sinh_arg * temp1 / (cosh_arg * s_t * temp2 - w0 * sinh_arg);

        let mut bjm1 = one::<T>() / ((two * w0) * (two * w0));
        let mut bjm2 = bjm1;

        let mu_tilde_1 = w1 * bjm1;
        coeffs.push(RkcStageCoeffs {
            mu: zero::<T>(),
            nu: zero::<T>(),
            mu_tilde: mu_tilde_1,
            gamma_tilde: zero::<T>(),
            c: mu_tilde_1,
        }); // Index 1

        let mut zjm1 = w0;
        let mut zjm2 = one::<T>();
        let mut dzjm1 = one::<T>();
        let mut dzjm2 = zero::<T>();
        let mut d2zjm1 = zero::<T>();
        let mut d2zjm2 = zero::<T>();

        for j in 2..=s {
            let zj: T = two * w0 * zjm1 - zjm2;
            let dzj: T = two * w0 * dzjm1 - dzjm2 + two * zjm1;
            let d2zj: T = two * w0 * d2zjm1 - d2zjm2 + four * dzjm1;

            let bj: T = d2zj / (dzj * dzj);
            let mu_j: T = two * w0 * bj / bjm1;
            let nu_j: T = -bj / bjm2;
            let mu_tilde_j: T = mu_j * w1 / w0;

            let ajm1: T = one::<T>() - zjm1 * bjm1;
            let gamma_tilde_j: T = -ajm1 * mu_tilde_j;

            let c_prev = coeffs[j - 1].c;
            let c_prev2 = coeffs[j - 2].c;
            let c_j: T = mu_j * c_prev + nu_j * c_prev2 + mu_tilde_j + gamma_tilde_j;

            coeffs.push(RkcStageCoeffs {
                mu: mu_j,
                nu: nu_j,
                mu_tilde: mu_tilde_j,
                gamma_tilde: gamma_tilde_j,
                c: c_j,
            });

            zjm2 = zjm1;
            zjm1 = zj;
            dzjm2 = dzjm1;
            dzjm1 = dzj;
            d2zjm2 = d2zjm1;
            d2zjm1 = d2zj;
            bjm2 = bjm1;
            bjm1 = bj;
        }

        coeffs
    }

    fn step_raw<F: RhsFunction<T>>(
        &self,
        rhs: &F,
        t0: T,
        y0: &TimeState<T>,
        dt: T,
    ) -> Result<TimeState<T>> {
        let n = state_len(y0);

        // Y_0 = y_n
        let mut y_prev2 = y0.clone(); // Y_{j-2} (initially Y_0)

        // Evaluate F_0 = F(t0, Y_0)
        let f0 = rhs.evaluate(t0, &y_prev2)?;

        // Y_1 = Y_0 + mu_tilde_1 * dt * F_0
        let mu_tilde_1 = self.coeffs[1].mu_tilde;
        let mut y_prev1 = state_zeros(n); // Y_{j-1} (initially Y_1)
        for i in 0..n {
            y_prev1[i] = y_prev2[i] + mu_tilde_1 * dt * f0[i];
        }

        // Stages 2 to s
        for j in 2..=self.config.num_stages {
            let coeff = &self.coeffs[j];

            let t_stage = t0 + self.coeffs[j - 1].c * dt;
            let f_prev = rhs.evaluate(t_stage, &y_prev1)?;

            // Y_j = (1 - mu - nu) Y_0 + mu Y_{j-1} + nu Y_{j-2} + mu_tilde dt F_{j-1} + gamma_tilde dt F_0
            let term_y0 = one::<T>() - coeff.mu - coeff.nu;

            let mut y_curr = state_zeros(n);
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

    fn integrate_adaptive<F: RhsFunction<T>>(
        &self,
        rhs: &F,
        t0: T,
        y0: &TimeState<T>,
        t_final: T,
        dt_initial: T,
    ) -> Result<(TimeState<T>, T)> {
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
                        cfd_core::error::ConvergenceErrorKind::MaxIterationsExceeded {
                            max: max_attempts,
                        },
                    ));
                }

                let y_full = self.step_raw(rhs, t, &y, dt)?;
                let dt_half = dt / from_f64(2.0);
                let y_half = self.step_raw(rhs, t, &y, dt_half)?;
                let y_half_full = self.step_raw(rhs, t + dt_half, &y_half, dt_half)?;

                let n = state_len(&y);
                let mut err_sum = zero::<T>();
                for i in 0..n {
                    let max_component = <T as NumericElement>::abs(y[i])
                        .max_scalar(<T as NumericElement>::abs(y_half_full[i]));
                    let scale = self.config.atol + self.config.rtol * max_component;
                    let denom = if scale > zero::<T>() {
                        scale
                    } else {
                        one::<T>()
                    };
                    let diff = <T as NumericElement>::abs(y_half_full[i] - y_full[i])
                        / (from_f64::<T>(3.0) * denom);
                    err_sum += diff * diff;
                }
                let n_t = from_usize(n);
                let err_norm = <T as NumericElement>::sqrt(err_sum / n_t);

                if err_norm <= one::<T>() {
                    y = y_half_full;
                    t += dt;
                    let factor = if err_norm == zero::<T>() {
                        from_f64(5.0)
                    } else {
                        let exponent: T = from_f64(1.0 / 3.0);
                        self.config.safety_factor * <T as FloatElement>::powf(err_norm, -exponent)
                    };
                    let factor = factor.max_scalar(from_f64(0.1)).min_scalar(from_f64(5.0));
                    dt *= factor;
                    accepted = true;
                } else {
                    let exponent: T = from_f64(1.0 / 3.0);
                    let factor =
                        self.config.safety_factor * <T as FloatElement>::powf(err_norm, -exponent);
                    let factor = factor.max_scalar(from_f64(0.1)).min_scalar(from_f64(0.5));
                    dt *= factor;
                    if t + dt > t_final {
                        dt = t_final - t;
                    }
                }
            }
        }

        Ok((y, dt))
    }

    /// Solve ODE system using RKC method
    pub fn step<F: RhsFunction<T>>(
        &self,
        rhs: &F,
        t0: T,
        y0: &TimeState<T>,
        dt: T,
    ) -> Result<TimeState<T>> {
        let t_final = t0 + dt;
        let (y, _) = self.integrate_adaptive(rhs, t0, y0, t_final, dt)?;
        Ok(y)
    }

    /// Solve ODE system with adaptive time stepping
    pub fn solve_adaptive<F: RhsFunction<T>>(
        &self,
        rhs: &F,
        t0: T,
        y0: &TimeState<T>,
        t_final: T,
        dt_initial: T,
    ) -> Result<(TimeState<T>, T)> {
        self.integrate_adaptive(rhs, t0, y0, t_final, dt_initial)
    }

    /// Get the configuration of the RKC solver
    pub fn config(&self) -> &RkcConfig<T> {
        &self.config
    }
}

impl<T: RealField + Copy + FloatElement> Default for RungeKuttaChebyshev<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time_stepping::traits::{state_from_vec, state_scale};
    use approx::assert_relative_eq;

    struct ExponentialDecay<T: RealField + Copy> {
        lambda: T,
    }

    impl<T: RealField + Copy> ExponentialDecay<T> {
        fn new(lambda: T) -> Self {
            Self { lambda }
        }
    }

    impl<T: RealField + Copy> RhsFunction<T> for ExponentialDecay<T> {
        fn evaluate(&self, _t: T, y: &TimeState<T>) -> Result<TimeState<T>> {
            Ok(state_scale(y, -self.lambda))
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

        let y0 = state_from_vec(vec![1.0]);
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

        let y0 = state_from_vec(vec![1.0]);
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

        let y0 = state_from_vec(vec![1.0]);
        let (y_final, _) = rkc.solve_adaptive(&rhs, 0.0, &y0, 1.0, 0.1).unwrap();

        let analytical = (-lambda).exp();
        assert_relative_eq!(y_final[0], analytical, epsilon = 1e-3);
    }
}
