//! Exponential integrators for stiff equation systems.
//!
//! The state and matrix boundary in this module is Leto-owned:
//! [`TimeState`] is the rank-1 ODE state vector and [`TimeMatrix`] is the
//! rank-2 linear operator. Dense matrix exponentials are delegated to
//! `leto-ops` so the matrix-function implementation has one provider owner.
//!
//! Evidence tier for the implementation: analytical ETD/RK formulae plus
//! value-semantic tests against closed-form scalar ODE solutions. This is not
//! a machine-checked proof of the Padé approximation inside `leto-ops`.

use super::traits::{
    add_scaled_in_place, from_f64, state_len, state_norm, state_zeros, zero, TimeMatrix, TimeState,
};
use crate::error::Result;
use cfd_core::error::{ConvergenceErrorKind, Error};
use eunomia::RealField;
use leto_ops::{matexp, RealScalar};

/// Configuration for exponential integrators.
#[derive(Debug, Clone)]
pub struct ExponentialConfig {
    /// Krylov subspace dimension reserved for future sparse exponential paths.
    pub krylov_dimension: usize,
    /// Tolerance reserved for future Krylov convergence checks.
    pub krylov_tolerance: f64,
    /// Maximum iterations for fixed-point and series convergence.
    pub krylov_max_iter: usize,
    /// Tolerance for phi-function series terms.
    pub phi_tolerance: f64,
    /// Whether future sparse Krylov paths may adapt the subspace dimension.
    pub adaptive_krylov: bool,
}

impl Default for ExponentialConfig {
    fn default() -> Self {
        Self {
            krylov_dimension: 50,
            krylov_tolerance: 1e-12,
            krylov_max_iter: 100,
            phi_tolerance: 1e-12,
            adaptive_krylov: true,
        }
    }
}

/// Exponential Time Differencing (ETD) schemes.
pub struct ExponentialTimeDifferencing<T: RealField + RealScalar + Copy> {
    config: ExponentialConfig,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + RealScalar + Copy> Default for ExponentialTimeDifferencing<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + RealScalar + Copy> ExponentialTimeDifferencing<T> {
    /// Create a new ETD integrator.
    pub fn new() -> Self {
        Self::with_config(ExponentialConfig::default())
    }

    /// Create an ETD integrator with custom configuration.
    pub fn with_config(config: ExponentialConfig) -> Self {
        Self {
            config,
            _phantom: std::marker::PhantomData,
        }
    }

    /// ETD1 scheme for `du/dt = A*u + f(u)`.
    ///
    /// `u_{n+1} = exp(h A) u_n + h phi_1(h A) f(u_n)`.
    pub fn etd1(
        &self,
        u: &TimeState<T>,
        matrix: &TimeMatrix<T>,
        rhs: &TimeState<T>,
        dt: T,
    ) -> Result<TimeState<T>> {
        validate_matrix_vector(matrix, u, "ETD1 linear state")?;
        validate_matrix_vector(matrix, rhs, "ETD1 RHS")?;

        let mut u_new = self.matrix_exponential_vector(matrix, u, dt)?;
        let phi1_f = self.phi_vector(matrix, rhs, dt, PhiFunction::Phi1)?;
        add_scaled_in_place(&mut u_new, &phi1_f, dt, "ETD1 phi accumulation")?;
        Ok(u_new)
    }

    /// ETD2 scheme for `du/dt = A*u + f(u)`.
    ///
    /// `u_{n+1} = exp(h A) u_n + h phi_2(h A) (f(u_{n+1}) - f(u_n))`.
    pub fn etd2<F>(
        &self,
        u: &TimeState<T>,
        matrix: &TimeMatrix<T>,
        rhs_function: F,
        dt: T,
        tolerance: T,
        max_iter: usize,
    ) -> Result<TimeState<T>>
    where
        F: Fn(&TimeState<T>) -> TimeState<T>,
    {
        validate_matrix_vector(matrix, u, "ETD2 linear state")?;

        let exp_au = self.matrix_exponential_vector(matrix, u, dt)?;
        let f_old = rhs_function(u);
        validate_state_len(&f_old, state_len(u), "ETD2 old RHS")?;

        let mut u_new = u.clone();
        for _ in 0..max_iter {
            let f_new = rhs_function(&u_new);
            validate_state_len(&f_new, state_len(u), "ETD2 new RHS")?;

            let f_diff = state_sub(&f_new, &f_old, "ETD2 RHS difference")?;
            let phi2_diff = self.phi_vector(matrix, &f_diff, dt, PhiFunction::Phi2)?;
            let u_updated = state_add_scaled(&exp_au, &phi2_diff, dt, "ETD2 update")?;
            let residual = state_norm(&state_sub(&u_updated, &u_new, "ETD2 residual")?);
            u_new = u_updated;

            if residual < tolerance {
                return Ok(u_new);
            }
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded { max: max_iter },
        ))
    }
}

/// Exponential Runge-Kutta methods.
pub struct ExponentialRungeKutta4<T: RealField + RealScalar + Copy> {
    config: ExponentialConfig,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + RealScalar + Copy> Default for ExponentialRungeKutta4<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + RealScalar + Copy> ExponentialRungeKutta4<T> {
    /// Create a new exponential RK4 integrator.
    pub fn new() -> Self {
        Self::with_config(ExponentialConfig::default())
    }

    /// Create an ERK4 integrator with custom configuration.
    pub fn with_config(config: ExponentialConfig) -> Self {
        Self {
            config,
            _phantom: std::marker::PhantomData,
        }
    }

    /// Perform one fourth-order explicit RK step over Leto state storage.
    pub fn step<F>(&self, u: &TimeState<T>, rhs_function: F, dt: T) -> Result<TimeState<T>>
    where
        F: Fn(&TimeState<T>) -> TimeState<T>,
    {
        let half = from_f64::<T>(0.5);
        let two = from_f64::<T>(2.0);
        let six = from_f64::<T>(6.0);
        let n = state_len(u);

        let f0 = rhs_function(u);
        validate_state_len(&f0, n, "ERK4 first RHS")?;

        let u1 = state_add_scaled(u, &f0, dt * half, "ERK4 first stage")?;
        let f1 = rhs_function(&u1);
        validate_state_len(&f1, n, "ERK4 second RHS")?;

        let u2 = state_add_scaled(u, &f1, dt * half, "ERK4 second stage")?;
        let f2 = rhs_function(&u2);
        validate_state_len(&f2, n, "ERK4 third RHS")?;

        let u3 = state_add_scaled(u, &f2, dt, "ERK4 third stage")?;
        let f3 = rhs_function(&u3);
        validate_state_len(&f3, n, "ERK4 fourth RHS")?;

        let mut u_new = state_zeros(n);
        for i in 0..n {
            u_new[i] = u[i] + (dt / six) * (f0[i] + two * f1[i] + two * f2[i] + f3[i]);
        }

        Ok(u_new)
    }

    /// Return the stored configuration.
    pub fn config(&self) -> &ExponentialConfig {
        &self.config
    }
}

impl<T: RealField + RealScalar + Copy> ExponentialTimeDifferencing<T> {
    /// Compute `exp(scale * A) * v` using the Leto-ops matrix exponential.
    fn matrix_exponential_vector(
        &self,
        matrix: &TimeMatrix<T>,
        vector: &TimeState<T>,
        scale: T,
    ) -> Result<TimeState<T>> {
        validate_matrix_vector(matrix, vector, "matrix exponential vector")?;
        let scaled_matrix = matrix_scale(matrix, scale);
        let exp_matrix = matexp(&scaled_matrix.view()).map_err(|error| {
            Error::InvalidInput(format!("Leto matrix exponential failed: {error}"))
        })?;
        matrix_vector(&exp_matrix, vector, "matrix exponential product")
    }

    /// Compute `phi_k(scale * A) * v` through its Taylor definition.
    fn phi_vector(
        &self,
        matrix: &TimeMatrix<T>,
        vector: &TimeState<T>,
        scale: T,
        function: PhiFunction,
    ) -> Result<TimeState<T>> {
        validate_matrix_vector(matrix, vector, "phi-function vector")?;
        let scaled_matrix = matrix_scale(matrix, scale);
        let n = state_len(vector);
        let mut result = match function {
            PhiFunction::Phi1 => vector.clone(),
            PhiFunction::Phi2 => state_scale(vector, from_f64(0.5)),
        };
        let mut term = result.clone();
        let tolerance = from_f64::<T>(self.config.phi_tolerance);

        for iteration in 0..self.config.krylov_max_iter {
            let denominator = function.next_denominator::<T>(iteration);
            term = matrix_vector(&scaled_matrix, &term, "phi-function recurrence")?;
            for i in 0..n {
                term[i] = term[i] / denominator;
                result[i] += term[i];
            }
            if state_norm(&term) < tolerance {
                return Ok(result);
            }
        }

        Err(Error::Convergence(
            ConvergenceErrorKind::MaxIterationsExceeded {
                max: self.config.krylov_max_iter,
            },
        ))
    }
}

#[derive(Clone, Copy)]
enum PhiFunction {
    Phi1,
    Phi2,
}

impl PhiFunction {
    fn next_denominator<T: RealScalar>(self, iteration: usize) -> T {
        match self {
            Self::Phi1 => from_f64((iteration + 2) as f64),
            Self::Phi2 => from_f64((iteration + 3) as f64),
        }
    }
}

fn validate_state_len<T>(state: &TimeState<T>, expected: usize, context: &str) -> Result<()> {
    let actual = state_len(state);
    if actual != expected {
        return Err(Error::InvalidInput(format!(
            "{context}: state has {actual} entries, expected {expected}"
        )));
    }
    Ok(())
}

fn validate_matrix_vector<T>(
    matrix: &TimeMatrix<T>,
    vector: &TimeState<T>,
    context: &str,
) -> Result<usize> {
    let [rows, cols] = matrix.shape();
    let vector_len = state_len(vector);
    if rows != cols {
        return Err(Error::InvalidInput(format!(
            "{context}: matrix has shape {rows}x{cols}, expected square"
        )));
    }
    if cols != vector_len {
        return Err(Error::InvalidInput(format!(
            "{context}: matrix has {cols} columns, vector has {vector_len} entries"
        )));
    }
    Ok(rows)
}

fn matrix_scale<T: RealScalar>(matrix: &TimeMatrix<T>, scale: T) -> TimeMatrix<T> {
    let [rows, cols] = matrix.shape();
    TimeMatrix::from_shape_vec(
        [rows, cols],
        (0..rows * cols)
            .map(|idx| {
                let row = idx / cols;
                let col = idx % cols;
                matrix[[row, col]] * scale
            })
            .collect(),
    )
    .expect("invariant: scaled matrix storage matches shape")
}

fn matrix_vector<T: RealScalar>(
    matrix: &TimeMatrix<T>,
    vector: &TimeState<T>,
    context: &str,
) -> Result<TimeState<T>> {
    let rows = validate_matrix_vector(matrix, vector, context)?;
    let cols = state_len(vector);
    let mut output = state_zeros(rows);
    for row in 0..rows {
        let mut sum = zero::<T>();
        for col in 0..cols {
            sum += matrix[[row, col]] * vector[col];
        }
        output[row] = sum;
    }
    Ok(output)
}

fn state_add_scaled<T: RealField + RealScalar + Copy>(
    base: &TimeState<T>,
    increment: &TimeState<T>,
    scale: T,
    context: &str,
) -> Result<TimeState<T>> {
    let mut output = base.clone();
    add_scaled_in_place(&mut output, increment, scale, context)?;
    Ok(output)
}

fn state_sub<T: RealField + RealScalar + Copy>(
    lhs: &TimeState<T>,
    rhs: &TimeState<T>,
    context: &str,
) -> Result<TimeState<T>> {
    let n = state_len(lhs);
    validate_state_len(rhs, n, context)?;
    let mut output = state_zeros(n);
    for i in 0..n {
        output[i] = lhs[i] - rhs[i];
    }
    Ok(output)
}

fn state_scale<T: RealScalar>(state: &TimeState<T>, scale: T) -> TimeState<T> {
    let n = state_len(state);
    TimeState::from_shape_vec([n], (0..n).map(|i| state[i] * scale).collect())
        .expect("invariant: scaled state storage matches shape")
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::time_stepping::traits::{state_from_vec, state_neg};
    use approx::assert_relative_eq;

    fn matrix_from_row_slice(rows: usize, cols: usize, values: &[f64]) -> TimeMatrix<f64> {
        TimeMatrix::from_shape_vec([rows, cols], values.to_vec())
            .expect("invariant: test matrix storage matches shape")
    }

    #[test]
    fn test_etd1_simple() {
        let etd = ExponentialTimeDifferencing::<f64>::new();

        let u0 = state_from_vec(vec![1.0]);
        let matrix = matrix_from_row_slice(1, 1, &[-1.0]);
        let rhs = state_from_vec(vec![0.0]);
        let dt = 0.1;

        let result = etd.etd1(&u0, &matrix, &rhs, dt).unwrap();

        let analytical = (-dt).exp();
        assert_relative_eq!(result[0], analytical, epsilon = 1e-12);
    }

    #[test]
    fn test_etd1_phi_source_term() {
        let etd = ExponentialTimeDifferencing::<f64>::new();

        let u0 = state_from_vec(vec![0.0]);
        let matrix = matrix_from_row_slice(1, 1, &[-2.0]);
        let rhs = state_from_vec(vec![3.0]);
        let dt = 0.25;

        let result = etd.etd1(&u0, &matrix, &rhs, dt).unwrap();

        let analytical = 3.0 * (1.0 - (-2.0 * dt).exp()) / 2.0;
        assert_relative_eq!(result[0], analytical, epsilon = 1e-12);
    }

    #[test]
    fn test_erk4_simple() {
        let erk4 = ExponentialRungeKutta4::<f64>::new();

        let rhs = |u: &TimeState<f64>| state_neg(u);
        let u0 = state_from_vec(vec![1.0]);
        let dt = 0.1;

        let result = erk4.step(&u0, rhs, dt).unwrap();

        let analytical = (-dt).exp();
        assert_relative_eq!(result[0], analytical, epsilon = 1e-6);
    }

    #[test]
    fn etd_rejects_dimension_mismatch() {
        let etd = ExponentialTimeDifferencing::<f64>::new();
        let u0 = state_from_vec(vec![1.0, 2.0]);
        let matrix = matrix_from_row_slice(1, 1, &[-1.0]);
        let rhs = state_from_vec(vec![0.0]);

        let result = etd.etd1(&u0, &matrix, &rhs, 0.1);

        assert!(matches!(result, Err(Error::InvalidInput(_))));
    }
}
