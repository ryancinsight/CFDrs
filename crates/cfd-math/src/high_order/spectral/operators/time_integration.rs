//! Spectral time integration methods.
//!
//! Split out of `operators.rs` so the differentiation, interpolation,
//! quadrature and filtering operators keep one file and the time integrators
//! keep another; the public path (`operators::time_integration`) is unchanged.

use leto::Array1;

fn add_scaled(lhs: &Array1<f64>, rhs: &Array1<f64>, scale: f64) -> Array1<f64> {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: add_scaled operands must have equal length"
    );

    Array1::from_shape_fn(lhs.shape(), |[i]| lhs[i] + scale * rhs[i])
}

fn add_assign_scaled(lhs: &mut Array1<f64>, rhs: &Array1<f64>, scale: f64) {
    assert_eq!(
        lhs.shape(),
        rhs.shape(),
        "invariant: add_assign_scaled operands must have equal length"
    );

    for i in 0..lhs.shape()[0] {
        lhs[i] += scale * rhs[i];
    }
}

/// Runge-Kutta time integration for spectral methods
pub struct RungeKutta {
    /// Butcher tableau, row-major over `stages × stages`: coefficient `(i, j)`
    /// is `a[i * stages + j]`.
    a: Vec<f64>,
    /// Time step coefficients
    b: Vec<f64>,
    /// Node coefficients
    c: Vec<f64>,
    /// Number of stages
    stages: usize,
}

impl RungeKutta {
    /// Create a new Runge-Kutta method from Butcher tableau
    ///
    /// Rows shorter than the stage count are zero-padded and longer rows are
    /// truncated, exactly as the per-row storage did: only the strict lower
    /// triangle `j < i < stages` is ever read.
    pub fn new(a: Vec<Vec<f64>>, b: Vec<f64>, c: Vec<f64>) -> Self {
        let stages = b.len();

        // One row-major block instead of a `Vec` per stage: the stepping loop
        // reads each stage's coefficients contiguously.
        let mut tableau = vec![0.0; stages * stages];
        for stage in 0..stages {
            let target = &mut tableau[stage * stages..(stage + 1) * stages];
            for (slot, &coefficient) in target.iter_mut().zip(&a[stage]) {
                *slot = coefficient;
            }
        }

        Self {
            a: tableau,
            b,
            c,
            stages,
        }
    }

    /// Create the classical 4th-order Runge-Kutta method
    pub fn rk4() -> Self {
        let a = vec![vec![], vec![0.5], vec![0.0, 0.5], vec![0.0, 0.0, 1.0]];

        let b = vec![1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0];
        let c = vec![0.0, 0.5, 0.5, 1.0];

        Self::new(a, b, c)
    }

    /// Perform one time step
    pub fn step<F>(&self, f: F, t: f64, dt: f64, u: &Array1<f64>) -> Array1<f64>
    where
        F: Fn(f64, &Array1<f64>) -> Array1<f64>,
    {
        let mut k = Vec::with_capacity(self.stages);
        let mut u_new = u.clone();

        // Compute stage values
        for i in 0..self.stages {
            let mut sum = Array1::zeros(u.shape());

            for j in 0..i {
                add_assign_scaled(&mut sum, &k[j], self.a[i * self.stages + j]);
            }

            let t_stage = t + self.c[i] * dt;
            let u_stage = add_scaled(u, &sum, dt);

            k.push(f(t_stage, &u_stage));
        }

        // Combine stages
        for i in 0..self.stages {
            add_assign_scaled(&mut u_new, &k[i], dt * self.b[i]);
        }

        u_new
    }
}

/// Adams-Bashforth time integration for spectral methods
pub struct AdamsBashforth {
    /// Coefficients
    coeffs: Vec<f64>,
    /// Order of the method
    order: usize,
}

impl AdamsBashforth {
    /// Create a new Adams-Bashforth method of given order
    pub fn new(order: usize) -> Self {
        // Coefficients for Adams-Bashforth methods
        let coeffs = match order {
            1 => vec![1.0],
            2 => vec![3.0 / 2.0, -1.0 / 2.0],
            3 => vec![23.0 / 12.0, -16.0 / 12.0, 5.0 / 12.0],
            4 => vec![55.0 / 24.0, -59.0 / 24.0, 37.0 / 24.0, -9.0 / 24.0],
            _ => panic!("Unsupported Adams-Bashforth order"),
        };

        Self { coeffs, order }
    }

    /// Perform one time step
    pub fn step<F>(
        &self,
        f: F,
        t: f64,
        dt: f64,
        u: &Array1<f64>,
        f_prev: &[Array1<f64>],
    ) -> (Array1<f64>, Vec<Array1<f64>>)
    where
        F: Fn(f64, &Array1<f64>) -> Array1<f64>,
    {
        let f_new = f(t, u);
        let mut sum = Array1::zeros(u.shape());
        add_assign_scaled(&mut sum, &f_new, self.coeffs[0]);

        for (i, f_prev_i) in f_prev.iter().take(self.order - 1).enumerate() {
            add_assign_scaled(&mut sum, f_prev_i, self.coeffs[i + 1]);
        }

        let u_new = add_scaled(u, &sum, dt);

        // Update history of derivatives
        let mut f_history = Vec::with_capacity(self.order);
        f_history.push(f_new);
        f_history.extend_from_slice(&f_prev[..self.order - 1]);

        (u_new, f_history)
    }
}

#[cfg(test)]
mod tests {
    use super::{RungeKutta, add_assign_scaled, add_scaled};
    use leto::Array1;

    /// The pre-flattening step, kept verbatim as the oracle: one `Vec<f64>` per
    /// stage, indexed `a[i][j]`.
    fn step_per_row(
        a: &[Vec<f64>],
        b: &[f64],
        c: &[f64],
        f: &dyn Fn(f64, &Array1<f64>) -> Array1<f64>,
        t: f64,
        dt: f64,
        u: &Array1<f64>,
    ) -> Array1<f64> {
        let stages = b.len();
        let mut k: Vec<Array1<f64>> = Vec::with_capacity(stages);
        let mut u_new = u.clone();

        for i in 0..stages {
            let mut sum = Array1::zeros(u.shape());
            for j in 0..i {
                add_assign_scaled(&mut sum, &k[j], a[i][j]);
            }
            let u_stage = add_scaled(u, &sum, dt);
            k.push(f(t + c[i] * dt, &u_stage));
        }

        for i in 0..stages {
            add_assign_scaled(&mut u_new, &k[i], dt * b[i]);
        }

        u_new
    }

    fn rk4_tableau() -> (Vec<Vec<f64>>, Vec<f64>, Vec<f64>) {
        (
            vec![vec![], vec![0.5], vec![0.0, 0.5], vec![0.0, 0.0, 1.0]],
            vec![1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0],
            vec![0.0, 0.5, 0.5, 1.0],
        )
    }

    #[test]
    fn the_tableau_is_one_row_major_block() {
        let rk = RungeKutta::rk4();
        assert_eq!(rk.a.len(), rk.stages * rk.stages);
        // First row is empty in the source tableau, so it is padded with zeros.
        assert_eq!(rk.a[0], 0.0);
        assert_eq!(rk.a[rk.stages], 0.5);
        assert_eq!(rk.a[2 * rk.stages + 1], 0.5);
        assert_eq!(rk.a[3 * rk.stages + 2], 1.0);
    }

    #[test]
    fn flattened_tableau_matches_the_per_row_step_bitwise() {
        let (a, b, c) = rk4_tableau();
        let rk = RungeKutta::new(a.clone(), b.clone(), c.clone());
        const VALUES: [f64; 5] = [1.0, 1.37, 1.74, 2.11, 2.48];
        let u = Array1::from_shape_fn([5], |[i]| VALUES[i]);
        let f = |_t: f64, v: &Array1<f64>| Array1::from_shape_fn(v.shape(), |[i]| 0.1 - 0.9 * v[i]);

        for step in 0..4_i32 {
            let t = f64::from(step) * 0.25;
            let dt = 0.1;
            let flat = rk.step(f, t, dt, &u);
            let per_row = step_per_row(&a, &b, &c, &f, t, dt, &u);
            assert_eq!(
                flat.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                per_row.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                "stage {step} diverged from the per-row tableau"
            );
        }
    }
}
