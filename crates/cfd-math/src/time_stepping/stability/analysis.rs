//! Runge-Kutta stability region computation and von Neumann analysis.
//!
//! Provides stability function evaluation for explicit RK methods via
//! $R(z) = 1 + z\,b^T (I - z\,A)^{-1} \mathbf{1}$ (Hairer et al. 1993, §IV.2).

use super::{
    ComplexPoint, MethodInfo, NumericalScheme, StabilityAnalyzer, StabilityRegion, StabilityType,
    VonNeumannAnalysis,
};
use cfd_core::error::{Error, Result};
use nalgebra::{DMatrix, DVector, RealField};
use num_complex::Complex as NumComplex;
use num_traits::ToPrimitive;
use std::f64::consts::PI;

impl<T: RealField + Copy + ToPrimitive> StabilityAnalyzer<T> {
    /// Compute stability region for explicit Runge-Kutta method
    ///
    /// For an s-stage explicit RK method with Butcher tableau (A, b, c),
    /// the stability function is R(z) = 1 + z*b^T * (I - z*A)^(-1) * 1
    ///
    /// # Arguments
    /// * `a` - A matrix from Butcher tableau
    /// * `b` - b vector from Butcher tableau
    /// * `c` - c vector from Butcher tableau
    ///
    /// # Returns
    /// Stability region boundary as vector of complex numbers
    pub fn compute_rk_stability_region(
        &self,
        a: &DMatrix<T>,
        b: &DVector<T>,
        c: &DVector<T>,
    ) -> Result<StabilityRegion<T>> {
        let s = b.len();

        // Validate Butcher tableau dimensions
        if a.nrows() != s || a.ncols() != s || c.len() != s {
            return Err(Error::InvalidInput(format!(
                "Invalid Butcher tableau dimensions: A is {}x{}, b has {}, c has {}",
                a.nrows(),
                a.ncols(),
                b.len(),
                c.len()
            )));
        }

        // Check if method is explicit (lower triangular A with zero diagonal)
        for i in 0..s {
            if a[(i, i)] != T::zero() {
                return Err(Error::InvalidInput(
                    "Not an explicit Runge-Kutta method".to_string(),
                ));
            }
            for j in (i + 1)..s {
                if a[(i, j)] != T::zero() {
                    return Err(Error::InvalidInput(
                        "Not a lower triangular A matrix".to_string(),
                    ));
                }
            }
        }

        let max_r = self
            .max_z
            .to_f64()
            .ok_or_else(|| Error::InvalidInput("max_z must be convertible to f64".to_string()))?;

        let mut boundary_points = Vec::with_capacity(self.resolution);
        let mut interior_points = Vec::with_capacity(self.resolution * 3);

        // Approximate the stability boundary by searching along rays for each angle θ.
        // For each θ, find the largest radius r such that |R(r e^{iθ})| <= 1.
        for i in 0..self.resolution {
            let theta = 2.0 * PI * i as f64 / self.resolution as f64;
            let r_boundary = self.find_rk_boundary_radius(a, b, c, theta, max_r)?;

            let z = NumComplex::new(r_boundary * theta.cos(), r_boundary * theta.sin());
            let z_real = T::from_f64(z.re).ok_or_else(|| {
                Error::InvalidInput("stability boundary real part not representable in T".to_string())
            })?;
            let z_imag = T::from_f64(z.im).ok_or_else(|| {
                Error::InvalidInput("stability boundary imag part not representable in T".to_string())
            })?;

            boundary_points.push(ComplexPoint {
                real: z_real,
                imag: z_imag,
                stability: true,
            });

            // Add a few interior samples along the same ray (useful for visualization).
            for frac in [0.25_f64, 0.5_f64, 0.75_f64] {
                let r = r_boundary * frac;
                let zi = NumComplex::new(r * theta.cos(), r * theta.sin());
                let zi_real = T::from_f64(zi.re).ok_or_else(|| {
                    Error::InvalidInput(
                        "stability interior sample real part not representable in T".to_string(),
                    )
                })?;
                let zi_imag = T::from_f64(zi.im).ok_or_else(|| {
                    Error::InvalidInput(
                        "stability interior sample imag part not representable in T".to_string(),
                    )
                })?;
                interior_points.push(ComplexPoint {
                    real: zi_real,
                    imag: zi_imag,
                    stability: true,
                });
            }
        }

        Ok(StabilityRegion {
            boundary: boundary_points,
            interior: interior_points,
            method_info: MethodInfo {
                name: "Explicit Runge-Kutta".to_string(),
                order: self.estimate_rk_order(a, b, c),
                stages: s,
                stability_type: StabilityType::Explicit,
            },
        })
    }

    /// Compute the absolute stability limit on the negative real axis for an explicit RK method.
    ///
    /// Returns the largest $r \ge 0$ such that $|R(-r)| \le 1$.
    pub fn compute_rk_absolute_stability_limit(
        &self,
        a: &DMatrix<T>,
        b: &DVector<T>,
        c: &DVector<T>,
    ) -> Result<T> {
        let max_r = self
            .max_z
            .to_f64()
            .ok_or_else(|| Error::InvalidInput("max_z must be convertible to f64".to_string()))?;

        // θ = π corresponds to the negative real axis.
        let r = self.find_rk_boundary_radius(a, b, c, PI, max_r)?;

        T::from_f64(r).ok_or_else(|| {
            Error::InvalidInput("absolute stability limit not representable in T".to_string())
        })
    }

    fn find_rk_boundary_radius(
        &self,
        a: &DMatrix<T>,
        b: &DVector<T>,
        c: &DVector<T>,
        theta: f64,
        max_r: f64,
    ) -> Result<f64> {
        const STABLE_TOL: f64 = 1e-12;
        const RADIAL_SCAN_STEPS: usize = 120;
        const BISECTION_ITERS: usize = 50;

        let mut prev_r = 0.0;
        let mut prev_stable = true;

        let mut bracket_low = 0.0;
        let mut bracket_high: Option<f64> = None;

        // Coarse scan from 0 -> max_r looking for first instability.
        for step in 1..=RADIAL_SCAN_STEPS {
            let r = max_r * (step as f64) / (RADIAL_SCAN_STEPS as f64);
            let z = NumComplex::new(r * theta.cos(), r * theta.sin());
            let r_z = self.compute_rk_stability_function(a, b, c, z)?;
            let stable = r_z.norm() <= 1.0 + STABLE_TOL;

            if prev_stable && !stable {
                bracket_low = prev_r;
                bracket_high = Some(r);
                break;
            }

            prev_r = r;
            prev_stable = stable;
        }

        // If we never found an instability up to max_r, return max_r.
        let Some(mut high) = bracket_high else {
            return Ok(max_r);
        };
        let mut low = bracket_low;

        // Refine boundary via bisection.
        for _ in 0..BISECTION_ITERS {
            let mid = 0.5 * (low + high);
            let z = NumComplex::new(mid * theta.cos(), mid * theta.sin());
            let r_z = self.compute_rk_stability_function(a, b, c, z)?;
            let stable = r_z.norm() <= 1.0 + STABLE_TOL;
            if stable {
                low = mid;
            } else {
                high = mid;
            }
        }

        Ok(0.5 * (low + high))
    }

    /// Compute stability function for explicit Runge-Kutta method
    fn compute_rk_stability_function(
        &self,
        a: &DMatrix<T>,
        b: &DVector<T>,
        _c: &DVector<T>,
        z: NumComplex<f64>,
    ) -> Result<NumComplex<f64>> {
        let s = b.len();

        // For explicit methods, R(z) = 1 + z·b^T (I - z·A)^{-1} 𝟏
        //
        // This matrix form evaluates the stability polynomial *exactly*
        // for any s-stage explicit RK method.  For an order-p method the
        // polynomial satisfies R(z) = Σ_{k=0}^{p} z^k/k! + O(z^{p+1}),
        // matching the Taylor expansion of exp(z) through order p.
        //
        // Reference: Hairer, Nørsett & Wanner (1993) §IV.2, Theorem 2.1.

        // Convert matrices to complex
        let mut a_complex = DMatrix::<NumComplex<f64>>::zeros(s, s);
        let mut b_complex = DVector::<NumComplex<f64>>::zeros(s);

        for i in 0..s {
            b_complex[i] = NumComplex::new(b[i].to_f64().unwrap(), 0.0);
            for j in 0..s {
                a_complex[(i, j)] = NumComplex::new(a[(i, j)].to_f64().unwrap(), 0.0);
            }
        }

        // Compute (I - z*A)
        let identity = DMatrix::<NumComplex<f64>>::identity(s, s);
        let z_a = &a_complex * z; // Matrix * scalar, not scalar * matrix
        let matrix = &identity - &z_a;

        // Compute inverse
        match matrix.try_inverse() {
            Some(inv_matrix) => {
                // Compute b^T * inv_matrix * 1 (ones vector)
                let ones = DVector::<NumComplex<f64>>::from_element(s, NumComplex::new(1.0, 0.0));
                let temp = &inv_matrix * &ones;
                let coeff = b_complex.dot(&temp);
                let r_z = NumComplex::new(1.0, 0.0) + z * coeff;
                Ok(r_z)
            }
            None => {
                // Singular matrix - method is unstable for this z
                Ok(NumComplex::new(f64::INFINITY, f64::INFINITY))
            }
        }
    }

    /// Estimate order of Runge-Kutta method from Butcher tableau
    fn estimate_rk_order(&self, _a: &DMatrix<T>, _b: &DVector<T>, _c: &DVector<T>) -> usize {
        let _s = _b.len();

        // Basic order conditions checks
        // Order 1: sum b_i = 1
        let sum_b = _b.iter().fold(T::zero(), |acc, &x| acc + x);
        if (sum_b - T::one()).abs() > T::from_f64(1e-10).unwrap() {
            return 0;
        }

        // Order 2: sum b_i * c_i = 1/2
        let sum_b_c = _b
            .iter()
            .zip(_c.iter())
            .fold(T::zero(), |acc, (&bi, &ci)| acc + bi * ci);
        if (sum_b_c - T::from_f64(0.5).unwrap()).abs() > T::from_f64(1e-10).unwrap() {
            return 1;
        }

        // Order 3: sum b_i * c_i^2 = 1/3
        let sum_b_c2 = _b
            .iter()
            .zip(_c.iter())
            .fold(T::zero(), |acc, (&bi, &ci)| acc + bi * ci * ci);
        if (sum_b_c2 - T::from_f64(1.0 / 3.0).unwrap()).abs() > T::from_f64(1e-10).unwrap() {
            return 2;
        }

        // Order 4: More complex conditions involving A matrix
        // For now, assume order 4 if basic conditions pass
        4
    }

    /// Perform von Neumann stability analysis for linear PDEs
    ///
    /// For a PDE ∂u/∂t = L u, where L is a linear spatial operator,
    /// the von Neumann method assumes solutions of the form u_j^n = g^n * e^(i k x_j)
    /// and analyzes the amplification factor |g| <= 1 for stability.
    ///
    /// # Arguments
    /// * `spatial_operator` - Function representing the spatial discretization L
    /// * `dt` - Time step size
    /// * `wave_numbers` - Range of wave numbers k to test
    ///
    /// # Returns
    /// Von Neumann stability analysis result
    pub fn von_neumann_analysis<F>(
        &self,
        spatial_operator: F,
        dt: T,
        wave_numbers: &[T],
    ) -> Result<VonNeumannAnalysis<T>>
    where
        F: Fn(NumComplex<f64>) -> NumComplex<f64>, // L_hat(k) - spatial operator in frequency domain
    {
        let mut amplification_factors = Vec::with_capacity(wave_numbers.len());
        let mut max_amplification = T::zero();
        let mut critical_wave_number = T::zero();

        for &k in wave_numbers {
            let k_complex = NumComplex::new(0.0, k.to_f64().unwrap());

            // Compute spatial operator in frequency domain
            let l_hat = spatial_operator(k_complex);

            // Amplification factor for forward Euler: g = 1 + dt * L_hat(k)
            let dt_f64 = dt.to_f64().unwrap();
            let g = NumComplex::new(1.0, 0.0) + NumComplex::new(dt_f64, 0.0) * l_hat;

            let amplification = T::from_f64(g.norm()).unwrap();
            amplification_factors.push(amplification);

            if amplification > max_amplification {
                max_amplification = amplification;
                critical_wave_number = k;
            }
        }

        let is_stable = max_amplification <= T::from_f64(1.0001).unwrap(); // Allow small numerical errors

        Ok(VonNeumannAnalysis {
            wave_numbers: wave_numbers.to_vec(),
            amplification_factors,
            max_amplification,
            critical_wave_number,
            is_stable,
            stability_margin: T::from_f64(1.0).unwrap() - max_amplification,
        })
    }

    /// Perform von Neumann stability analysis using an explicit Runge-Kutta method.
    ///
    /// For a linear semi-discrete system $u_t = \hat L(k) u$, an explicit RK method has
    /// amplification factor $g(k) = R(\Delta t\,\hat L(k))$.
    pub fn von_neumann_analysis_explicit_rk<F>(
        &self,
        a: &DMatrix<T>,
        b: &DVector<T>,
        c: &DVector<T>,
        spatial_operator: F,
        dt: T,
        wave_numbers: &[T],
    ) -> Result<VonNeumannAnalysis<T>>
    where
        F: Fn(NumComplex<f64>) -> NumComplex<f64>,
    {
        let dt_f64 = dt
            .to_f64()
            .ok_or_else(|| Error::InvalidInput("dt must be convertible to f64".to_string()))?;

        let mut amplification_factors = Vec::with_capacity(wave_numbers.len());
        let mut max_amplification = T::zero();
        let mut critical_wave_number = T::zero();

        for &k in wave_numbers {
            let k_complex = NumComplex::new(0.0, k.to_f64().unwrap());
            let l_hat = spatial_operator(k_complex);
            let z = NumComplex::new(dt_f64, 0.0) * l_hat;

            let g = self.compute_rk_stability_function(a, b, c, z)?;
            let amplification = T::from_f64(g.norm()).unwrap();
            amplification_factors.push(amplification);

            if amplification > max_amplification {
                max_amplification = amplification;
                critical_wave_number = k;
            }
        }

        let is_stable = max_amplification <= T::from_f64(1.0001).unwrap();

        Ok(VonNeumannAnalysis {
            wave_numbers: wave_numbers.to_vec(),
            amplification_factors,
            max_amplification,
            critical_wave_number,
            is_stable,
            stability_margin: T::from_f64(1.0).unwrap() - max_amplification,
        })
    }

    /// Convenience wrapper for von Neumann analysis using a built-in scheme.
    ///
    /// `ForwardEuler` uses the classic $g = 1 + \Delta t\,\hat L(k)$.
    /// `RK3` and `RK4` use their explicit RK stability functions.
    pub fn von_neumann_analysis_with_scheme<F>(
        &self,
        scheme: NumericalScheme,
        spatial_operator: F,
        dt: T,
        wave_numbers: &[T],
    ) -> Result<VonNeumannAnalysis<T>>
    where
        F: Fn(NumComplex<f64>) -> NumComplex<f64>,
    {
        match scheme {
            NumericalScheme::ForwardEuler => self.von_neumann_analysis(spatial_operator, dt, wave_numbers),
            NumericalScheme::RK3 => {
                // Heun/Kutta 3rd-order as used in validation
                let a = DMatrix::from_row_slice(
                    3,
                    3,
                    &[
                        T::zero(),
                        T::zero(),
                        T::zero(),
                        T::from_f64(1.0 / 3.0).unwrap(),
                        T::zero(),
                        T::zero(),
                        T::zero(),
                        T::from_f64(2.0 / 3.0).unwrap(),
                        T::zero(),
                    ],
                );
                let b = DVector::from_vec(vec![
                    T::from_f64(0.25).unwrap(),
                    T::zero(),
                    T::from_f64(0.75).unwrap(),
                ]);
                let c = DVector::from_vec(vec![
                    T::zero(),
                    T::from_f64(1.0 / 3.0).unwrap(),
                    T::from_f64(2.0 / 3.0).unwrap(),
                ]);
                self.von_neumann_analysis_explicit_rk(&a, &b, &c, spatial_operator, dt, wave_numbers)
            }
            NumericalScheme::RK4 => {
                let one_half = T::from_f64(0.5).unwrap();
                let a = DMatrix::from_row_slice(
                    4,
                    4,
                    &[
                        T::zero(),
                        T::zero(),
                        T::zero(),
                        T::zero(),
                        one_half,
                        T::zero(),
                        T::zero(),
                        T::zero(),
                        T::zero(),
                        one_half,
                        T::zero(),
                        T::zero(),
                        T::zero(),
                        T::zero(),
                        T::one(),
                        T::zero(),
                    ],
                );
                let b = DVector::from_vec(vec![
                    T::from_f64(1.0 / 6.0).unwrap(),
                    T::from_f64(1.0 / 3.0).unwrap(),
                    T::from_f64(1.0 / 3.0).unwrap(),
                    T::from_f64(1.0 / 6.0).unwrap(),
                ]);
                let c = DVector::from_vec(vec![T::zero(), one_half, one_half, T::one()]);
                self.von_neumann_analysis_explicit_rk(&a, &b, &c, spatial_operator, dt, wave_numbers)
            }
            _ => Err(Error::InvalidInput(
                "von Neumann analysis is only implemented for explicit schemes".to_string(),
            )),
        }
    }
}
