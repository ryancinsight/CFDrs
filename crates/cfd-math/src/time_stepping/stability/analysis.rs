//! Runge-Kutta stability region computation and von Neumann analysis.
//!
//! Provides stability function evaluation for explicit RK methods via
//! $R(z) = 1 + z\,b^T (I - z\,A)^{-1} \mathbf{1}$ (Hairer et al. 1993, §IV.2).

use super::{
    abs, from_f64, matrix_from_row_slice, one, to_f64, vector_from_vec, zero, ComplexPoint,
    MethodInfo, NumericalScheme, StabilityAnalyzer, StabilityRegion, StabilityType,
    VonNeumannAnalysis,
};
use cfd_core::error::{Error, Result};
use eunomia::Complex as AtlasComplex;
use eunomia::{FloatElement, RealField};
use leto::{Array1, Array2};
use std::f64::consts::PI;

impl<T: RealField + Copy + FloatElement> StabilityAnalyzer<T> {
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
        a: &Array2<T>,
        b: &Array1<T>,
        c: &Array1<T>,
    ) -> Result<StabilityRegion<T>> {
        let s = b.shape()[0];

        // Validate Butcher tableau dimensions
        if a.shape()[0] != s || a.shape()[1] != s || c.shape()[0] != s {
            return Err(Error::InvalidInput(format!(
                "Invalid Butcher tableau dimensions: A is {}x{}, b has {}, c has {}",
                a.shape()[0],
                a.shape()[1],
                b.shape()[0],
                c.shape()[0]
            )));
        }

        // Check if method is explicit (lower triangular A with zero diagonal)
        for i in 0..s {
            if a[[i, i]] != zero::<T>() {
                return Err(Error::InvalidInput(
                    "Not an explicit Runge-Kutta method".to_string(),
                ));
            }
            for j in (i + 1)..s {
                if a[[i, j]] != zero::<T>() {
                    return Err(Error::InvalidInput(
                        "Not a lower triangular A matrix".to_string(),
                    ));
                }
            }
        }

        let max_r = to_f64(self.max_z);

        let mut boundary_points = Vec::with_capacity(self.resolution);
        let mut interior_points = Vec::with_capacity(self.resolution * 3);

        // Approximate the stability boundary by searching along rays for each angle θ.
        // For each θ, find the largest radius r such that |R(r e^{iθ})| <= 1.
        for i in 0..self.resolution {
            let theta = 2.0 * PI * i as f64 / self.resolution as f64;
            let r_boundary = self.find_rk_boundary_radius(a, b, c, theta, max_r)?;

            let z = AtlasComplex::new(r_boundary * theta.cos(), r_boundary * theta.sin());
            let z_real = from_f64(z.re);
            let z_imag = from_f64(z.im);

            boundary_points.push(ComplexPoint {
                real: z_real,
                imag: z_imag,
                stability: true,
            });

            // Add a few interior samples along the same ray (useful for visualization).
            for frac in [0.25_f64, 0.5_f64, 0.75_f64] {
                let r = r_boundary * frac;
                let zi = AtlasComplex::new(r * theta.cos(), r * theta.sin());
                let zi_real = from_f64(zi.re);
                let zi_imag = from_f64(zi.im);
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
        a: &Array2<T>,
        b: &Array1<T>,
        c: &Array1<T>,
    ) -> Result<T> {
        let max_r = to_f64(self.max_z);

        // θ = π corresponds to the negative real axis.
        let r = self.find_rk_boundary_radius(a, b, c, PI, max_r)?;

        Ok(from_f64(r))
    }

    fn find_rk_boundary_radius(
        &self,
        a: &Array2<T>,
        b: &Array1<T>,
        c: &Array1<T>,
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
            let z = AtlasComplex::new(r * theta.cos(), r * theta.sin());
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
            let z = AtlasComplex::new(mid * theta.cos(), mid * theta.sin());
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
        a: &Array2<T>,
        b: &Array1<T>,
        _c: &Array1<T>,
        z: AtlasComplex<f64>,
    ) -> Result<AtlasComplex<f64>> {
        let s = b.shape()[0];

        // For explicit methods, R(z) = 1 + z·b^T (I - z·A)^{-1} 𝟏
        //
        // Since explicit RK tableaux have strictly lower-triangular A, the
        // system (I - z·A)x = 𝟏 has unit diagonal and is solved exactly by
        // forward substitution. For an order-p method the polynomial satisfies
        // R(z) = Σ_{k=0}^{p} z^k/k! + O(z^{p+1}), matching the Taylor
        // expansion of exp(z) through order p.
        //
        // Reference: Hairer, Nørsett & Wanner (1993) §IV.2, Theorem 2.1.

        let one = AtlasComplex::new(1.0, 0.0);
        let mut stage_values = vec![one; s];
        for i in 0..s {
            let mut lower_sum = AtlasComplex::new(0.0, 0.0);
            for j in 0..i {
                let a_ij = to_f64(a[[i, j]]);
                lower_sum += stage_values[j] * a_ij;
            }
            stage_values[i] = one + z * lower_sum;
        }

        let mut coeff = AtlasComplex::new(0.0, 0.0);
        for i in 0..s {
            let b_i = to_f64(b[i]);
            coeff += stage_values[i] * b_i;
        }
        Ok(one + z * coeff)
    }

    /// Estimate order of Runge-Kutta method from Butcher tableau
    fn estimate_rk_order(&self, _a: &Array2<T>, _b: &Array1<T>, _c: &Array1<T>) -> usize {
        let _s = _b.shape()[0];

        // Basic order conditions checks
        // Order 1: sum b_i = 1
        let sum_b = _b.iter().fold(zero::<T>(), |acc, &x| acc + x);
        if abs(sum_b - one::<T>()) > from_f64(1e-10) {
            return 0;
        }

        // Order 2: sum b_i * c_i = 1/2
        let sum_b_c = _b
            .iter()
            .zip(_c.iter())
            .fold(zero::<T>(), |acc, (&bi, &ci)| acc + bi * ci);
        if abs(sum_b_c - from_f64(0.5)) > from_f64(1e-10) {
            return 1;
        }

        // Order 3: sum b_i * c_i^2 = 1/3
        let sum_b_c2 = _b
            .iter()
            .zip(_c.iter())
            .fold(zero::<T>(), |acc, (&bi, &ci)| acc + bi * ci * ci);
        if abs(sum_b_c2 - from_f64(1.0 / 3.0)) > from_f64(1e-10) {
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
        F: Fn(AtlasComplex<f64>) -> AtlasComplex<f64>, // L_hat(k) - spatial operator in frequency domain
    {
        let mut amplification_factors = Vec::with_capacity(wave_numbers.len());
        let mut max_amplification = zero::<T>();
        let mut critical_wave_number = zero::<T>();

        for &k in wave_numbers {
            let k_complex = AtlasComplex::new(0.0, to_f64(k));

            // Compute spatial operator in frequency domain
            let l_hat = spatial_operator(k_complex);

            // Amplification factor for forward Euler: g = 1 + dt * L_hat(k)
            let dt_f64 = to_f64(dt);
            let g = AtlasComplex::new(1.0, 0.0) + AtlasComplex::new(dt_f64, 0.0) * l_hat;

            let amplification = from_f64(g.norm());
            amplification_factors.push(amplification);

            if amplification > max_amplification {
                max_amplification = amplification;
                critical_wave_number = k;
            }
        }

        let is_stable = max_amplification <= from_f64(1.0001); // Allow small numerical errors

        Ok(VonNeumannAnalysis {
            wave_numbers: wave_numbers.to_vec(),
            amplification_factors,
            max_amplification,
            critical_wave_number,
            is_stable,
            stability_margin: from_f64::<T>(1.0) - max_amplification,
        })
    }

    /// Perform von Neumann stability analysis using an explicit Runge-Kutta method.
    ///
    /// For a linear semi-discrete system $u_t = \hat L(k) u$, an explicit RK method has
    /// amplification factor $g(k) = R(\Delta t\,\hat L(k))$.
    pub fn von_neumann_analysis_explicit_rk<F>(
        &self,
        a: &Array2<T>,
        b: &Array1<T>,
        c: &Array1<T>,
        spatial_operator: F,
        dt: T,
        wave_numbers: &[T],
    ) -> Result<VonNeumannAnalysis<T>>
    where
        F: Fn(AtlasComplex<f64>) -> AtlasComplex<f64>,
    {
        let dt_f64 = to_f64(dt);

        let mut amplification_factors = Vec::with_capacity(wave_numbers.len());
        let mut max_amplification = zero::<T>();
        let mut critical_wave_number = zero::<T>();

        for &k in wave_numbers {
            let k_complex = AtlasComplex::new(0.0, to_f64(k));
            let l_hat = spatial_operator(k_complex);
            let z = AtlasComplex::new(dt_f64, 0.0) * l_hat;

            let g = self.compute_rk_stability_function(a, b, c, z)?;
            let amplification = from_f64(g.norm());
            amplification_factors.push(amplification);

            if amplification > max_amplification {
                max_amplification = amplification;
                critical_wave_number = k;
            }
        }

        let is_stable = max_amplification <= from_f64(1.0001);

        Ok(VonNeumannAnalysis {
            wave_numbers: wave_numbers.to_vec(),
            amplification_factors,
            max_amplification,
            critical_wave_number,
            is_stable,
            stability_margin: from_f64::<T>(1.0) - max_amplification,
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
        F: Fn(AtlasComplex<f64>) -> AtlasComplex<f64>,
    {
        match scheme {
            NumericalScheme::ForwardEuler => {
                self.von_neumann_analysis(spatial_operator, dt, wave_numbers)
            }
            NumericalScheme::RK3 => {
                // Heun/Kutta 3rd-order as used in validation
                let a = matrix_from_row_slice(
                    3,
                    3,
                    &[
                        zero(),
                        zero(),
                        zero(),
                        from_f64(1.0 / 3.0),
                        zero(),
                        zero(),
                        zero(),
                        from_f64(2.0 / 3.0),
                        zero(),
                    ],
                );
                let b = vector_from_vec(vec![from_f64(0.25), zero(), from_f64(0.75)]);
                let c = vector_from_vec(vec![zero(), from_f64(1.0 / 3.0), from_f64(2.0 / 3.0)]);
                self.von_neumann_analysis_explicit_rk(
                    &a,
                    &b,
                    &c,
                    spatial_operator,
                    dt,
                    wave_numbers,
                )
            }
            NumericalScheme::RK4 => {
                let one_half = from_f64(0.5);
                let a = matrix_from_row_slice(
                    4,
                    4,
                    &[
                        zero(),
                        zero(),
                        zero(),
                        zero(),
                        one_half,
                        zero(),
                        zero(),
                        zero(),
                        zero(),
                        one_half,
                        zero(),
                        zero(),
                        zero(),
                        zero(),
                        one(),
                        zero(),
                    ],
                );
                let b = vector_from_vec(vec![
                    from_f64(1.0 / 6.0),
                    from_f64(1.0 / 3.0),
                    from_f64(1.0 / 3.0),
                    from_f64(1.0 / 6.0),
                ]);
                let c = vector_from_vec(vec![zero(), one_half, one_half, one()]);
                self.von_neumann_analysis_explicit_rk(
                    &a,
                    &b,
                    &c,
                    spatial_operator,
                    dt,
                    wave_numbers,
                )
            }
            _ => Err(Error::InvalidInput(
                "von Neumann analysis is only implemented for explicit schemes".to_string(),
            )),
        }
    }
}
