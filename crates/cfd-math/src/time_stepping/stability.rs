//! Stability analysis for time-stepping schemes in CFD
//!
//! This module provides comprehensive stability analysis including:
//! - Stability region computation for Runge-Kutta methods
//! - CFL condition verification for advection-dominated flows
//! - Von Neumann stability analysis for linear PDEs
//! - Stability region plotting and visualization
//!
//! References:
//! - Hairer & Nørsett (1993): Solving Ordinary Differential Equations I
//! - Trefthen (1996): Finite Difference and Spectral Methods for Ordinary and Partial Differential Equations
//! - LeVeque (2002): Finite Volume Methods for Hyperbolic Problems

use cfd_core::error::{Error, Result};
use nalgebra::{DMatrix, DVector, RealField};
use num_complex::Complex as NumComplex;
use num_traits::ToPrimitive;
use std::f64::consts::PI;

/// Stability analysis for time-stepping schemes
#[derive(Debug, Clone)]
pub struct StabilityAnalyzer<T: RealField + Copy> {
    /// Number of points for stability region boundary computation
    resolution: usize,
    /// Maximum |z| value for stability region analysis
    max_z: T,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + ToPrimitive> StabilityAnalyzer<T> {
    /// Create new stability analyzer with default parameters
    pub fn new() -> Self {
        Self {
            resolution: 200,
            max_z: T::from_f64(10.0).unwrap(),
            _phantom: std::marker::PhantomData,
        }
    }

    /// Create with custom parameters
    pub fn with_params(resolution: usize, max_z: T) -> Self {
        Self {
            resolution,
            max_z,
            _phantom: std::marker::PhantomData,
        }
    }

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
                Error::InvalidInput(
                    "stability boundary real part not representable in T".to_string(),
                )
            })?;
            let z_imag = T::from_f64(z.im).ok_or_else(|| {
                Error::InvalidInput(
                    "stability boundary imag part not representable in T".to_string(),
                )
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

        // Check if method is explicit (strictly lower triangular A)
        let mut is_explicit = true;
        for i in 0..s {
            for j in i..s {
                if a[(i, j)] != T::zero() {
                    is_explicit = false;
                    break;
                }
            }
            if !is_explicit {
                break;
            }
        }

        if is_explicit {
            // For explicit methods, R(z) = 1 + z*b^T * (I - z*A)^(-1) * 1
            // We can compute w = (I - z*A)^(-1) * 1 using forward substitution
            // w = 1 + z*A*w
            // w_i = 1 + z * sum_{j=0}^{i-1} a_{ij} * w_j

            let mut w = vec![NumComplex::new(0.0, 0.0); s];
            let one = NumComplex::new(1.0, 0.0);

            for i in 0..s {
                let mut sum_aw = NumComplex::new(0.0, 0.0);
                for j in 0..i {
                    // a[(i, j)] is T, convert to f64 then complex
                    let a_val = a[(i, j)].to_f64().unwrap();
                    sum_aw += NumComplex::new(a_val, 0.0) * w[j];
                }
                w[i] = one + z * sum_aw;
            }

            // R(z) = 1 + z * sum_{i=0}^{s-1} b_i * w_i
            let mut sum_bw = NumComplex::new(0.0, 0.0);
            for i in 0..s {
                let b_val = b[i].to_f64().unwrap();
                sum_bw += NumComplex::new(b_val, 0.0) * w[i];
            }

            let r_z = one + z * sum_bw;
            Ok(r_z)
        } else {
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
                    let ones =
                        DVector::<NumComplex<f64>>::from_element(s, NumComplex::new(1.0, 0.0));
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

    /// Analyze CFL condition for advection equation: ∂u/∂t + a ∂u/∂x = 0
    ///
    /// CFL condition: |a| * Δt / Δx <= C (where C depends on scheme)
    ///
    /// # Arguments
    /// * `velocity` - Maximum advection velocity |a|
    /// * `dt` - Time step size
    /// * `dx` - Spatial grid spacing
    /// * `scheme` - Numerical scheme type
    ///
    /// # Returns
    /// CFL analysis result
    pub fn analyze_cfl_condition(
        &self,
        velocity: T,
        dt: T,
        dx: T,
        scheme: NumericalScheme,
    ) -> CFLAnalysis<T> {
        let cfl_number = velocity.abs() * dt / dx;
        let max_cfl = scheme.max_cfl_number();

        let stability = if cfl_number <= max_cfl {
            StabilityStatus::Stable
        } else if cfl_number <= max_cfl * T::from_f64(1.5).unwrap() {
            StabilityStatus::MarginallyStable
        } else {
            StabilityStatus::Unstable
        };

        CFLAnalysis {
            cfl_number,
            max_cfl,
            stability,
            scheme,
            recommendations: self.generate_cfl_recommendations(cfl_number, max_cfl),
        }
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
            NumericalScheme::ForwardEuler => {
                self.von_neumann_analysis(spatial_operator, dt, wave_numbers)
            }
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

    /// Generate CFL condition recommendations
    fn generate_cfl_recommendations(&self, cfl: T, max_cfl: T) -> Vec<String> {
        let mut recommendations = Vec::new();

        if cfl > max_cfl {
            let ratio = (cfl / max_cfl).to_f64().unwrap();
            recommendations.push(format!(
                "CFL number ({:.3}) exceeds stability limit ({:.3}) by {:.1}x",
                cfl.to_f64().unwrap(),
                max_cfl.to_f64().unwrap(),
                ratio
            ));
            recommendations.push("Reduce time step or increase grid resolution".to_string());
            recommendations.push("Consider implicit methods for stiff problems".to_string());
        } else if cfl > max_cfl * T::from_f64(0.8).unwrap() {
            recommendations.push(format!(
                "CFL number ({:.3}) approaching stability limit ({:.3})",
                cfl.to_f64().unwrap(),
                max_cfl.to_f64().unwrap()
            ));
            recommendations.push("Monitor solution for oscillations".to_string());
        } else {
            recommendations.push(format!(
                "CFL number ({:.3}) well within stability limit ({:.3})",
                cfl.to_f64().unwrap(),
                max_cfl.to_f64().unwrap()
            ));
        }

        recommendations
    }
}

/// Stability region representation
#[derive(Debug, Clone)]
pub struct StabilityRegion<T: RealField + Copy> {
    /// Boundary points of stability region
    pub boundary: Vec<ComplexPoint<T>>,
    /// Interior test points
    pub interior: Vec<ComplexPoint<T>>,
    /// Method information
    pub method_info: MethodInfo,
}

/// Complex plane point with stability information
#[derive(Debug, Clone)]
pub struct ComplexPoint<T: RealField + Copy> {
    /// Real part
    pub real: T,
    /// Imaginary part
    pub imag: T,
    /// Whether this point is in the stability region
    pub stability: bool,
}

/// Method information
#[derive(Debug, Clone)]
pub struct MethodInfo {
    /// Method name
    pub name: String,
    /// Accuracy order
    pub order: usize,
    /// Number of stages
    pub stages: usize,
    /// Stability type
    pub stability_type: StabilityType,
}

/// Stability type classification
#[derive(Debug, Clone, PartialEq)]
pub enum StabilityType {
    /// Explicit methods (conditionally stable)
    Explicit,
    /// Implicit methods (unconditionally stable)
    Implicit,
    /// A-stable methods (stable for all |z|)
    AStable,
}

/// CFL condition analysis result
#[derive(Debug, Clone)]
pub struct CFLAnalysis<T: RealField + Copy> {
    /// Computed CFL number
    pub cfl_number: T,
    /// Maximum stable CFL number for this scheme
    pub max_cfl: T,
    /// Stability status
    pub stability: StabilityStatus,
    /// Numerical scheme used
    pub scheme: NumericalScheme,
    /// Recommendations for stability
    pub recommendations: Vec<String>,
}

/// Stability status
#[derive(Debug, Clone, PartialEq)]
pub enum StabilityStatus {
    /// Within stability region
    Stable,
    /// Near stability boundary
    MarginallyStable,
    /// Outside stability region
    Unstable,
}

/// Numerical scheme types for CFL analysis
#[derive(Debug, Clone, PartialEq)]
pub enum NumericalScheme {
    /// Forward Euler (CFL ≤ 1 for advection)
    ForwardEuler,
    /// Explicit Runge-Kutta 3 (CFL ≈ 1.7)
    RK3,
    /// Classic Runge-Kutta 4 (CFL ≈ 2.8)
    RK4,
    /// Lax-Wendroff (CFL ≤ 1)
    LaxWendroff,
    /// Upwind (CFL ≤ 1)
    Upwind,
    /// Central difference (unconditionally unstable)
    CentralDifference,
}

impl NumericalScheme {
    /// Maximum CFL number for stability
    pub fn max_cfl_number<T: RealField + Copy>(&self) -> T {
        match self {
            NumericalScheme::ForwardEuler
            | NumericalScheme::LaxWendroff
            | NumericalScheme::Upwind => T::from_f64(1.0).unwrap(),
            NumericalScheme::RK3 => T::from_f64(1.7).unwrap(),
            NumericalScheme::RK4 => T::from_f64(2.8).unwrap(),
            NumericalScheme::CentralDifference => T::from_f64(0.0).unwrap(), // Unstable
        }
    }
}

/// Von Neumann stability analysis result
#[derive(Debug, Clone)]
pub struct VonNeumannAnalysis<T: RealField + Copy> {
    /// Wave numbers tested
    pub wave_numbers: Vec<T>,
    /// Amplification factors for each wave number
    pub amplification_factors: Vec<T>,
    /// Maximum amplification factor
    pub max_amplification: T,
    /// Wave number with maximum amplification
    pub critical_wave_number: T,
    /// Whether the scheme is stable
    pub is_stable: bool,
    /// Stability margin (1.0 - max_amplification)
    pub stability_margin: T,
}

impl<T: RealField + Copy + ToPrimitive> Default for StabilityAnalyzer<T> {
    fn default() -> Self {
        Self::new()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::DVector;

    #[test]
    fn test_cfl_analysis() {
        let analyzer = StabilityAnalyzer::<f64>::new();

        // Test stable case
        let result = analyzer.analyze_cfl_condition(
            1.0, // velocity
            0.1, // dt
            0.2, // dx
            NumericalScheme::ForwardEuler,
        );

        assert_eq!(result.cfl_number, 0.5);
        assert_eq!(result.max_cfl, 1.0);
        assert!(matches!(result.stability, StabilityStatus::Stable));
    }

    #[test]
    fn test_unstable_cfl() {
        let analyzer = StabilityAnalyzer::<f64>::new();

        // Test unstable case
        let result = analyzer.analyze_cfl_condition(
            1.0, // velocity
            0.3, // dt
            0.1, // dx
            NumericalScheme::ForwardEuler,
        );

        use approx::assert_relative_eq;
        assert_relative_eq!(result.cfl_number, 3.0, epsilon = 1e-10);
        assert!(matches!(result.stability, StabilityStatus::Unstable));
    }

    #[test]
    fn test_rk4_butcher_tableau() {
        let analyzer = StabilityAnalyzer::<f64>::new();

        // Classic RK4 Butcher tableau
        let a = DMatrix::from_row_slice(
            4,
            4,
            &[
                0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            ],
        );

        let b = DVector::from_vec(vec![1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]);
        let c = DVector::from_vec(vec![0.0, 0.5, 0.5, 1.0]);

        let result = analyzer.compute_rk_stability_region(&a, &b, &c);
        assert!(result.is_ok());

        let region = result.unwrap();
        assert_eq!(region.method_info.order, 4);
        assert_eq!(region.method_info.stages, 4);
        assert!(!region.boundary.is_empty());
    }

    #[test]
    fn test_absolute_stability_limits_rk1_rk3_rk4() {
        use approx::assert_relative_eq;

        let analyzer = StabilityAnalyzer::<f64>::with_params(128, 10.0);

        // RK1 (Forward Euler)
        let a1 = DMatrix::from_row_slice(1, 1, &[0.0]);
        let b1 = DVector::from_vec(vec![1.0]);
        let c1 = DVector::from_vec(vec![0.0]);
        let lim1 = analyzer
            .compute_rk_absolute_stability_limit(&a1, &b1, &c1)
            .unwrap();
        assert_relative_eq!(lim1, 2.0, epsilon = 1e-6);

        // RK3 (Kutta/Heun 3rd order as used in validation)
        let a3 = DMatrix::from_row_slice(
            3,
            3,
            &[0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 2.0 / 3.0, 0.0],
        );
        let b3 = DVector::from_vec(vec![0.25, 0.0, 0.75]);
        let c3 = DVector::from_vec(vec![0.0, 1.0 / 3.0, 2.0 / 3.0]);
        let lim3 = analyzer
            .compute_rk_absolute_stability_limit(&a3, &b3, &c3)
            .unwrap();
        assert!(lim3 > 2.45 && lim3 < 2.58);

        // Classic RK4
        let a4 = DMatrix::from_row_slice(
            4,
            4,
            &[
                0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            ],
        );
        let b4 = DVector::from_vec(vec![1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]);
        let c4 = DVector::from_vec(vec![0.0, 0.5, 0.5, 1.0]);
        let lim4 = analyzer
            .compute_rk_absolute_stability_limit(&a4, &b4, &c4)
            .unwrap();
        assert!(lim4 > 2.7 && lim4 < 2.9);
    }

    #[test]
    fn test_von_neumann_rk4_constant_negative_operator() {
        // If L_hat = -1, z = -dt, and RK4's amplification should be close to exp(-dt)
        // for modest dt. We only assert a loose bound to avoid overfitting.
        let analyzer = StabilityAnalyzer::<f64>::with_params(64, 10.0);
        let dt = 1.0;
        let wave_numbers = vec![1.0, 2.0, 3.0];

        let spatial_operator = |_k: NumComplex<f64>| NumComplex::new(-1.0, 0.0);

        let analysis = analyzer
            .von_neumann_analysis_with_scheme(
                NumericalScheme::RK4,
                spatial_operator,
                dt,
                &wave_numbers,
            )
            .unwrap();

        assert!(analysis.is_stable);
        for &amp in &analysis.amplification_factors {
            assert!(amp > 0.30 && amp < 0.50);
        }
    }
}
