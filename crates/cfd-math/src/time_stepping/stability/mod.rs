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

mod analysis;

use eunomia::{FloatElement, NumericElement, RealField};
use leto::{Array1, Array2};

fn from_f64<T: FloatElement>(value: f64) -> T {
    <T as FloatElement>::from_f64(value)
}

fn to_f64<T: NumericElement>(value: T) -> f64 {
    <T as NumericElement>::to_f64(value)
}

fn abs<T: NumericElement>(value: T) -> T {
    <T as NumericElement>::abs(value)
}

fn zero<T: NumericElement>() -> T {
    <T as NumericElement>::ZERO
}

fn one<T: NumericElement>() -> T {
    <T as NumericElement>::ONE
}

fn vector_from_vec<T>(values: Vec<T>) -> Array1<T> {
    Array1::from_shape_vec([values.len()], values)
        .expect("invariant: vector length matches Leto rank-1 shape")
}

fn matrix_from_row_slice<T: Copy>(rows: usize, cols: usize, values: &[T]) -> Array2<T> {
    Array2::from_shape_vec([rows, cols], values.to_vec())
        .expect("invariant: row-major slice length matches Leto rank-2 shape")
}

/// Stability analysis for time-stepping schemes
#[derive(Debug, Clone)]
pub struct StabilityAnalyzer<T: RealField + Copy> {
    /// Number of points for stability region boundary computation
    resolution: usize,
    /// Maximum |z| value for stability region analysis
    max_z: T,
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FloatElement> StabilityAnalyzer<T> {
    /// Create new stability analyzer with default parameters
    pub fn new() -> Self {
        Self {
            resolution: 200,
            max_z: from_f64(10.0),
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
        let cfl_number = abs(velocity) * dt / dx;
        let max_cfl = scheme.max_cfl_number();

        let stability = if cfl_number <= max_cfl {
            StabilityStatus::Stable
        } else if cfl_number <= max_cfl * from_f64(1.5) {
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

    /// Generate CFL condition recommendations
    fn generate_cfl_recommendations(&self, cfl: T, max_cfl: T) -> Vec<String> {
        let mut recommendations = Vec::new();

        if cfl > max_cfl {
            let ratio = to_f64(cfl / max_cfl);
            recommendations.push(format!(
                "CFL number ({:.3}) exceeds stability limit ({:.3}) by {:.1}x",
                to_f64(cfl),
                to_f64(max_cfl),
                ratio
            ));
            recommendations.push("Reduce time step or increase grid resolution".to_string());
            recommendations.push("Consider implicit methods for stiff problems".to_string());
        } else if cfl > max_cfl * from_f64(0.8) {
            recommendations.push(format!(
                "CFL number ({:.3}) approaching stability limit ({:.3})",
                to_f64(cfl),
                to_f64(max_cfl)
            ));
            recommendations.push("Monitor solution for oscillations".to_string());
        } else {
            recommendations.push(format!(
                "CFL number ({:.3}) well within stability limit ({:.3})",
                to_f64(cfl),
                to_f64(max_cfl)
            ));
        }

        recommendations
    }
}

impl<T: RealField + Copy + FloatElement> Default for StabilityAnalyzer<T> {
    fn default() -> Self {
        Self::new()
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
    pub fn max_cfl_number<T: RealField + Copy + FloatElement>(&self) -> T {
        match self {
            NumericalScheme::RK3 => from_f64(1.7),
            NumericalScheme::RK4 => from_f64(2.8),
            NumericalScheme::ForwardEuler
            | NumericalScheme::LaxWendroff
            | NumericalScheme::Upwind => from_f64(1.0),
            NumericalScheme::CentralDifference => from_f64(0.0), // Unstable
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

#[cfg(test)]
mod tests {
    use super::*;

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

        use eunomia::assert_relative_eq;
        assert_relative_eq!(result.cfl_number, 3.0, epsilon = 1e-10);
        assert!(matches!(result.stability, StabilityStatus::Unstable));
    }

    #[test]
    fn test_rk4_butcher_tableau() {
        let analyzer = StabilityAnalyzer::<f64>::new();

        // Classic RK4 Butcher tableau
        let a = matrix_from_row_slice(
            4,
            4,
            &[
                0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            ],
        );

        let b = vector_from_vec(vec![1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]);
        let c = vector_from_vec(vec![0.0, 0.5, 0.5, 1.0]);

        let result = analyzer.compute_rk_stability_region(&a, &b, &c);
        assert!(result.is_ok());

        let region = result.unwrap();
        assert_eq!(region.method_info.order, 4);
        assert_eq!(region.method_info.stages, 4);
        assert!(!region.boundary.is_empty());
    }

    #[test]
    fn test_absolute_stability_limits_rk1_rk3_rk4() {
        use eunomia::assert_relative_eq;

        let analyzer = StabilityAnalyzer::<f64>::with_params(128, 10.0);

        // RK1 (Forward Euler)
        let a1 = matrix_from_row_slice(1, 1, &[0.0]);
        let b1 = vector_from_vec(vec![1.0]);
        let c1 = vector_from_vec(vec![0.0]);
        let lim1 = analyzer
            .compute_rk_absolute_stability_limit(&a1, &b1, &c1)
            .unwrap();
        assert_relative_eq!(lim1, 2.0, epsilon = 1e-6);

        // RK3 (Kutta/Heun 3rd order as used in validation)
        let a3 = matrix_from_row_slice(
            3,
            3,
            &[0.0, 0.0, 0.0, 1.0 / 3.0, 0.0, 0.0, 0.0, 2.0 / 3.0, 0.0],
        );
        let b3 = vector_from_vec(vec![0.25, 0.0, 0.75]);
        let c3 = vector_from_vec(vec![0.0, 1.0 / 3.0, 2.0 / 3.0]);
        let lim3 = analyzer
            .compute_rk_absolute_stability_limit(&a3, &b3, &c3)
            .unwrap();
        assert!(lim3 > 2.45 && lim3 < 2.58);

        // Classic RK4
        let a4 = matrix_from_row_slice(
            4,
            4,
            &[
                0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 0.5, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0,
            ],
        );
        let b4 = vector_from_vec(vec![1.0 / 6.0, 1.0 / 3.0, 1.0 / 3.0, 1.0 / 6.0]);
        let c4 = vector_from_vec(vec![0.0, 0.5, 0.5, 1.0]);
        let lim4 = analyzer
            .compute_rk_absolute_stability_limit(&a4, &b4, &c4)
            .unwrap();
        assert!(lim4 > 2.7 && lim4 < 2.9);
    }

    #[test]
    fn test_von_neumann_rk4_constant_negative_operator() {
        // If L_hat = -1, z = -dt, and RK4's amplification should be close to exp(-dt)
        // for modest dt. We only assert a loose bound to avoid overfitting.
        use eunomia::Complex as AtlasComplex;
        let analyzer = StabilityAnalyzer::<f64>::with_params(64, 10.0);
        let dt = 1.0;
        let wave_numbers = vec![1.0, 2.0, 3.0];

        let spatial_operator = |_k: AtlasComplex<f64>| AtlasComplex::new(-1.0, 0.0);

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
