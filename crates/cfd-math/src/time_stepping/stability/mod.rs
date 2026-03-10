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

use nalgebra::RealField;
use num_traits::ToPrimitive;

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
            max_z: T::from_f64(10.0).unwrap_or_else(num_traits::Zero::zero),
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
        let cfl_number = velocity.abs() * dt / dx;
        let max_cfl = scheme.max_cfl_number();

        let stability = if cfl_number <= max_cfl {
            StabilityStatus::Stable
        } else if cfl_number <= max_cfl * T::from_f64(1.5).unwrap_or_else(num_traits::Zero::zero) {
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
            let ratio = (cfl / max_cfl).to_f64().unwrap();
            recommendations.push(format!(
                "CFL number ({:.3}) exceeds stability limit ({:.3}) by {:.1}x",
                cfl.to_f64().unwrap(),
                max_cfl.to_f64().unwrap(),
                ratio
            ));
            recommendations.push("Reduce time step or increase grid resolution".to_string());
            recommendations.push("Consider implicit methods for stiff problems".to_string());
        } else if cfl > max_cfl * T::from_f64(0.8).unwrap_or_else(num_traits::Zero::zero) {
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

impl<T: RealField + Copy + ToPrimitive> Default for StabilityAnalyzer<T> {
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
    pub fn max_cfl_number<T: RealField + Copy>(&self) -> T {
        match self {
            NumericalScheme::RK3 => T::from_f64(1.7).unwrap_or_else(num_traits::Zero::zero),
            NumericalScheme::RK4 => T::from_f64(2.8).unwrap_or_else(num_traits::Zero::zero),
            NumericalScheme::ForwardEuler
            | NumericalScheme::LaxWendroff
            | NumericalScheme::Upwind => T::from_f64(1.0).unwrap_or_else(num_traits::Zero::zero),
            NumericalScheme::CentralDifference => {
                T::from_f64(0.0).unwrap_or_else(num_traits::Zero::zero)
            } // Unstable
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
        let a = nalgebra::DMatrix::from_row_slice(
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
        let a1 = nalgebra::DMatrix::from_row_slice(1, 1, &[0.0]);
        let b1 = DVector::from_vec(vec![1.0]);
        let c1 = DVector::from_vec(vec![0.0]);
        let lim1 = analyzer
            .compute_rk_absolute_stability_limit(&a1, &b1, &c1)
            .unwrap();
        assert_relative_eq!(lim1, 2.0, epsilon = 1e-6);

        // RK3 (Kutta/Heun 3rd order as used in validation)
        let a3 = nalgebra::DMatrix::from_row_slice(
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
        let a4 = nalgebra::DMatrix::from_row_slice(
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
        use num_complex::Complex as NumComplex;
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
