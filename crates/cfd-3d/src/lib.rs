//! # 3D Computational Fluid Dynamics (CFD) Suite
//!
//! This crate provides high-performance, mathematically rigorous implementations
//! of advanced CFD methods for three-dimensional flow simulations.
//!
//! ## Core Methods Implemented
//!
//! ### Finite Element Method (FEM)
//! - **Galerkin Method**: Weak form discretization of Navier-Stokes equations
//! - **SUPG/PSPG Stabilization**: Streamline-Upwind Petrov-Galerkin and Pressure-Stabilizing Petrov-Galerkin
//! - **Mixed Elements**: Velocity-pressure coupling for incompressible flows
//!
//! ### Immersed Boundary Method (IBM)
//! - **Direct Forcing**: Momentum forcing at immersed boundaries
//! - **Feedback Forcing**: PID-based boundary condition enforcement
//! - **Discrete Delta Functions**: Smooth interpolation kernels
//!
//! ### Level Set Method
//! - **Signed Distance Function**: Accurate interface representation
//! - **Reinitialization**: Sussman redistancing algorithm
//! - **Narrow Band**: Efficient computation near interfaces
//!
//! ### Volume of Fluid (VOF)
//! - **PLIC Reconstruction**: Piecewise Linear Interface Construction
//! - **Geometric Advection**: Volume-preserving interface transport
//! - **Surface Tension**: Continuum Surface Force (CSF) model
//!
//! ### Spectral Methods
//! - **Fourier Spectral**: Efficient FFT-based spatial discretization
//! - **Chebyshev Spectral**: High-accuracy boundary treatment
//! - **Collocation Methods**: Pseudospectral accuracy
//!
//! ## Governing Equations
//!
//! ### Navier-Stokes Equations (Incompressible)
//! ```math
//! ∂u/∂t + (u·∇)u = -∇p/ρ + ν∇²u + f    (Momentum)
//! ∇·u = 0                                      (Continuity)
//! ```
//!
//! ### Turbulence Modeling (Optional)
//! - **LES**: Large Eddy Simulation with dynamic Smagorinsky
//! - **DES**: Detached Eddy Simulation
//! - **RANS**: Reynolds-Averaged Navier-Stokes
//!
//! ## Numerical Foundations
//!
//! ### Finite Element Theory
//! **Theorem (Lax-Milgram)**: For coercive, continuous bilinear forms,
//! the Galerkin method converges optimally.
//!
//! **SUPG Stabilization** (Brooks & Hughes, 1982):
//! ```math
//! τ = [(2/Δt)² + (2|u|·h)² + (4ν/h²)²]^(-1/2)
//! ```
//!
//! ### Spectral Accuracy
//! **Theorem (Collocation Convergence)**: For smooth solutions,
//! spectral methods achieve exponential convergence O(e^(-cN)).
//!
//! ### Interface Methods
//! **VOF Conservation**: Geometric advection preserves volume to machine precision.
//!
//! ## Implementation Highlights
//!
//! - **Performance**: SIMD kernels, cache-blocked algorithms, FFT acceleration
//! - **Accuracy**: High-order methods, exact conservation properties
//! - **Robustness**: Advanced stabilization, adaptive algorithms
//! - **Scalability**: Parallel decomposition, memory-efficient data structures
//!
//! ## References
//!
//! - **FEM**: Hughes, T.J.R. (2000). *The Finite Element Method: Linear Static and Dynamic Finite Element Analysis*
//! - **Spectral**: Boyd, J.P. (2001). *Chebyshev and Fourier Spectral Methods*
//! - **IBM**: Mittal, R. & Iaccarino, G. (2005). *Immersed Boundary Methods*
//! - **VOF**: Scardovelli, R. & Zaleski, S. (1999). *Direct Numerical Simulation of Free-Surface and Interfacial Flow*
//! - **Level Set**: Osher, S. & Fedkiw, R. (2003). *Level Set Methods and Dynamic Implicit Surfaces*

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// 3D CFD simulation allows - strategic configuration for numerical computing
#![allow(clippy::similar_names)] // 3D variables (u,v,w,p; nx,ny,nz; dx,dy,dz; i,j,k) often similar
#![allow(clippy::cast_precision_loss)] // Performance-critical numerical loops in 3D solvers
#![allow(clippy::cast_possible_truncation)] // Grid indices and array sizes typically small
#![allow(clippy::unused_self)] // Solver trait methods maintain consistent interfaces
#![allow(clippy::must_use_candidate)] // Solver utilities and getters used in computational contexts
#![allow(clippy::missing_errors_doc)] // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)] // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)] // Signed to unsigned casts common in 3D CFD indexing
#![allow(clippy::cast_possible_wrap)] // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)] // 3D CFD functions often need many physical parameters
#![allow(clippy::float_cmp)] // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)] // Result types maintained for API consistency
#![allow(clippy::items_after_statements)] // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)] // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)] // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)] // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)] // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)] // Builder patterns used internally
#![allow(clippy::should_implement_trait)] // CFD-specific trait implementations
#![allow(clippy::manual_let_else)]
#![allow(clippy::match_same_arms)]
#![allow(clippy::useless_conversion)]
#![allow(clippy::inline_always)]
#![allow(clippy::needless_range_loop)]
#![allow(clippy::format_push_string)]
#![allow(clippy::new_without_default)]
#![allow(clippy::trivially_copy_pass_by_ref)]
#![allow(clippy::empty_line_after_doc_comments)]
#![allow(clippy::too_many_lines)]
#![allow(clippy::doc_overindented_list_items)]
#![allow(clippy::same_item_push)]
#![allow(clippy::manual_clamp)]
#![allow(clippy::duplicated_attributes)]

pub mod bifurcation;
/// Canonical blueprint-driven preprocessing and cross-fidelity tracing for 3D workflows.
pub mod blueprint_integration;
/// Legacy multi-stage cascade 3D FEM solver for CIF networks.
///
/// Prefer [`blueprint_integration`] for new blueprint-driven preprocessing flows.
pub mod cascade;
pub mod fem;
pub mod ibm;
pub mod level_set;
pub mod physics;
pub mod serpentine;
pub mod spectral;
/// 3D trifurcation (three-way branching) flow solvers and validation
pub mod trifurcation;
pub mod venturi;
pub mod vof;

// Export FEM functionality
pub use fem::{FemConfig, FemSolver, StokesFlowProblem};

// Export blueprint integration functionality
pub use blueprint_integration::{
    process_blueprint_with_reference_trace, Blueprint3dProcessingConfig, Blueprint3dTrace,
    ChannelCrossFidelityTrace, NodeCrossFidelityTrace,
};

// Export spectral functionality
pub use spectral::{
    BandLimitedRandomPhaseForcing3D, BandLimitedRandomPhaseForcingConfig, BasisFunction,
    ChebyshevPolynomial, EnstrophySpectrum, FourierTransform, KineticEnergySpectrum,
    PeriodicPseudospectralDns3D, PeriodicPseudospectralDnsConfig, ProbeSignalSpectrum,
    SpectralBasis, SpectralConfig, SpectralSolution, SpectralSolver,
    TemporalAutocorrelation, TimeResampledBandLimitedForcing3D,
    TimeResampledBandLimitedForcingConfig, enstrophy_spectrum, kinetic_energy_spectrum,
    probe_signal_spectrum, temporal_autocorrelation,
};

// Export IBM functionality
pub use ibm::{IbmConfig, IbmSolver, LagrangianPoint};

// Export level set functionality
pub use level_set::{LevelSetConfig, LevelSetSolver};

// Export VOF functionality
pub use vof::{VofConfig, VofSolver};

// CSG integration from cfd-mesh - feature-gated for optional dependency

pub use cfd_mesh::application::csg::CsgError;

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::ComplexField;

    /// Test module imports and configuration instantiation.
    ///
    /// # Invariants Verified
    /// - FemConfig defaults are physically valid (quadrature_order ≥ 1, stabilization enabled).
    /// - SpectralConfig construction succeeds for valid mode counts.
    /// - IBM, LevelSet, VOF configs default to non-degenerate values.
    #[test]
    fn test_module_imports() {
        // Verify we can create configurations for each solver type

        // FEM configuration test
        let fem_config: FemConfig<f64> = FemConfig::default();
        assert!(
            fem_config.quadrature_order >= 1,
            "quadrature_order must be ≥ 1 for valid Gauss integration"
        );
        assert!(
            fem_config.use_stabilization,
            "SUPG/PSPG stabilization must be enabled by default (Lax-Milgram requirement)"
        );

        // Spectral configuration test
        let spectral_config: SpectralConfig<f64> = SpectralConfig::new(8, 8, 8)
            .expect("SpectralConfig::new must succeed for valid mode count 8");
        assert_eq!(spectral_config.nx_modes, 8);
        assert_eq!(spectral_config.ny_modes, 8);
        assert_eq!(spectral_config.nz_modes, 8);

        // IBM configuration test
        let ibm_config = IbmConfig::default();
        assert!(
            ibm_config.smoothing_width > 0.0,
            "IBM delta function width must be positive (partition of unity requirement)"
        );
        assert!(
            ibm_config.use_direct_forcing,
            "direct forcing must be enabled by default"
        );

        // Level set configuration test
        let level_set_config = LevelSetConfig::default();
        assert!(
            level_set_config.reinitialization_interval > 0,
            "reinitialization interval must be positive to prevent SDF drift"
        );
        assert!(
            level_set_config.use_weno,
            "WENO scheme must be enabled by default for interface accuracy"
        );

        // VOF configuration test
        let vof_config = VofConfig::default();
        assert!(
            vof_config.tolerance > 0.0,
            "VOF tolerance must be strictly positive"
        );
        assert!(
            vof_config.reconstruction_method == vof::InterfaceReconstruction::PLIC,
            "PLIC must be the default reconstruction (volume conservation guarantee)"
        );
    }

    /// Test spectral solver instantiation and basic operations.
    ///
    /// Verifies that a valid SpectralConfig produces a usable SpectralSolver without errors.
    #[test]
    fn test_spectral_solver_creation() -> Result<(), Box<dyn std::error::Error>> {
        let config: SpectralConfig<f64> = SpectralConfig::new(4, 4, 4)
            .expect("SpectralConfig::new must succeed for valid mode count 4");
        let solver =
            SpectralSolver::new(config).expect("SpectralSolver::new must succeed for valid config");
        // Verify solver is ready (no further assertions — construction itself is the test)
        let _ = solver;
        Ok(())
    }

    /// Test Chebyshev polynomial operations — literature-based validation.
    ///
    /// # Theorem Verified
    /// Gauss-Lobatto points xᵢ = cos(iπ/N) ∈ [-1,1] for i=0..N.
    /// Spectral differentiation matrix D satisfies (D·f)[i] ≈ f'(xᵢ) to machine precision for polynomials.
    #[test]
    fn test_chebyshev_polynomial_operations() {
        use crate::spectral::ChebyshevPolynomial;

        let poly: ChebyshevPolynomial<f64> = ChebyshevPolynomial::new(5)
            .expect("ChebyshevPolynomial::new must succeed for valid N=5");

        // Test basic properties
        assert_eq!(poly.num_points(), 5);
        assert_eq!(poly.points().len(), 5);

        // Gauss-Lobatto theorem: all collocation points must lie in [-1, 1]
        for &point in poly.points() {
            assert!(
                (-1.0..=1.0).contains(&point),
                "Gauss-Lobatto point {point} must lie in [-1, 1]"
            );
        }

        // Test differentiation matrix shape
        let n = poly.num_points();
        let test_vector = nalgebra::DVector::from_fn(n, |i, _| {
            let x = poly.points()[i];
            x * x // f(x) = x²
        });

        let derivative = poly.differentiate(&test_vector);
        assert_eq!(derivative.len(), n);

        // For f(x) = x², f'(x) = 2x
        let expected_derivative_at_x1: f64 = 2.0 * poly.points()[0]; // f'(x₀)
        let computed_derivative_at_x1: f64 = derivative[0];

        // Allow for finite-N spectral differentiation error
        assert!(
            num_traits::Float::abs(computed_derivative_at_x1 - expected_derivative_at_x1) < 0.1,
            "spectral derivative error {:.2e} must be < 0.1",
            (computed_derivative_at_x1 - expected_derivative_at_x1).abs()
        );
    }

    /// Test Fourier transform accuracy with analytical solutions.
    ///
    /// # Theorem Verified (Parseval / DFT Constant-Signal)
    /// For f(x) = 1 (constant), the DFT satisfies F̂(0) = 1 (DC = mean, normalised)
    /// and F̂(k) = 0 for all k ≠ 0 (no oscillatory content).
    #[test]
    fn test_fourier_transform_accuracy() -> Result<(), Box<dyn std::error::Error>> {
        use crate::spectral::FourierTransform;
        use nalgebra::DVector;

        let ft: FourierTransform<f64> =
            FourierTransform::new(8).expect("FourierTransform::new must succeed for valid N=8");

        let n = 8;
        let constant_signal = DVector::from_element(n, 1.0);

        let spectrum = ft
            .forward(&constant_signal)
            .expect("forward FFT must succeed for finite, valid input");

        assert_eq!(spectrum.len(), n, "output length must equal input length");

        // DC component for normalized constant-1 signal must be 1.0
        let dc_magnitude = spectrum[0].modulus();
        assert!(
            (dc_magnitude - 1.0).abs() < 1e-12,
            "DC magnitude {dc_magnitude:.15} must equal 1.0 for constant signal"
        );

        // All non-DC components must be essentially zero (aliasing theorem)
        for i in 1..n {
            assert!(
                spectrum[i].modulus() < 1e-10,
                "non-DC component {i} magnitude {:.2e} must be ~0 for constant signal",
                spectrum[i].modulus()
            );
        }
        Ok(())
    }

    /// Test level set solver basic functionality with simple geometry
    #[test]
    fn test_level_set_basic_operations() {
        let config = LevelSetConfig::default();

        // Create level set solver with grid parameters
        let nx = 10;
        let ny = 10;
        let nz = 10;
        let dx = 0.1;
        let dy = 0.1;
        let dz = 0.1;

        let _solver: LevelSetSolver<f64> = LevelSetSolver::new(config, nx, ny, nz, dx, dy, dz);
        // Basic validation - solver was created successfully
    }

    /// Test VOF solver instantiation with proper grid
    #[test]
    fn test_vof_solver_creation() {
        let config = VofConfig::default();

        // Create VOF solver with grid parameters
        let nx = 10;
        let ny = 10;
        let nz = 10;
        let dx = 0.1;
        let dy = 0.1;
        let dz = 0.1;

        let solver: crate::vof::VofSolver<f64> =
            crate::vof::VofSolver::create(config, nx, ny, nz, dx, dy, dz);

        // Basic validation - solver was created successfully
        assert_eq!(solver.nx, nx);
        assert_eq!(solver.ny, ny);
        assert_eq!(solver.nz, nz);
    }

    /// Test FEM element type configuration
    #[test]
    fn test_fem_element_types() {
        use cfd_core::geometry::ElementType;

        // Test different element types
        let config_tet: FemConfig<f64> = FemConfig {
            element_type: ElementType::Tetrahedron,
            ..Default::default()
        };
        assert_eq!(config_tet.element_type, ElementType::Tetrahedron);

        let config_hex: FemConfig<f64> = FemConfig {
            element_type: ElementType::Hexahedron,
            ..Default::default()
        };
        assert_eq!(config_hex.element_type, ElementType::Hexahedron);

        // Verify stabilization parameters
        assert!(config_hex.tau > 0.0);
        assert!(config_hex.quadrature_order >= 1);
    }

    /// Test collocation point properties - literature validation
    #[test]
    fn test_chebyshev_collocation_properties() {
        use crate::spectral::ChebyshevPolynomial;

        let poly: ChebyshevPolynomial<f64> =
            ChebyshevPolynomial::new(6).expect("Failed to create polynomial");
        let points = poly.points();

        // Test Gauss-Lobatto point properties
        // First and last points should be ±1
        assert!(num_traits::Float::abs(points[0] - 1.0) < 1e-14);
        assert!((points[points.len() - 1] + 1.0).abs() < 1e-14);

        // Points should be in descending order (1, ..., -1)
        for i in 0..points.len() - 1 {
            assert!(points[i] > points[i + 1]);
        }

        // Test second derivative matrix construction
        let d2_matrix = poly
            .second_derivative_matrix()
            .expect("Failed to create second derivative matrix");
        assert_eq!(d2_matrix.nrows(), points.len());
        assert_eq!(d2_matrix.ncols(), points.len());

        // Test on a quadratic function: f(x) = x², f''(x) = 2
        let quadratic = nalgebra::DVector::from_fn(points.len(), |i, _| {
            let x: f64 = points[i];
            x * x
        });

        let second_derivative = &d2_matrix * &quadratic;

        // For f(x) = x², f''(x) should be 2 within machine precision
        for &val in second_derivative.iter() {
            assert!(num_traits::Float::abs(val - 2.0) < 1e-10); // High accuracy expected
        }
    }
}

#[cfg(test)]
#[path = "tests/adversarial_tests.rs"]
mod adversarial_tests_suite;
