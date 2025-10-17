//! 3D CFD simulations.

#![warn(missing_docs)]
#![warn(clippy::all)]
#![warn(clippy::pedantic)]
#![allow(clippy::module_name_repetitions)]
// 3D CFD simulation allows - strategic configuration for numerical computing
#![allow(clippy::similar_names)]           // 3D variables (u,v,w,p; nx,ny,nz; dx,dy,dz; i,j,k) often similar
#![allow(clippy::cast_precision_loss)]     // Performance-critical numerical loops in 3D solvers
#![allow(clippy::cast_possible_truncation)] // Grid indices and array sizes typically small
#![allow(clippy::unused_self)]             // Solver trait methods maintain consistent interfaces
#![allow(clippy::must_use_candidate)]      // Solver utilities and getters used in computational contexts
#![allow(clippy::missing_errors_doc)]      // Error documentation deferred for internal APIs
#![allow(clippy::missing_panics_doc)]      // Panic documentation deferred for internal APIs
#![allow(clippy::cast_sign_loss)]          // Signed to unsigned casts common in 3D CFD indexing
#![allow(clippy::cast_possible_wrap)]      // Wrap-around acceptable for grid indices
#![allow(clippy::too_many_arguments)]      // 3D CFD functions often need many physical parameters
#![allow(clippy::float_cmp)]               // Float comparisons necessary in numerical algorithms
#![allow(clippy::unnecessary_wraps)]        // Result types maintained for API consistency
#![allow(clippy::items_after_statements)]  // Helper functions after statements improve readability
#![allow(clippy::many_single_char_names)]       // Mathematical notation (i,j,k,x,y,z) is standard
#![allow(clippy::unreadable_literal)]      // Long literals used for precise physical constants
#![allow(clippy::redundant_closure_for_method_calls)] // Closures improve readability in numerical pipelines
#![allow(clippy::doc_markdown)]            // Math notation doesn't need backticks
#![allow(clippy::needless_pass_by_value)]  // Pass by value for Copy types is idiomatic
#![allow(clippy::return_self_not_must_use)]  // Builder patterns used internally
#![allow(clippy::ptr_arg)]                 // &Vec used for API compatibility
#![allow(clippy::should_implement_trait)]  // CFD-specific trait implementations

pub mod fem;
pub mod ibm;
pub mod level_set;
pub mod spectral;
pub mod vof;

// Export FEM functionality
pub use fem::{FemConfig, FemSolver, StokesFlowProblem};

// Export spectral functionality
pub use spectral::{
    BasisFunction, ChebyshevPolynomial, FourierTransform, SpectralBasis, SpectralConfig,
    SpectralSolution, SpectralSolver,
};

// Export IBM functionality
pub use ibm::{IbmConfig, IbmSolver, LagrangianPoint};

// Export level set functionality
pub use level_set::{LevelSetConfig, LevelSetSolver};

// Export VOF functionality
pub use vof::{VofConfig, VofSolver};

// CSG integration from cfd-mesh - feature-gated for optional dependency
#[cfg(feature = "csg")]
pub use cfd_mesh::csg::CsgMeshAdapter;

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::ComplexField;

    /// Test module imports and configuration instantiation
    #[test]
    fn test_module_imports() {
        // Verify we can create configurations for each solver type

        // FEM configuration test
        let fem_config: FemConfig<f64> = FemConfig::default();
        assert!(fem_config.quadrature_order >= 1);
        assert!(fem_config.use_stabilization);

        // Spectral configuration test
        let spectral_config: SpectralConfig<f64> =
            SpectralConfig::new(8, 8, 8).expect("Failed to create spectral config");
        assert_eq!(spectral_config.nx_modes, 8);
        assert_eq!(spectral_config.ny_modes, 8);
        assert_eq!(spectral_config.nz_modes, 8);

        // IBM configuration test
        let ibm_config = IbmConfig::default();
        assert!(ibm_config.smoothing_width > 0.0);
        assert!(ibm_config.use_direct_forcing);

        // Level set configuration test
        let level_set_config = LevelSetConfig::default();
        assert!(level_set_config.reinitialization_interval > 0);
        assert!(level_set_config.use_weno);

        // VOF configuration test
        let vof_config = VofConfig::default();
        assert!(vof_config.tolerance > 0.0);
        assert!(vof_config.use_plic);
    }

    /// Test spectral solver instantiation and basic operations
    #[test]
    fn test_spectral_solver_creation() {
        let config: SpectralConfig<f64> =
            SpectralConfig::new(4, 4, 4).expect("Failed to create config");
        let solver = SpectralSolver::new(config);

        assert!(solver.is_ok());
        let _solver = solver.unwrap();
        // Basic validation - solver was created successfully
    }

    /// Test Chebyshev polynomial operations - literature based validation
    #[test]
    fn test_chebyshev_polynomial_operations() {
        use crate::spectral::ChebyshevPolynomial;

        let poly: ChebyshevPolynomial<f64> =
            ChebyshevPolynomial::new(5).expect("Failed to create Chebyshev polynomial");

        // Test basic properties
        assert_eq!(poly.num_points(), 5);
        assert_eq!(poly.points().len(), 5);

        // Test that collocation points are in [-1, 1] (Gauss-Lobatto property)
        for &point in poly.points() {
            assert!((-1.0..=1.0).contains(&point));
        }

        // Test differentiation matrix is square
        let n = poly.num_points();
        let test_vector = nalgebra::DVector::from_fn(n, |i, _| {
            let x = poly.points()[i];
            x * x // Test quadratic function x²
        });

        let derivative = poly.differentiate(&test_vector);
        assert_eq!(derivative.len(), n);

        // For f(x) = x², f'(x) = 2x
        // Check at endpoints where we expect specific values
        let expected_derivative_at_x1: f64 = 2.0 * poly.points()[0]; // f'(-1) = -2
        let computed_derivative_at_x1: f64 = derivative[0];

        // Allow for numerical differentiation error
        assert!((computed_derivative_at_x1 - expected_derivative_at_x1).abs() < 0.1);
    }

    /// Test Fourier transform accuracy with analytical solutions
    #[test]
    fn test_fourier_transform_accuracy() {
        use crate::spectral::FourierTransform;
        use nalgebra::DVector;

        let ft: FourierTransform<f64> =
            FourierTransform::new(8).expect("Failed to create Fourier transform");

        // Test with a simple constant function
        let n = 8;
        let constant_signal = DVector::from_element(n, 1.0);

        let spectrum = ft.forward(&constant_signal).expect("Forward FFT failed");

        // For a constant signal, energy should be concentrated at k=0
        assert!(spectrum.len() == n);

        // The DC component (k=0) for a constant signal = 1 with normalization in forward DFT
        let dc_magnitude = spectrum[0].modulus();

        // This implementation normalizes by dividing by n in forward DFT
        // So for constant signal of 1, DC magnitude should be 1
        let expected_magnitude = 1.0;
        assert!((dc_magnitude - expected_magnitude).abs() < 1e-12);

        // All other components should be essentially zero (relaxed tolerance for naive DFT)
        for i in 1..n {
            assert!(spectrum[i].modulus() < 1e-10);
        }
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
        use cfd_core::domains::mesh_operations::ElementType;

        let mut config: FemConfig<f64> = FemConfig::default();

        // Test different element types
        config.element_type = ElementType::Tetrahedron;
        assert_eq!(config.element_type, ElementType::Tetrahedron);

        config.element_type = ElementType::Hexahedron;
        assert_eq!(config.element_type, ElementType::Hexahedron);

        // Verify stabilization parameters
        assert!(config.tau > 0.0);
        assert!(config.quadrature_order >= 1);
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
        assert!((points[0] - 1.0).abs() < 1e-14);
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

        // For f(x) = x², f''(x) should be approximately 2 everywhere
        for &val in second_derivative.iter() {
            assert!((val - 2.0).abs() < 1e-10); // High accuracy expected
        }
    }
}
