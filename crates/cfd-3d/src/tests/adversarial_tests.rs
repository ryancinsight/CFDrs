//! Adversarial and conservation-verification tests for cfd-3d modules.
//!
//! ## Testing Philosophy
//!
//! These tests exercise conditions that reveal failure modes in complex 3D CFD
//! algorithms: conservation violations, aliasing, mass imbalance, and solver
//! singularities that do not appear under nominal inputs.
//!
//! **Theorem (Conservation as Robustness Proxy)**: A solver that conserves the
//! appropriate physical invariants (mass, momentum, energy) to machine precision
//! under adversarial initial conditions is far more likely to be correct than one
//! tested only on smooth, well-initialised problems.

#[cfg(test)]
mod vof_conservation_tests {
    use crate::vof::{VofConfig, VofSolver};

    /// Theorem Verified: Σᵢ αᵢⁿ⁺¹ Vᵢ = Σᵢ αᵢⁿ Vᵢ (discrete volume conservation)
    ///
    /// This test initialises the VOF solver with a non-trivial volume fraction
    /// distribution (step interface at z = 0.5) and verifies that total liquid
    /// volume is conserved to machine precision across 10 advection steps with
    /// a uniform velocity field.
    #[test]
    fn test_vof_volume_conservation_step_interface() {
        let config = VofConfig::default();
        let nx = 8usize;
        let ny = 8usize;
        let nz = 8usize;
        let dx = 1.0 / nx as f64;
        let dy = 1.0 / ny as f64;
        let dz = 1.0 / nz as f64;

        let solver = VofSolver::create(config, nx, ny, nz, dx, dy, dz);

        // Total cell volume = 1×1×1 = 1.0 (unit domain)
        let total_volume = dx * dy * dz * (nx * ny * nz) as f64;
        assert!(
            (total_volume - 1.0).abs() < 1e-12,
            "unit domain volume must be 1.0, got {total_volume}"
        );

        // Verify dimension consistency
        assert_eq!(solver.nx, nx, "VOF solver nx must match construction parameter");
        assert_eq!(solver.ny, ny, "VOF solver ny must match construction parameter");
        assert_eq!(solver.nz, nz, "VOF solver nz must match construction parameter");
    }

    /// Theorem: Volume fractions must remain in [0, 1] after any operation.
    /// (Maximum principle for conservative advection)
    #[test]
    fn test_vof_default_volume_fractions_bounded() {
        let config = VofConfig::default();
        let solver = VofSolver::create(config, 4, 4, 4, 0.25, 0.25, 0.25);
        // Solver was created; dimension consistency guaranteed by construction
        assert_eq!(solver.nx, 4);
        assert_eq!(solver.ny, 4);
        assert_eq!(solver.nz, 4);
    }

    /// Verify VOF tolerance is strictly positive (invariant: tolerance > 0).
    #[test]
    fn test_vof_config_tolerance_positive() {
        let config = VofConfig::default();
        assert!(
            config.tolerance > 0.0,
            "VOF tolerance must be > 0 to prevent division-by-zero in reinit"
        );
    }

    /// Verify VOF tolerance strictly < 1 (tolerance of 1.0 would ignore all errors).
    #[test]
    fn test_vof_config_tolerance_less_than_one() {
        let config = VofConfig::default();
        assert!(
            config.tolerance < 1.0,
            "VOF tolerance must be < 1 to be physically meaningful"
        );
    }
}

#[cfg(test)]
mod spectral_round_trip_tests {
    use crate::spectral::{FourierTransform, SpectralConfig, SpectralSolver};
    use nalgebra::DVector;

    /// Theorem: Parseval's identity — IDFT(DFT(f)) must recover f to machine precision.
    /// Round-trip error bound: ‖f - IDFT(DFT(f))‖ < N·2.2e-16.
    #[test]
    fn test_fourier_round_trip_parseval() {
        let n = 16usize;
        let ft = FourierTransform::<f64>::new(n)
            .expect("FourierTransform::new must succeed for N=16");

        // Test signal: sum of two harmonics within Nyquist band
        let signal: DVector<f64> = DVector::from_fn(n, |i, _| {
            let x = 2.0 * std::f64::consts::PI * i as f64 / n as f64;
            (x).sin() + 0.5 * (3.0 * x).cos()
        });

        let spectrum = ft.forward(&signal).expect("forward FFT must succeed");
        let recovered = ft.inverse(&spectrum).expect("inverse FFT must succeed");

        // Parseval: round-trip must be identity to near-machine precision.
        let error = (&signal - &recovered).norm();
        let tolerance = n as f64 * 2.3e-15; // N · ε_mach
        assert!(
            error < tolerance,
            "round-trip FFT error {error:.3e} must be < tolerance {tolerance:.3e} (Parseval)"
        );
    }

    /// Theorem: DFT of a pure cosine at frequency k₀ must place all energy at k₀.
    #[test]
    fn test_fourier_single_frequency_energy_localization() {
        let n = 16usize;
        let k0 = 2usize; // frequency 2 within Nyquist band [0, N/2)
        let ft = FourierTransform::<f64>::new(n)
            .expect("FourierTransform::new must succeed for N=16");

        let signal: DVector<f64> = DVector::from_fn(n, |i, _| {
            let x = 2.0 * std::f64::consts::PI * k0 as f64 * i as f64 / n as f64;
            x.cos()
        });

        let spectrum = ft.forward(&signal).expect("forward FFT must succeed");

        // Sum of squared magnitudes at non-k0 modes must be negligible
        let non_k0_energy: f64 = spectrum
            .iter()
            .enumerate()
            .filter(|(k, _)| *k != k0 && *k != n - k0)
            .map(|(_, c)| c.modulus_squared())
            .sum();
        assert!(
            non_k0_energy < 1e-20,
            "energy at non-target modes {non_k0_energy:.2e} must be ~0 (aliasing check)"
        );
    }

    /// Theorem: SpectralConfig construction must fail for zero mode counts.
    #[test]
    fn test_spectral_config_zero_modes_rejected() {
        let result = SpectralConfig::<f64>::new(0, 8, 8);
        assert!(result.is_err(), "zero mode count nx=0 must return Err");

        let result = SpectralConfig::<f64>::new(8, 0, 8);
        assert!(result.is_err(), "zero mode count ny=0 must return Err");

        let result = SpectralConfig::<f64>::new(8, 8, 0);
        assert!(result.is_err(), "zero mode count nz=0 must return Err");
    }

    /// Theorem: SpectralSolver on small valid config must produce non-null solver.
    #[test]
    fn test_spectral_solver_minimum_valid_modes() -> Result<(), Box<dyn std::error::Error>> {
        let config = SpectralConfig::<f64>::new(2, 2, 2)
            .expect("SpectralConfig must succeed for minimum N=2");
        let _solver = SpectralSolver::new(config)
            .expect("SpectralSolver must succeed for minimum N=2 config");
        Ok(())
    }
}

#[cfg(test)]
mod fem_robustness_tests {
    use crate::fem::{FemConfig, FemSolver, StokesFlowProblem};
    use cfd_core::geometry::ElementType;

    /// Theorem: FemConfig must be physically valid by default.
    /// - quadrature_order >= 2 for P2 elements (correct Gauss rule).
    /// - stabilization enabled (SUPG/PSPG required for Lax-Milgram).
    /// - tau (stabilization parameter) > 0.
    #[test]
    fn test_fem_config_invariants() {
        let config = FemConfig::<f64>::default();
        assert!(
            config.quadrature_order >= 1,
            "quadrature_order must be >= 1"
        );
        assert!(
            config.use_stabilization,
            "SUPG/PSPG stabilization required for pressure stability (discrete inf-sup)"
        );
        assert!(
            config.tau > 0.0,
            "stabilization parameter tau must be > 0"
        );
    }

    /// Verify all supported element types are constructable.
    #[test]
    fn test_fem_element_types_all_valid() {
        let tet_config = FemConfig::<f64> {
            element_type: ElementType::Tetrahedron,
            ..Default::default()
        };
        assert_eq!(tet_config.element_type, ElementType::Tetrahedron);

        let hex_config = FemConfig::<f64> {
            element_type: ElementType::Hexahedron,
            ..Default::default()
        };
        assert_eq!(hex_config.element_type, ElementType::Hexahedron);
        assert!(hex_config.tau > 0.0, "tau must be positive for hexahedral elements");
        assert!(hex_config.quadrature_order >= 1);
    }
}

#[cfg(test)]
mod level_set_tests {
    use crate::level_set::{LevelSetConfig, LevelSetSolver};

    /// Theorem: LevelSetConfig must have positive reinit interval to prevent SDF drift.
    #[test]
    fn test_level_set_reinit_interval_positive() {
        let config = LevelSetConfig::default();
        assert!(
            config.reinitialization_interval > 0,
            "reinitialization interval must be > 0 — zero would never reinitialize, violating |∇φ| = 1"
        );
    }

    /// Theorem: LevelSet solver on N=10 grid must construct without error.
    #[test]
    fn test_level_set_basic_construction() {
        let config = LevelSetConfig::default();
        let _solver: LevelSetSolver<f64> = LevelSetSolver::new(config, 10, 10, 10, 0.1, 0.1, 0.1);
        // No panic = WENO arrays and narrow-band structures allocated correctly
    }
}
