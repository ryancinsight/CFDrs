//! Property-based tests for convergence monitoring and criteria
//!
//! Uses proptest for generative testing of convergence behavior under various conditions.

use cfd_validation::convergence::{
    ConvergenceMonitor, ConvergenceStatus, GridConvergenceIndex,
};
use proptest::prelude::*;

proptest! {
    /// Test that convergence monitor correctly identifies converging sequences
    #[test]
    fn test_converging_sequence_detected(
        base_error in 0.1f64..10.0,
        reduction_rate in 0.1f64..0.9
    ) {
        let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 100);
        
        let mut error = base_error;
        // Run enough iterations to ensure we reach absolute tolerance
        // For worst case (base=10, rate=0.9), we need: 10 * 0.9^n < 1e-6
        // => n > log(1e-7) / log(0.9) ≈ 152 iterations
        for _ in 0..200 {
            monitor.update(error);
            error *= reduction_rate;
        }
        
        let status = monitor.check_status();
        prop_assert!(
            status.is_converged(),
            "Monotonically decreasing sequence should converge after sufficient iterations"
        );
    }

    /// Test that diverging sequences are detected
    #[test]
    fn test_diverging_sequence_detected(
        base_error in 0.01f64..1.0,
        growth_rate in 1.2f64..2.0
    ) {
        let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 100);
        
        let mut error = base_error;
        for _ in 0..10 {
            monitor.update(error);
            error *= growth_rate;
        }
        
        let status = monitor.check_status();
        prop_assert!(
            matches!(status, ConvergenceStatus::Diverging { .. }),
            "Growing error sequence should be detected as diverging"
        );
    }

    /// Test that stalled convergence is detected
    #[test]
    fn test_stalled_convergence_detected(
        stall_error in 0.01f64..1.0
    ) {
        let mut monitor = ConvergenceMonitor::<f64>::new(1e-7, 1e-3, 10);
        
        // Simulate stalling at a particular error level
        for _ in 0..15 {
            monitor.update(stall_error * (1.0 + 1e-8));
        }
        
        let status = monitor.check_status();
        prop_assert!(
            matches!(status, ConvergenceStatus::Stalled { .. }),
            "Unchanging error should be detected as stalled"
        );
    }

    /// Test GCI calculation produces positive uncertainty estimates
    #[test]
    fn test_gci_positive_uncertainty(
        order in 1.0f64..4.0,
        refinement_ratio in 1.5f64..4.0,
        f_fine in 0.1f64..100.0,
        relative_change in 0.001f64..0.5
    ) {
        let gci = GridConvergenceIndex::<f64>::new(3, order, refinement_ratio);
        let f_coarse = f_fine * (1.0 + relative_change);
        
        let gci_fine = gci.compute_fine(f_fine, f_coarse);
        
        prop_assert!(
            gci_fine >= 0.0,
            "GCI should produce non-negative uncertainty: got {}", gci_fine
        );
        prop_assert!(
            gci_fine.is_finite(),
            "GCI should produce finite values"
        );
    }

    /// Test that relative convergence criterion is scale-invariant
    #[test]
    fn test_relative_convergence_scale_invariant(
        scale in 0.01f64..1000.0,
        reduction_rate in 0.1f64..0.9
    ) {
        // Use scale-dependent absolute tolerance to ensure scale invariance
        let abs_tol_1 = 1e-6;
        let abs_tol_2 = 1e-6 * scale;
        let rel_tol = 1e-3;
        
        let mut monitor1 = ConvergenceMonitor::<f64>::new(abs_tol_1, rel_tol, 100);
        let mut monitor2 = ConvergenceMonitor::<f64>::new(abs_tol_2, rel_tol, 100);
        
        let mut error1 = 1.0;
        let mut error2 = scale;
        
        // Run enough iterations to reach convergence
        for _ in 0..200 {
            monitor1.update(error1);
            monitor2.update(error2);
            error1 *= reduction_rate;
            error2 *= reduction_rate;
        }
        
        let status1 = monitor1.check_status();
        let status2 = monitor2.check_status();
        
        prop_assert_eq!(
            status1.is_converged(),
            status2.is_converged(),
            "Relative convergence should be scale-invariant with scale-dependent tolerances"
        );
    }

    /// Test convergence monitor behavior with noisy convergence
    #[test]
    fn test_noisy_convergence_robustness(
        base_error in 0.1f64..10.0,
        reduction_rate in 0.5f64..0.95,
        noise_amplitude in 0.0f64..0.1
    ) {
        let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 100);
        
        let mut error = base_error;
        for i in 0..30 {
            // Add oscillating noise to convergence
            let noise = noise_amplitude * (i as f64 * 0.5).sin();
            let noisy_error = error * (1.0 + noise);
            monitor.update(noisy_error);
            error *= reduction_rate;
        }
        
        let status = monitor.check_status();
        // With sufficient iterations and strong reduction, should eventually converge
        // even with noise
        prop_assert!(
            status.is_converged() || matches!(status, ConvergenceStatus::NotConverged { .. }),
            "Noisy but converging sequence should eventually converge or be in progress"
        );
    }
}

/// Test that the GCI asymptotic range indicator works correctly
#[test]
fn test_gci_asymptotic_range() {
    let gci = GridConvergenceIndex::<f64>::new(3, 2.0, 2.0);
    
    // For true asymptotic second-order convergence, we need:
    // f(h) = f_exact + C*h²
    // But GCI uses relative error: epsilon = |f_coarse - f_fine| / |f_fine|
    // 
    // Let's use a large enough base value to minimize relative error effects
    let f_exact = 100.0;
    let c = 1.0;
    let h1 = 0.1;
    let h2 = 0.2;
    let h3 = 0.4;
    
    let f1 = f_exact + c * h1 * h1; // 100.01
    let f2 = f_exact + c * h2 * h2; // 100.04
    let f3 = f_exact + c * h3 * h3; // 100.16
    
    // Compute GCI between grids 1-2 (fine grid GCI)
    let gci_fine_12 = gci.compute_fine(f1, f2);
    // Compute GCI between grids 2-3 (fine grid GCI for coarser pair)
    let gci_fine_23 = gci.compute_fine(f2, f3);
    
    // According to Roache (1998), for asymptotic convergence:
    // GCI_23 / (r^p * GCI_12) ≈ 1
    let r_p = 2.0f64.powi(2); // refinement_ratio^order = 2^2 = 4
    let ratio = gci_fine_23 / (r_p * gci_fine_12);
    
    // In asymptotic range, this ratio should be close to 1
    // With larger base values, relative error effects are minimized
    assert!(
        (ratio - 1.0).abs() < 0.02,
        "GCI ratio indicates asymptotic convergence: {} (expected ~1.0)",
        ratio
    );
}

/// Test convergence history tracking
#[test]
fn test_convergence_history() {
    let mut monitor = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 100);
    
    let errors = vec![1.0, 0.5, 0.25, 0.125, 0.0625];
    for &error in &errors {
        monitor.update(error);
    }
    
    assert_eq!(monitor.history.len(), errors.len());
    
    for (i, &error) in errors.iter().enumerate() {
        assert_eq!(monitor.history[i], error);
    }
}
