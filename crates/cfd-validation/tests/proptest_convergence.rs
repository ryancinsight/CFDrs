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
        for _ in 0..20 {
            monitor.update(error);
            error *= reduction_rate;
        }
        
        let status = monitor.check_status();
        prop_assert!(
            status.is_converged(),
            "Monotonically decreasing sequence should converge"
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
        let mut monitor1 = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 100);
        let mut monitor2 = ConvergenceMonitor::<f64>::new(1e-6, 1e-3, 100);
        
        let mut error1 = 1.0;
        let mut error2 = scale;
        
        for _ in 0..15 {
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
            "Relative convergence should be scale-invariant"
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
    
    // Perfect second-order convergence
    let f1 = 1.0;
    let f2 = 1.0 + 0.04; // Error proportional to h²
    let f3 = 1.0 + 0.16; // Error proportional to (2h)² = 4h²
    
    let gci12 = gci.compute_fine(f1, f2);
    let gci23 = gci.compute_fine(f2, f3);
    
    let ratio = gci23 / (gci12 * 2.0f64.powi(2));
    
    // In asymptotic range, this ratio should be close to 1
    assert!(
        (ratio - 1.0).abs() < 0.1,
        "GCI ratio indicates asymptotic convergence: {}",
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
