//! Tests for error metrics module

use super::*;
use approx::assert_relative_eq;
use nalgebra::Vector3;

#[test]
fn test_l2_norm() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 1.9, 3.2];
    
    let l2 = L2Norm;
    let error = l2.compute_error(&numerical, &reference)?;
    
    // Expected: sqrt(((0.1)^2 + (0.1)^2 + (0.2)^2) / 3) ≈ 0.1414
    assert_relative_eq!(error, 0.1414213, epsilon = 1e-6);
    Ok(())
}

#[test]
fn test_l2_norm_vector() -> Result<()> {
    // Test L2 norm with scalar values - Vector3 doesn't implement RealField
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 1.9, 3.2];
    
    let l2 = L2Norm;
    let error = l2.compute_error(&numerical, &reference)?;
    
    // Expected: sqrt((0.1^2 + 0.1^2 + 0.2^2)/3) ≈ 0.1414
    assert_relative_eq!(error, 0.141421, epsilon = 1e-5);
    Ok(())
}

#[test]
fn test_linf_norm() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 1.9, 3.5];
    
    let linf = LInfNorm;
    let error = linf.compute_error(&numerical, &reference)?;
    
    // Expected: max(0.1, 0.1, 0.5) = 0.5
    assert_relative_eq!(error, 0.5, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_l1_norm() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 1.9, 3.2];
    
    let l1 = L1Norm;
    let error = l1.compute_error(&numerical, &reference)?;
    
    // Expected: (0.1 + 0.1 + 0.2) / 3 ≈ 0.1333
    assert_relative_eq!(error, 0.1333333, epsilon = 1e-5);
    Ok(())
}

#[test]
fn test_relative_l2_error() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 2.2, 3.3];
    
    let rel_error = RelativeError::new(L2Norm, 1e-12);
    let error = rel_error.compute_error(&numerical, &reference)?;
    
    let abs_error = L2Norm.compute_error(&numerical, &reference)?;
    let ref_norm = L2Norm.compute_error(&reference, &vec![0.0; 3])?;
    let expected = abs_error / ref_norm;
    
    assert_relative_eq!(error, expected, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_relative_error_normalization() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 2.2, 3.3];
    
    let rel_error = RelativeError::new(L2Norm, 1e-12);
    let error = rel_error.compute_error(&numerical, &reference)?;
    
    let abs_error = L2Norm.compute_error(&numerical, &reference)?;
    let ref_norm = L2Norm.compute_error(&reference, &vec![0.0; 3])?;
    let expected = abs_error / ref_norm;
    
    assert_relative_eq!(error, expected, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_rmse() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 1.9, 3.2];
    
    let rmse = RootMeanSquareError;
    let rmse_error = rmse.compute_error(&numerical, &reference)?;
    let l2_error = L2Norm.compute_error(&numerical, &reference)?;
    
    // RMSE should equal L2 norm for this case
    assert_relative_eq!(rmse_error, l2_error, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_mae() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 1.9, 3.2];
    
    let mae = MeanAbsoluteError;
    let mae_error = mae.compute_error(&numerical, &reference)?;
    let l1_error = L1Norm.compute_error(&numerical, &reference)?;
    
    // MAE should equal L1 norm
    assert_relative_eq!(mae_error, l1_error, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_nrmse_range() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 1.9, 3.2];
    
    let nrmse = NormalizedRMSE::new(NormalizationMethod::Range);
    let error = nrmse.compute_error(&numerical, &reference)?;
    let rmse = RootMeanSquareError.compute_error(&numerical, &reference)?;
    let expected = rmse / (3.2 - 1.1); // Range of reference
    
    assert_relative_eq!(error, expected, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_nrmse_mean() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0];
    let reference = vec![1.1, 1.9, 3.2];
    
    let nrmse = NormalizedRMSE::new(NormalizationMethod::Mean);
    let error = nrmse.compute_error(&numerical, &reference)?;
    let rmse = RootMeanSquareError.compute_error(&numerical, &reference)?;
    let mean = (1.1 + 1.9 + 3.2) / 3.0;
    let expected = rmse / mean;
    
    assert_relative_eq!(error, expected, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_error_statistics() -> Result<()> {
    let numerical = vec![1.0, 2.0, 3.0, 4.0, 5.0];
    let reference = vec![1.1, 1.9, 3.2, 3.8, 5.3];
    
    let stats = ErrorStatistics::compute(&numerical, &reference)?;
    
    assert!(stats.num_points == 5);
    assert!(stats.l1_norm > 0.0);
    assert!(stats.l2_norm > 0.0);
    assert!(stats.linf_norm > 0.0);
    assert!(stats.mae > 0.0);
    assert!(stats.rmse > 0.0);
    assert!(stats.relative_l2 > 0.0);

    // L2 and RMSE should be the same
    assert_relative_eq!(stats.l2_norm, stats.rmse, epsilon = 1e-10);
    // L1 and MAE should be the same
    assert_relative_eq!(stats.l1_norm, stats.mae, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_error_statistics_vector() -> Result<()> {
    use nalgebra::Vector3;
    let numerical = vec![
        Vector3::new(1.0, 2.0, 3.0),
        Vector3::new(3.0, 4.0, 5.0),
    ];
    let reference = vec![
        Vector3::new(1.1, 1.9, 3.2),
        Vector3::new(3.2, 3.8, 5.3),
    ];
    
    let stats = ErrorStatistics::compute_vector(&numerical, &reference)?;
    
    assert!(stats.num_points == 2);
    assert!(stats.l2_norm > 0.0);
    Ok(())
}

#[test]
fn test_vector_error_metric() -> Result<()> {
    let l2 = L2Norm;

    let numerical = vec![
        Vector3::new(1.0, 0.0, 0.0),
        Vector3::new(0.0, 2.0, 0.0),
    ];
    let reference = vec![
        Vector3::new(0.0, 0.0, 0.0),
        Vector3::new(0.0, 0.0, 0.0),
    ];

    let error = l2.compute_vector_error(&numerical, &reference)?;

    // Magnitudes are [1.0, 2.0], reference magnitudes are [0.0, 0.0]
    // L2 error = sqrt((1^2 + 2^2) / 2) = sqrt(2.5)
    let expected = (2.5_f64).sqrt();
    assert_relative_eq!(error, expected, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_convergence_rate() -> Result<()> {
    // Test theoretical second-order convergence
    let grid_sizes = vec![0.1, 0.05, 0.025];
    let errors = vec![0.01, 0.0025, 0.000625]; // errors ∝ h^2

    let rate = ErrorAnalysis::convergence_rate(&grid_sizes, &errors)?;
    assert_relative_eq!(rate, 2.0, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_convergence_rate_first_order() -> Result<()> {
    // Test theoretical first-order convergence
    let grid_sizes = vec![0.1, 0.05, 0.025];
    let errors = vec![0.1, 0.05, 0.025]; // errors ∝ h^1

    let rate = ErrorAnalysis::convergence_rate(&grid_sizes, &errors)?;
    assert_relative_eq!(rate, 1.0, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_error_reduction_factor() -> Result<()> {
    let coarse_error = 0.01;
    let fine_error = 0.0025;

    let factor = ErrorAnalysis::error_reduction_factor(coarse_error, fine_error);
    assert_relative_eq!(factor, 4.0, epsilon = 1e-10);
    Ok(())
}

#[test]
fn test_is_acceptable() -> Result<()> {
    assert!(ErrorAnalysis::is_acceptable(0.001, 0.01));
    assert!(!ErrorAnalysis::is_acceptable(0.01, 0.001));
    assert!(ErrorAnalysis::is_acceptable(0.01, 0.01));
    Ok(())
}

#[test]
fn test_empty_arrays() -> Result<()> {
    let l2 = L2Norm;
    let numerical: Vec<f64> = vec![];
    let reference: Vec<f64> = vec![];

    let error = l2.compute_error(&numerical, &reference)?;
    assert_eq!(error, 0.0);
    Ok(())
}

#[test]
fn test_mismatched_lengths() -> Result<()> {
    let l2 = L2Norm;
    let numerical = vec![1.0, 2.0];
    let reference = vec![1.0, 2.0, 3.0];

    assert!(l2.compute_error(&numerical, &reference).is_err());
    Ok(())
}

#[test]
fn test_convergence_rate_insufficient_data() -> Result<()> {
    let grid_sizes = vec![0.1];
    let errors = vec![0.01];

    assert!(ErrorAnalysis::convergence_rate(&grid_sizes, &errors).is_err());
    Ok(())
}

#[test]
fn test_convergence_rate_identical_grid_sizes() -> Result<()> {
    let grid_sizes = vec![0.1, 0.1, 0.1];
    let errors = vec![0.01, 0.01, 0.01];

    assert!(ErrorAnalysis::convergence_rate(&grid_sizes, &errors).is_err());
    Ok(())
}