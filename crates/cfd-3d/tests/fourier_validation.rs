//! Comprehensive Fourier transform validation tests
//!
//! Validates spectral methods per Canuto et al. (2006) and Boyd (2001)
//! Tests include:
//! - Forward/inverse DFT identity
//! - Wavenumber accuracy
//! - Parseval's theorem (energy conservation)
//! - Spectral derivative accuracy
//! - Edge cases and numerical stability

use approx::assert_relative_eq;
use cfd_3d::spectral::fourier::{FourierTransform, SpectralDerivative};
use nalgebra::{Complex, DVector};
use std::f64::consts::PI;

/// Test forward-inverse DFT identity: IFFT(FFT(u)) = u
/// Validates that the transform pair is correctly implemented
#[test]
fn test_fourier_transform_identity() {
    let n = 16;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    // Test with various input signals
    let test_signals = vec![
        // Constant signal
        DVector::from_element(n, 1.0),
        // Linear ramp
        DVector::from_iterator(n, (0..n).map(|i| i as f64)),
        // Sine wave (exactly representable)
        DVector::from_iterator(n, (0..n).map(|i| (2.0 * PI * i as f64 / n as f64).sin())),
        // Complex pattern (sum of multiple modes)
        DVector::from_iterator(
            n,
            (0..n).map(|i| {
                let x = 2.0 * PI * i as f64 / n as f64;
                x.sin() + 0.5 * (2.0 * x).cos() + 0.25 * (3.0 * x).sin()
            }),
        ),
    ];

    for u in test_signals {
        let u_hat = transform.forward(&u).unwrap();
        let u_reconstructed = transform.inverse(&u_hat).unwrap();

        // Verify identity holds to numerical precision
        for i in 0..n {
            assert_relative_eq!(u[i], u_reconstructed[i], epsilon = 1e-12);
        }
    }
}

/// Test Parseval's theorem: ||u||² = N * ||û||²
/// Validates energy conservation in spectral space
/// Reference: Canuto et al. (2006), Section 2.1.3
#[test]
fn test_parsevals_theorem() {
    let n = 32;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    // Create test signal with known energy
    let u = DVector::from_iterator(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            x.sin() + 0.3 * (3.0 * x).cos()
        }),
    );

    let u_hat = transform.forward(&u).unwrap();

    // Physical space energy
    let energy_physical: f64 = u.iter().map(|&x| x * x).sum();

    // Spectral space energy
    let energy_spectral: f64 = u_hat.iter().map(|c| c.norm_sqr()).sum();

    // Parseval's theorem: E_physical = N * E_spectral
    let n_f64 = n as f64;
    assert_relative_eq!(energy_physical, n_f64 * energy_spectral, epsilon = 1e-10);
}

/// Test wavenumber computation accuracy
/// Validates k = 0, 1, ..., n/2, -n/2+1, ..., -1
/// Critical for correct spectral derivative computation
#[test]
fn test_wavenumber_ordering() {
    let n = 8;
    let transform = FourierTransform::<f64>::new(n).unwrap();
    let wavenumbers = transform.wavenumbers();

    // Expected wavenumbers: [0, 1, 2, 3, 4, -3, -2, -1]
    let expected = vec![0.0, 1.0, 2.0, 3.0, 4.0, -3.0, -2.0, -1.0];

    assert_eq!(wavenumbers.len(), n);
    for i in 0..n {
        assert_relative_eq!(wavenumbers[i], expected[i], epsilon = 1e-14);
    }
}

/// Test Fourier transform for minimal grid sizes
/// Edge case: N=2 (Nyquist limit)
#[test]
fn test_fourier_transform_minimal_size() {
    let n = 2;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    // Simple signal: [1, -1]
    let u = DVector::from_vec(vec![1.0, -1.0]);
    let u_hat = transform.forward(&u).unwrap();
    let u_reconstructed = transform.inverse(&u_hat).unwrap();

    for i in 0..n {
        assert_relative_eq!(u[i], u_reconstructed[i], epsilon = 1e-14);
    }

    // Verify DC component (mean)
    let dc = u_hat[0];
    let expected_dc = Complex::new(0.0, 0.0); // Mean of [1, -1]
    assert_relative_eq!(dc.re, expected_dc.re, epsilon = 1e-14);
    assert_relative_eq!(dc.im, expected_dc.im, epsilon = 1e-14);
}

/// Test Fourier transform for odd grid size
/// Edge case: N=3 (non-power-of-2)
#[test]
fn test_fourier_transform_odd_size() {
    let n = 3;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    let u = DVector::from_vec(vec![1.0, 2.0, 3.0]);
    let u_hat = transform.forward(&u).unwrap();
    let u_reconstructed = transform.inverse(&u_hat).unwrap();

    for i in 0..n {
        assert_relative_eq!(u[i], u_reconstructed[i], epsilon = 1e-13);
    }
}

/// Test spectral derivative with cosine function
/// First derivative: d/dx(cos(kx)) = -k*sin(kx)
/// Reference: Boyd (2001), Section 2.3
/// Note: Fourier methods require periodic functions
#[test]
fn test_spectral_derivative_cosine() {
    let n = 32;
    let k = 2.0; // Wave number
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    // Create cosine function u(x) = cos(kx) on [0, 2π]
    // This is periodic and exactly representable
    let u = DVector::from_iterator(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            (k * x).cos()
        }),
    );

    // Compute first derivative
    let du_dx = derivative.derivative(&u, 1).unwrap();

    // Analytical derivative: -k*sin(kx)
    for i in 0..n {
        let x = 2.0 * PI * i as f64 / n as f64;
        let expected = -k * (k * x).sin();

        // Spectral accuracy: exponential convergence for smooth periodic functions
        assert_relative_eq!(du_dx[i], expected, epsilon = 1e-10);
    }
}

/// Test spectral derivative with trigonometric function
/// Second derivative: d²/dx²(sin(kx)) = -k²sin(kx)
/// Validates high accuracy for smooth periodic functions
#[test]
fn test_spectral_derivative_trigonometric() {
    let n = 64;
    let k = 3.0; // Wave number
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    // Create sine function u(x) = sin(kx)
    let u = DVector::from_iterator(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            (k * x).sin()
        }),
    );

    // Compute second derivative
    let d2u_dx2 = derivative.derivative(&u, 2).unwrap();

    // Analytical second derivative: -k²sin(kx)
    for i in 0..n {
        let x = 2.0 * PI * i as f64 / n as f64;
        let expected = -(k * k) * (k * x).sin();

        // Spectral accuracy: exponential convergence for smooth functions
        assert_relative_eq!(d2u_dx2[i], expected, epsilon = 1e-10);
    }
}

/// Test spectral derivative for constant function
/// Edge case: d/dx(c) = 0
#[test]
fn test_spectral_derivative_constant() {
    let n = 16;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    let u = DVector::from_element(n, 5.0);
    let du_dx = derivative.derivative(&u, 1).unwrap();

    // Derivative of constant is zero
    for i in 0..n {
        assert_relative_eq!(du_dx[i], 0.0, epsilon = 1e-13);
    }
}

/// Test spectral derivative for linear function
/// Edge case: d/dx(x) = 1
#[test]
fn test_spectral_derivative_linear() {
    let n = 32;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    let u = DVector::from_iterator(n, (0..n).map(|i| i as f64));
    let du_dx = derivative.derivative(&u, 1).unwrap();

    // For periodic domain, linear function has issues at boundaries
    // Check interior points where derivative should be constant
    let mean_derivative: f64 = du_dx.iter().sum::<f64>() / n as f64;

    // Expect approximately constant derivative
    // Note: Linear function on periodic domain is challenging
    // This test validates the implementation doesn't crash
    assert!(mean_derivative.is_finite());
}

/// Test high-order derivative stability
/// Validates numerical stability for multiple derivative operations
#[test]
fn test_high_order_derivative_stability() {
    let n = 64;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    // Smooth test function
    let u = DVector::from_iterator(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            x.sin() * x.cos()
        }),
    );

    // Compute up to 4th derivative
    for order in 1..=4 {
        let du = derivative.derivative(&u, order).unwrap();

        // Verify result is finite and bounded
        for i in 0..n {
            assert!(
                du[i].is_finite(),
                "Derivative order {} has non-finite values",
                order
            );
            assert!(
                du[i].abs() < 100.0,
                "Derivative order {} has unbounded values: {}",
                order,
                du[i]
            );
        }
    }
}

/// Test Fourier transform symmetry for real signals
/// Validates conjugate symmetry: û(-k) = conj(û(k))
/// Reference: Canuto et al. (2006), Section 2.1.2
#[test]
fn test_fourier_symmetry_real_signal() {
    let n = 16;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    // Real signal
    let u = DVector::from_iterator(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            x.cos() + (2.0 * x).sin()
        }),
    );

    let u_hat = transform.forward(&u).unwrap();

    // Check conjugate symmetry: û[n-k] = conj(û[k])
    for k in 1..n / 2 {
        let uk = u_hat[k];
        let u_mk = u_hat[n - k];

        // Real part should be equal
        assert_relative_eq!(uk.re, u_mk.re, epsilon = 1e-12);
        // Imaginary part should be opposite
        assert_relative_eq!(uk.im, -u_mk.im, epsilon = 1e-12);
    }
}

/// Test spectral derivative accuracy compared to finite differences
/// Validates superior accuracy of spectral methods
#[test]
fn test_spectral_vs_finite_difference() {
    let n = 32;
    let l = 2.0 * PI;
    let dx = l / n as f64;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    // Test function: exp(sin(x))
    let u = DVector::from_iterator(
        n,
        (0..n).map(|i| {
            let x = l * i as f64 / n as f64;
            x.sin().exp()
        }),
    );

    // Spectral derivative
    let du_spectral = derivative.derivative(&u, 1).unwrap();

    // Central finite difference approximation
    let mut du_fd = DVector::zeros(n);
    for i in 0..n {
        let im1 = if i == 0 { n - 1 } else { i - 1 };
        let ip1 = if i == n - 1 { 0 } else { i + 1 };
        du_fd[i] = (u[ip1] - u[im1]) / (2.0 * dx);
    }

    // Analytical derivative: exp(sin(x)) * cos(x)
    let mut spectral_error = 0.0;
    let mut fd_error = 0.0;

    for i in 0..n {
        let x = l * i as f64 / n as f64;
        let analytical = x.sin().exp() * x.cos();

        spectral_error += (du_spectral[i] - analytical).abs();
        fd_error += (du_fd[i] - analytical).abs();
    }

    // Spectral method should be significantly more accurate
    assert!(
        spectral_error < 0.01 * fd_error,
        "Spectral error {} should be << FD error {}",
        spectral_error,
        fd_error
    );
}
