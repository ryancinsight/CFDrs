//! Comprehensive Fourier transform validation tests
//!
//! Validates spectral methods per Canuto et al. (2006) and Boyd (2001)
//! using the Leto/Eunomia public Fourier seam:
//! - Forward/inverse DFT identity
//! - Wavenumber accuracy
//! - Parseval's theorem (energy conservation)
//! - Spectral derivative accuracy
//! - Edge cases and numerical stability

use eunomia::assert_relative_eq;
use cfd_3d::spectral::fourier::{FourierTransform, SpectralDerivative};
use eunomia::Complex;
use leto::Array1;
use std::f64::consts::PI;

fn array_from_iter<T>(n: usize, values: impl IntoIterator<Item = T>) -> Array1<T> {
    Array1::from_vec([n], values.into_iter().collect())
        .expect("invariant: test vector length matches declared Leto shape")
}

fn array_from_vec<T>(values: Vec<T>) -> Array1<T> {
    Array1::from_vec([values.len()], values)
        .expect("invariant: test vector length matches declared Leto shape")
}

fn value<T: Copy>(array: &Array1<T>, index: usize) -> T {
    *array
        .get([index])
        .expect("invariant: generated test index is in bounds")
}

fn set_value<T>(array: &mut Array1<T>, index: usize, value: T) {
    *array
        .get_mut([index])
        .expect("invariant: generated test index is in bounds") = value;
}

/// Test forward-inverse DFT identity: IFFT(FFT(u)) = u
/// Validates that the transform pair is correctly implemented.
#[test]
fn test_fourier_transform_identity() {
    let n = 16;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    let test_signals = vec![
        Array1::from_elem([n], 1.0),
        array_from_iter(n, (0..n).map(|i| i as f64)),
        array_from_iter(n, (0..n).map(|i| (2.0 * PI * i as f64 / n as f64).sin())),
        array_from_iter(
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

        for i in 0..n {
            assert_relative_eq!(value(&u, i), value(&u_reconstructed, i), epsilon = 1e-12);
        }
    }
}

/// Test Parseval's theorem: ||u||² = N * ||û||².
/// Reference: Canuto et al. (2006), Section 2.1.3.
#[test]
fn test_parsevals_theorem() {
    let n = 32;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    let u = array_from_iter(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            x.sin() + 0.3 * (3.0 * x).cos()
        }),
    );

    let u_hat = transform.forward(&u).unwrap();

    let energy_physical: f64 = u.iter().map(|&x| x * x).sum();
    let energy_spectral: f64 = u_hat.iter().map(|c| c.norm_sqr()).sum();

    assert_relative_eq!(energy_physical, n as f64 * energy_spectral, epsilon = 1e-10);
}

/// Test wavenumber computation accuracy.
#[test]
fn test_wavenumber_ordering() {
    let n = 8;
    let transform = FourierTransform::<f64>::new(n).unwrap();
    let wavenumbers = transform.wavenumbers();
    let expected = [0.0, 1.0, 2.0, 3.0, 4.0, -3.0, -2.0, -1.0];

    assert_eq!(wavenumbers.len(), n);
    for i in 0..n {
        assert_relative_eq!(wavenumbers[i], expected[i], epsilon = 1e-14);
    }
}

/// Test Fourier transform for minimal grid sizes.
#[test]
fn test_fourier_transform_minimal_size() {
    let n = 2;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    let u = array_from_vec(vec![1.0, -1.0]);
    let u_hat = transform.forward(&u).unwrap();
    let u_reconstructed = transform.inverse(&u_hat).unwrap();

    for i in 0..n {
        assert_relative_eq!(value(&u, i), value(&u_reconstructed, i), epsilon = 1e-14);
    }

    let dc = value(&u_hat, 0);
    let expected_dc = Complex::new(0.0, 0.0);
    assert_relative_eq!(dc.re, expected_dc.re, epsilon = 1e-14);
    assert_relative_eq!(dc.im, expected_dc.im, epsilon = 1e-14);
}

/// Test Fourier transform for odd grid size.
#[test]
fn test_fourier_transform_odd_size() {
    let n = 3;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    let u = array_from_vec(vec![1.0, 2.0, 3.0]);
    let u_hat = transform.forward(&u).unwrap();
    let u_reconstructed = transform.inverse(&u_hat).unwrap();

    for i in 0..n {
        assert_relative_eq!(value(&u, i), value(&u_reconstructed, i), epsilon = 1e-13);
    }
}

/// Test spectral derivative with cosine function.
#[test]
fn test_spectral_derivative_cosine() {
    let n = 32;
    let k = 2.0;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    let u = array_from_iter(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            (k * x).cos()
        }),
    );

    let du_dx = derivative.derivative(&u, 1).unwrap();

    for i in 0..n {
        let x = 2.0 * PI * i as f64 / n as f64;
        let expected = -k * (k * x).sin();
        assert_relative_eq!(value(&du_dx, i), expected, epsilon = 1e-10);
    }
}

/// Test spectral derivative with trigonometric function.
#[test]
fn test_spectral_derivative_trigonometric() {
    let n = 64;
    let k = 3.0;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    let u = array_from_iter(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            (k * x).sin()
        }),
    );

    let d2u_dx2 = derivative.derivative(&u, 2).unwrap();

    for i in 0..n {
        let x = 2.0 * PI * i as f64 / n as f64;
        let expected = -(k * k) * (k * x).sin();
        assert_relative_eq!(value(&d2u_dx2, i), expected, epsilon = 1e-10);
    }
}

/// Test spectral derivative for constant function.
#[test]
fn test_spectral_derivative_constant() {
    let n = 16;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    let u = Array1::from_elem([n], 5.0);
    let du_dx = derivative.derivative(&u, 1).unwrap();

    for i in 0..n {
        assert_relative_eq!(value(&du_dx, i), 0.0, epsilon = 1e-13);
    }
}

/// Test spectral derivative for linear function.
#[test]
fn test_spectral_derivative_linear() {
    let n = 32;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    let u = array_from_iter(n, (0..n).map(|i| i as f64));
    let du_dx = derivative.derivative(&u, 1).unwrap();
    let mean_derivative: f64 = du_dx.iter().sum::<f64>() / n as f64;

    assert!(mean_derivative.is_finite());
}

/// Test high-order derivative stability.
#[test]
fn test_high_order_derivative_stability() {
    let n = 64;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    let u = array_from_iter(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            x.sin() * x.cos()
        }),
    );

    for order in 1..=4 {
        let du = derivative.derivative(&u, order).unwrap();

        for i in 0..n {
            let derivative_value = value(&du, i);
            assert!(
                derivative_value.is_finite(),
                "Derivative order {} has non-finite values",
                order
            );
            assert!(
                derivative_value.abs() < 100.0,
                "Derivative order {} has unbounded values: {}",
                order,
                derivative_value
            );
        }
    }
}

/// Test Fourier transform symmetry for real signals.
#[test]
fn test_fourier_symmetry_real_signal() {
    let n = 16;
    let transform = FourierTransform::<f64>::new(n).unwrap();

    let u = array_from_iter(
        n,
        (0..n).map(|i| {
            let x = 2.0 * PI * i as f64 / n as f64;
            x.cos() + (2.0 * x).sin()
        }),
    );

    let u_hat = transform.forward(&u).unwrap();

    for k in 1..n / 2 {
        let uk = value(&u_hat, k);
        let u_mk = value(&u_hat, n - k);

        assert_relative_eq!(uk.re, u_mk.re, epsilon = 1e-12);
        assert_relative_eq!(uk.im, -u_mk.im, epsilon = 1e-12);
    }
}

/// Test spectral derivative accuracy compared to finite differences.
#[test]
fn test_spectral_vs_finite_difference() {
    let n = 32;
    let l = 2.0 * PI;
    let dx = l / n as f64;
    let derivative = SpectralDerivative::<f64>::new(n).unwrap();

    let u = array_from_iter(
        n,
        (0..n).map(|i| {
            let x = l * i as f64 / n as f64;
            x.sin().exp()
        }),
    );

    let du_spectral = derivative.derivative(&u, 1).unwrap();

    let mut du_fd = Array1::zeros([n]);
    for i in 0..n {
        let im1 = if i == 0 { n - 1 } else { i - 1 };
        let ip1 = if i == n - 1 { 0 } else { i + 1 };
        set_value(
            &mut du_fd,
            i,
            (value(&u, ip1) - value(&u, im1)) / (2.0 * dx),
        );
    }

    let mut spectral_error = 0.0;
    let mut fd_error = 0.0;

    for i in 0..n {
        let x = l * i as f64 / n as f64;
        let analytical = x.sin().exp() * x.cos();

        spectral_error += (value(&du_spectral, i) - analytical).abs();
        fd_error += (value(&du_fd, i) - analytical).abs();
    }

    assert!(
        spectral_error < 0.01 * fd_error,
        "Spectral error {} should be << FD error {}",
        spectral_error,
        fd_error
    );
}
