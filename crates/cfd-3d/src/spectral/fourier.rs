//! Fourier transform operations for spectral methods.
//!
//! This module delegates plan management and transform execution to the
//! tracked Apollo FFT workspace submodule. Apollo owns the reusable CPU FFT
//! cache and backend logic; this wrapper exposes Leto arrays and Eunomia
//! complex scalars at the CFDrs spectral boundary.
//!
//! # Theorem — Parseval Scaling Under Legacy CFDrs Normalization
//!
//! Apollo's 1D CPU plans use an unnormalized forward transform and a normalized
//! inverse. CFDrs keeps the historical convention where the stored spectrum is
//! scaled by $1/n$ and the inverse rescales by $n$, so that round-tripping
//! still satisfies `inverse(forward(u)) = u` and the Parseval check in the test
//! suite remains `‖u‖² = n‖\hat{u}‖²`.
//!
//! # Theorem — 3/2 Dealiasing Rule (Orszag 1971)
//!
//! For a quadratic nonlinearity computed on $N$ modes, the aliased
//! contribution can be eliminated by zero-padding to $M \geq 3N/2$ modes
//! before transforming to physical space.

use crate::atlas_array::{fft_1d_array, ifft_1d_array, values};
use apollo_fft::Complex64 as ApolloComplex64;
use cfd_core::error::Result;
use core::ops::Neg;
use eunomia::{Complex, FloatElement, NumericElement};
use leto::Array1;

fn invalid_configuration(message: impl Into<String>) -> cfd_core::error::Error {
    cfd_core::error::Error::InvalidConfiguration(message.into())
}

fn convert_to_f64<T>(value: T, _context: &str) -> Result<f64>
where
    T: NumericElement,
{
    Ok(<T as NumericElement>::to_f64(value))
}

fn convert_from_f64<T>(value: f64, _context: &str) -> Result<T>
where
    T: FloatElement,
{
    Ok(<T as FloatElement>::from_f64(value))
}

fn scalar_from_usize<T>(value: usize, context: &str) -> Result<T>
where
    T: FloatElement,
{
    convert_from_f64(value as f64, context)
}

fn ensure_length(actual: usize, expected: usize, context: &str) -> Result<()> {
    if actual == expected {
        Ok(())
    } else {
        Err(invalid_configuration(format!(
            "{context}: expected length {expected}, got {actual}"
        )))
    }
}

/// Fourier transform operations
pub struct FourierTransform<
    T: cfd_mesh::domain::core::Scalar + FloatElement + Copy + Neg<Output = T>,
> {
    /// Number of modes
    n: usize,
    /// Wavenumbers
    wavenumbers: Vec<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + FloatElement + Copy + Neg<Output = T>>
    FourierTransform<T>
{
    /// Create new Fourier transform operator
    pub fn new(n: usize) -> Result<Self> {
        if n == 0 {
            return Err(invalid_configuration(
                "FourierTransform: number of modes must be greater than zero",
            ));
        }

        let mut wavenumbers = Vec::with_capacity(n);

        // Compute wavenumbers k = 0, 1, ..., n/2, -n/2+1, ..., -1
        for i in 0..n {
            let k = if i <= n / 2 {
                i as f64
            } else {
                (i as f64) - (n as f64)
            };

            let wn = convert_from_f64(k, "FourierTransform::new wavenumber")?;
            wavenumbers.push(wn);
        }

        Ok(Self { n, wavenumbers })
    }

    /// Forward Fast Fourier Transform (FFT) using Cooley-Tukey algorithm
    pub fn forward(&self, u: &Array1<T>) -> Result<Array1<Complex<T>>> {
        ensure_length(u.size(), self.n, "FourierTransform::forward")?;

        let scale = scalar_from_usize::<T>(self.n, "FourierTransform::forward scale")?;
        let real_values = u
            .iter()
            .copied()
            .map(|value| convert_to_f64(value, "FourierTransform::forward input"))
            .collect::<Result<Vec<_>>>()?;
        let real_signal = Array1::from_vec([self.n], real_values).map_err(|error| {
            invalid_configuration(format!("FourierTransform::forward input shape: {error}"))
        })?;

        let spectrum = fft_1d_array(&real_signal)?;
        let normalized = values(&spectrum)
            .into_iter()
            .map(|value| {
                Ok(Complex::new(
                    convert_from_f64::<T>(value.re, "FourierTransform::forward output real part")?
                        / scale,
                    convert_from_f64::<T>(
                        value.im,
                        "FourierTransform::forward output imaginary part",
                    )? / scale,
                ))
            })
            .collect::<Result<Vec<_>>>()?;

        Array1::from_vec([self.n], normalized).map_err(|error| {
            invalid_configuration(format!("FourierTransform::forward output shape: {error}"))
        })
    }

    /// Inverse Fast Fourier Transform (IFFT)
    pub fn inverse(&self, u_hat: &Array1<Complex<T>>) -> Result<Array1<T>> {
        ensure_length(u_hat.size(), self.n, "FourierTransform::inverse")?;

        let scale = scalar_from_usize::<T>(self.n, "FourierTransform::inverse scale")?;
        let spectral_values = u_hat
            .iter()
            .copied()
            .map(|value| {
                Ok(ApolloComplex64::new(
                    convert_to_f64(value.re, "FourierTransform::inverse input real part")?,
                    convert_to_f64(value.im, "FourierTransform::inverse input imaginary part")?,
                ))
            })
            .collect::<Result<Vec<_>>>()?;
        let spectrum = Array1::from_vec([self.n], spectral_values).map_err(|error| {
            invalid_configuration(format!("FourierTransform::inverse input shape: {error}"))
        })?;

        let spatial = ifft_1d_array(&spectrum)?;
        let recovered = values(&spatial)
            .into_iter()
            .map(|value| {
                Ok(convert_from_f64::<T>(value, "FourierTransform::inverse output")? * scale)
            })
            .collect::<Result<Vec<_>>>()?;

        Array1::from_vec([self.n], recovered).map_err(|error| {
            invalid_configuration(format!("FourierTransform::inverse output shape: {error}"))
        })
    }

    /// Get wavenumbers
    #[must_use]
    pub fn wavenumbers(&self) -> &[T] {
        &self.wavenumbers
    }
}

/// Spectral derivative computation
pub struct SpectralDerivative<
    T: cfd_mesh::domain::core::Scalar + FloatElement + Copy + Neg<Output = T>,
> {
    transform: FourierTransform<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + FloatElement + Copy + Neg<Output = T>>
    SpectralDerivative<T>
{
    /// Create a new spectral derivative operator for the given grid size
    ///
    /// # Arguments
    /// * `n` - Number of grid points (any positive size is supported by Apollo)
    ///
    /// # Returns
    /// * `Result<Self>` - New spectral derivative operator or error if invalid size
    pub fn new(n: usize) -> Result<Self> {
        Ok(Self {
            transform: FourierTransform::new(n)?,
        })
    }

    /// Compute spectral derivative
    /// d^n u / dx^n = IFFT(i*k)^n * FFT(u)
    pub fn derivative(&self, u: &Array1<T>, order: usize) -> Result<Array1<T>> {
        // Transform to spectral space
        let u_hat = self.transform.forward(u)?;

        // Multiply by (ik)^order
        let mut du_hat = Array1::zeros([self.transform.n]);
        let i = Complex::new(<T as NumericElement>::ZERO, <T as NumericElement>::ONE);

        for (k_idx, &k) in self.transform.wavenumbers().iter().enumerate() {
            let ik = i * k;
            let factor = ik.powi(order as i32);
            du_hat[[k_idx]] = factor * u_hat[[k_idx]];
        }

        // Transform back to physical space
        self.transform.inverse(&du_hat)
    }
}
