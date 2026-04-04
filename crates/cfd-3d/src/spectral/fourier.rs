//! Fourier transform operations for spectral methods.
//!
//! This module delegates plan management and transform execution to the
//! tracked Apollo FFT workspace submodule. Apollo owns the reusable CPU FFT
//! cache and backend logic; this wrapper preserves the legacy `DVector`
//! interface and scaling convention used by the existing CFDrs spectral tests.
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

use apollofft::{fft_1d_array, ifft_1d_array, Complex64 as ApolloComplex64};
use cfd_core::error::Result;
use nalgebra::{Complex, DVector, RealField};
use ndarray::Array1;
use num_traits::{FromPrimitive, ToPrimitive};

fn invalid_configuration(message: impl Into<String>) -> cfd_core::error::Error {
    cfd_core::error::Error::InvalidConfiguration(message.into())
}

fn convert_to_f64<T>(value: T, context: &str) -> Result<f64>
where
    T: ToPrimitive,
{
    value.to_f64().ok_or_else(|| {
        invalid_configuration(format!("{context}: value cannot be represented as f64"))
    })
}

fn convert_from_f64<T>(value: f64, context: &str) -> Result<T>
where
    T: FromPrimitive,
{
    T::from_f64(value).ok_or_else(|| {
        invalid_configuration(format!("{context}: value cannot be represented in the target scalar type"))
    })
}

fn scalar_from_usize<T>(value: usize, context: &str) -> Result<T>
where
    T: FromPrimitive,
{
    T::from_usize(value).ok_or_else(|| {
        invalid_configuration(format!("{context}: size {value} cannot be represented in the target scalar type"))
    })
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
pub struct FourierTransform<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + ToPrimitive + Copy> {
    /// Number of modes
    n: usize,
    /// Wavenumbers
    wavenumbers: Vec<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + ToPrimitive + Copy>
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

            let wn = <T as FromPrimitive>::from_f64(k).ok_or_else(|| {
                cfd_core::error::Error::InvalidConfiguration("Cannot convert wavenumber".into())
            })?;
            wavenumbers.push(wn);
        }

        Ok(Self { n, wavenumbers })
    }

    /// Forward Fast Fourier Transform (FFT) using Cooley-Tukey algorithm
    pub fn forward(&self, u: &DVector<T>) -> Result<DVector<Complex<T>>> {
        ensure_length(u.len(), self.n, "FourierTransform::forward")?;

        let scale = scalar_from_usize::<T>(self.n, "FourierTransform::forward scale")?;
        let real_signal = Array1::from_vec(
            u.iter()
                .copied()
                .map(|value| convert_to_f64(value, "FourierTransform::forward input"))
                .collect::<Result<Vec<_>>>()?,
        );

        let spectrum = fft_1d_array(&real_signal);
        let normalized = spectrum
            .into_iter()
            .map(|value| {
                Ok(Complex::new(
                    convert_from_f64::<T>(value.re, "FourierTransform::forward output real part")?
                        / scale,
                    convert_from_f64::<T>(value.im, "FourierTransform::forward output imaginary part")?
                        / scale,
                ))
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(DVector::from_vec(normalized))
    }

    /// Inverse Fast Fourier Transform (IFFT)
    pub fn inverse(&self, u_hat: &DVector<Complex<T>>) -> Result<DVector<T>> {
        ensure_length(u_hat.len(), self.n, "FourierTransform::inverse")?;

        let scale = scalar_from_usize::<T>(self.n, "FourierTransform::inverse scale")?;
        let spectrum = Array1::from_vec(
            u_hat
                .iter()
                .copied()
                .map(|value| {
                    Ok(ApolloComplex64::new(
                        convert_to_f64(value.re, "FourierTransform::inverse input real part")?,
                        convert_to_f64(value.im, "FourierTransform::inverse input imaginary part")?,
                    ))
                })
                .collect::<Result<Vec<_>>>()?,
        );

        let spatial = ifft_1d_array(&spectrum);
        let recovered = spatial
            .into_iter()
            .map(|value| {
                Ok(convert_from_f64::<T>(
                    value,
                    "FourierTransform::inverse output",
                )? * scale)
            })
            .collect::<Result<Vec<_>>>()?;

        Ok(DVector::from_vec(recovered))
    }

    /// Get wavenumbers
    #[must_use]
    pub fn wavenumbers(&self) -> &[T] {
        &self.wavenumbers
    }
}

/// Spectral derivative computation
pub struct SpectralDerivative<
    T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + ToPrimitive + Copy,
> {
    transform: FourierTransform<T>,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + FromPrimitive + ToPrimitive + Copy>
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
    pub fn derivative(&self, u: &DVector<T>, order: usize) -> Result<DVector<T>> {
        // Transform to spectral space
        let u_hat = self.transform.forward(u)?;

        // Multiply by (ik)^order
        let mut du_hat = DVector::zeros(self.transform.n);
        let i = Complex::new(T::zero(), T::one());

        for (k_idx, &k) in self.transform.wavenumbers().iter().enumerate() {
            let ik = i * k;
            let factor = ik.powi(order as i32);
            du_hat[k_idx] = factor * u_hat[k_idx];
        }

        // Transform back to physical space
        self.transform.inverse(&du_hat)
    }
}
