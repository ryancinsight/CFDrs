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

use crate::atlas_array::values;
use apollo_fft::{fft_1d_array, ifft_1d_array, PlanCacheProvider, PlanScratch, RealFftData};
use cfd_core::error::Result;
use core::ops::Neg;
use eunomia::{Complex, FloatElement, NumericElement};
use leto::Array1;

fn invalid_configuration(message: impl Into<String>) -> cfd_core::error::Error {
    cfd_core::error::Error::InvalidConfiguration(message.into())
}

fn scalar_from_usize<T>(value: usize) -> T
where
    T: FloatElement,
{
    <T as FloatElement>::from_f64(value as f64)
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
    T: cfd_mesh::domain::core::Scalar
        + FloatElement
        + RealFftData<PlanScalar = T>
        + PlanCacheProvider
        + Copy
        + Neg<Output = T>,
> {
    /// Number of modes
    n: usize,
    /// Wavenumbers
    wavenumbers: Vec<T>,
}

impl<T> FourierTransform<T>
where
    T: cfd_mesh::domain::core::Scalar
        + FloatElement
        + RealFftData<PlanScalar = T>
        + PlanCacheProvider
        + Copy
        + Neg<Output = T>,
    Complex<T>: Copy + PlanScratch,
{
    /// Create new Fourier transform operator
    pub fn new(n: usize) -> Result<Self> {
        if n == 0 {
            return Err(invalid_configuration(
                "FourierTransform: number of modes must be greater than zero",
            ));
        }

        let mut wavenumbers = Vec::with_capacity(n);

        let n_scalar = scalar_from_usize::<T>(n);

        // Compute wavenumbers k = 0, 1, ..., n/2, -n/2+1, ..., -1 in T.
        for i in 0..n {
            let index = scalar_from_usize::<T>(i);
            let k = if i <= n / 2 { index } else { index - n_scalar };

            wavenumbers.push(k);
        }

        Ok(Self { n, wavenumbers })
    }

    /// Forward Fast Fourier Transform (FFT) using Cooley-Tukey algorithm
    pub fn forward(&self, u: &Array1<T>) -> Result<Array1<Complex<T>>> {
        ensure_length(u.size(), self.n, "FourierTransform::forward")?;

        let scale = scalar_from_usize::<T>(self.n);
        let spectrum = fft_1d_array::<T>(u);
        let normalized = values(&spectrum)
            .into_iter()
            .map(|value| Complex::new(value.re / scale, value.im / scale))
            .collect::<Vec<_>>();

        Array1::from_vec([self.n], normalized).map_err(|error| {
            invalid_configuration(format!("FourierTransform::forward output shape: {error}"))
        })
    }

    /// Inverse Fast Fourier Transform (IFFT)
    pub fn inverse(&self, u_hat: &Array1<Complex<T>>) -> Result<Array1<T>> {
        ensure_length(u_hat.size(), self.n, "FourierTransform::inverse")?;

        let scale = scalar_from_usize::<T>(self.n);
        let spectrum =
            Array1::from_vec([self.n], u_hat.iter().copied().collect()).map_err(|error| {
                invalid_configuration(format!("FourierTransform::inverse input shape: {error}"))
            })?;

        let spatial = ifft_1d_array::<T>(&spectrum);
        let recovered = values(&spatial)
            .into_iter()
            .map(|value| value * scale)
            .collect::<Vec<_>>();

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
    T: cfd_mesh::domain::core::Scalar
        + FloatElement
        + RealFftData<PlanScalar = T>
        + PlanCacheProvider
        + Copy
        + Neg<Output = T>,
> {
    transform: FourierTransform<T>,
}

impl<T> SpectralDerivative<T>
where
    T: cfd_mesh::domain::core::Scalar
        + FloatElement
        + RealFftData<PlanScalar = T>
        + PlanCacheProvider
        + Copy
        + Neg<Output = T>,
    Complex<T>: Copy + PlanScratch,
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
