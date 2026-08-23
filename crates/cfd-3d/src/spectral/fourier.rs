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
        let zero = Complex::new(<T as NumericElement>::ZERO, <T as NumericElement>::ZERO);

        for (k_idx, &k) in self.transform.wavenumbers().iter().enumerate() {
            // Even-length transforms: the Nyquist mode has no representable
            // derivative in the N-point trigonometric interpolant space, so
            // its factor is zeroed (Trefethen, Spectral Methods in MATLAB,
            // ch. 2). Multiplying the real Nyquist coefficient by ik would
            // inject a component outside the approximation space and break
            // conjugate symmetry.
            let is_nyquist = self.transform.n % 2 == 0 && k_idx == self.transform.n / 2;
            let factor = if order > 0 && is_nyquist {
                zero
            } else {
                i * k
            };
            let factor = factor.powi(order as i32);
            du_hat[[k_idx]] = factor * u_hat[[k_idx]];
        }

        // Transform back to physical space
        self.transform.inverse(&du_hat)
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use leto::Array1;

    fn array_from_fn(n: usize, mut f: impl FnMut(usize) -> f64) -> Array1<f64> {
        Array1::from_vec([n], (0..n).map(&mut f).collect())
            .expect("signal length matches declared shape")
    }

    /// The even-length Nyquist mode has no representable derivative in the
    /// N-point trigonometric interpolant space; its factor must be zeroed
    /// (Trefethen, Spectral Methods in MATLAB, ch. 2).
    #[test]
    fn nyquist_mode_derivative_is_zeroed() {
        let n = 4_usize;
        let derivative = SpectralDerivative::<f64>::new(n).expect("valid size");

        // Pure Nyquist mode: cos(2x) sampled at x_j = πj/2 → [1, −1, 1, −1].
        let signal = array_from_fn(n, |j| {
            if j % 2 == 0 {
                1.0
            } else {
                -1.0
            }
        });

        let d = derivative
            .derivative(&signal, 1)
            .expect("first derivative must succeed");

        // Exact derivative −2·sin(2x) vanishes at every sample point.
        for (idx, &value) in d.iter().enumerate() {
            assert!(
                value.abs() < 1e-12,
                "du/dx[{idx}] = {} must vanish for the Nyquist mode",
                value
            );
        }
    }

    /// A non-Nyquist mode must still differentiate correctly.
    #[test]
    fn first_derivative_of_sine_matches_cosine() {
        let n = 8_usize;
        let derivative = SpectralDerivative::<f64>::new(n).expect("valid size");

        let signal = array_from_fn(n, |j| {
            let x = 2.0 * std::f64::consts::PI * j as f64 / n as f64;
            x.sin()
        });

        let d = derivative
            .derivative(&signal, 1)
            .expect("first derivative must succeed");

        for (idx, &value) in d.iter().enumerate() {
            let expected = {
                let x = 2.0 * std::f64::consts::PI * idx as f64 / n as f64;
                x.cos()
            };
            assert!(
                (value - expected).abs() < 1e-12,
                "du/dx[{idx}] = {value}, expected {expected}"
            );
        }
    }
}
