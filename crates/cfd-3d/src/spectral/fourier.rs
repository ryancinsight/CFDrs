//! Fourier transform operations for spectral methods
//!
//! Implements efficient FFT using Cooley-Tukey algorithm (O(n log n))
//! Reference: Canuto et al. (2006). "Spectral Methods: Fundamentals in Single Domains"
//! Cooley-Tukey: Cooley & Tukey (1965). "An algorithm for the machine calculation of complex Fourier series"

use cfd_core::error::Result;
use nalgebra::{Complex, DVector, RealField};
use num_traits::FromPrimitive;
use std::f64::consts::PI;

/// Fourier transform operations
pub struct FourierTransform<T: RealField + Copy> {
    /// Number of modes
    n: usize,
    /// Wavenumbers
    wavenumbers: Vec<T>,
}

impl<T: RealField + FromPrimitive + Copy> FourierTransform<T> {
    /// Create new Fourier transform operator
    pub fn new(n: usize) -> Result<Self> {
        let mut wavenumbers = Vec::with_capacity(n);

        // Compute wavenumbers k = 0, 1, ..., n/2, -n/2+1, ..., -1
        for i in 0..n {
            let k = if i <= n / 2 {
                i as f64
            } else {
                (i as f64) - (n as f64)
            };

            let wn = T::from_f64(k).ok_or_else(|| {
                cfd_core::error::Error::InvalidConfiguration("Cannot convert wavenumber".into())
            })?;
            wavenumbers.push(wn);
        }

        Ok(Self { n, wavenumbers })
    }

    /// Forward Fast Fourier Transform (FFT) using Cooley-Tukey algorithm
    ///
    /// Implements efficient O(n log n) FFT with fallback to DFT for non-power-of-2 sizes
    /// Requires n to be a power of 2 for optimal performance using radix-2 algorithm
    ///
    /// # Algorithm
    /// 1. Check if n is power of 2
    /// 2. If yes: Use in-place Cooley-Tukey FFT (bit-reversal + butterfly operations)
    /// 3. If no: Fall back to O(n²) DFT with proper scaling
    ///
    /// # References
    /// - Cooley, J.W. and Tukey, J.W. (1965). "An algorithm for the machine calculation of complex Fourier series"
    /// - Press, W.H. et al. (1992). "Numerical Recipes in C" §12.2
    pub fn forward(&self, u: &DVector<T>) -> Result<DVector<Complex<T>>> {
        let n = self.n;
        let scale = T::from_usize(n).unwrap_or(T::one());

        // Check if n is power of 2 (required for efficient radix-2 FFT)
        if n.is_power_of_two() {
            // Use efficient FFT for power-of-2 sizes
            let mut data: Vec<Complex<T>> = u.iter().map(|&x| Complex::new(x, T::zero())).collect();
            self.fft_inplace(&mut data, false);

            // Normalize by 1/n to get physical amplitudes
            for val in &mut data {
                *val /= scale;
            }

            Ok(DVector::from_vec(data))
        } else {
            // Fall back to DFT for non-power-of-2 sizes
            let mut u_hat = self.dft_forward(u)?;
            for val in u_hat.iter_mut() {
                *val /= scale;
            }
            Ok(u_hat)
        }
    }

    /// In-place Cooley-Tukey FFT implementation
    fn fft_inplace(&self, data: &mut [Complex<T>], inverse: bool) {
        let n = data.len();
        let sign = if inverse { T::one() } else { -T::one() };

        // Bit-reversal permutation
        let mut j = 0usize;
        for i in 1..n {
            let mut bit = n >> 1;
            while j & bit != 0 {
                j ^= bit;
                bit >>= 1;
            }
            j ^= bit;

            if i < j {
                data.swap(i, j);
            }
        }

        // Iterative FFT using Danielson-Lanczos lemma
        let mut length = 2;
        while length <= n {
            let angle = sign * T::from_f64(2.0 * PI).unwrap_or(T::one())
                / T::from_usize(length).unwrap_or(T::one());
            let wlen = Complex::new(angle.cos(), angle.sin());

            for i in (0..n).step_by(length) {
                let mut w = Complex::new(T::one(), T::zero());
                for j in 0..length / 2 {
                    let u_val = data[i + j];
                    let v_val = data[i + j + length / 2] * w;
                    data[i + j] = u_val + v_val;
                    data[i + j + length / 2] = u_val - v_val;
                    w *= wlen;
                }
            }
            length <<= 1;
        }
    }

    /// Fallback DFT for non-power-of-2 sizes (O(n²) complexity)
    fn dft_forward(&self, u: &DVector<T>) -> Result<DVector<Complex<T>>> {
        let n = self.n;
        let mut u_hat = DVector::zeros(n);
        let two_pi = T::from_f64(2.0 * PI).ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration("Cannot convert constant".into())
        })?;

        for k in 0..n {
            let mut sum = Complex::new(T::zero(), T::zero());
            for j in 0..n {
                let phase =
                    -two_pi * self.wavenumbers[k] * T::from_usize(j).unwrap_or_else(|| T::zero())
                        / T::from_usize(n).unwrap_or_else(|| T::zero());
                let exp = Complex::new(phase.cos(), phase.sin());
                sum += exp * Complex::new(u[j], T::zero());
            }
            u_hat[k] = sum;
        }

        Ok(u_hat)
    }

    /// Inverse Fast Fourier Transform (IFFT)
    ///
    /// Implements efficient inverse FFT with proper scaling (1/n)
    /// Uses FFT for power-of-2 sizes, DFT fallback for others
    pub fn inverse(&self, u_hat: &DVector<Complex<T>>) -> Result<DVector<T>> {
        let n = self.n;

        if n.is_power_of_two() {
            // Use efficient IFFT for power-of-2 sizes
            let mut data = u_hat.iter().copied().collect::<Vec<_>>();
            self.fft_inplace(&mut data, true);

            // Extract real part (normalization already done in forward)
            let real_part: Vec<T> = data.into_iter().map(|x| x.re).collect();

            Ok(DVector::from_vec(real_part))
        } else {
            // Fall back to IDFT for non-power-of-2 sizes
            self.dft_inverse_real(u_hat)
        }
    }

    /// Fallback inverse DFT returning real values
    fn dft_inverse_real(&self, u_hat: &DVector<Complex<T>>) -> Result<DVector<T>> {
        let n = self.n;
        let mut u = DVector::zeros(n);
        let two_pi = T::from_f64(2.0 * PI).ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration("Cannot convert constant".into())
        })?;

        for j in 0..n {
            let mut sum = Complex::new(T::zero(), T::zero());
            for k in 0..n {
                let phase =
                    two_pi * self.wavenumbers[k] * T::from_usize(j).unwrap_or_else(|| T::zero())
                        / T::from_usize(n).unwrap_or_else(|| T::zero());
                let exp = Complex::new(phase.cos(), phase.sin());
                sum += exp * u_hat[k];
            }
            u[j] = sum.re;
        }

        Ok(u)
    }

    /// Get wavenumbers
    #[must_use]
    pub fn wavenumbers(&self) -> &[T] {
        &self.wavenumbers
    }
}

/// Spectral derivative computation
pub struct SpectralDerivative<T: RealField + Copy> {
    transform: FourierTransform<T>,
}

impl<T: RealField + FromPrimitive + Copy> SpectralDerivative<T> {
    /// Create a new spectral derivative operator for the given grid size
    ///
    /// # Arguments
    /// * `n` - Number of grid points (must be power of 2 for optimal FFT performance)
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
