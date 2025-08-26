//! Fourier transform operations for spectral methods
//!
//! Reference: Canuto et al. (2006). "Spectral Methods: Fundamentals in Single Domains"

use cfd_core::error::Result;
use cfd_core::numeric;
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
    /// Forward discrete Fourier transform (DFT)
    /// Uses naive O(n²) algorithm - production code should use FFT
    pub fn forward(&self, u: &DVector<T>) -> Result<DVector<Complex<T>>> {
        let n = self.n;
        let mut u_hat = DVector::zeros(n);
        let two_pi = T::from_f64(2.0 * PI).ok_or_else(|| {
            cfd_core::error::Error::InvalidConfiguration("Cannot convert constant".into())
        })?;
        for k in 0..n {
            let mut sum = Complex::new(T::zero(), T::zero());
            for j in 0..n {
                let phase =
                    -two_pi * self.wavenumbers[k] * cfd_core::numeric::from_usize(j)?
                        / cfd_core::numeric::from_usize(n)?;
                let exp = Complex::new(phase.cos(), phase.sin());
                sum += exp * u[j];
            }
            u_hat[k] = sum / cfd_core::numeric::from_usize(n)?;
        Ok(u_hat)
    /// Inverse discrete Fourier transform (IDFT)
    pub fn inverse(&self, u_hat: &DVector<Complex<T>>) -> Result<DVector<T>> {
        let mut u = DVector::zeros(n);
        for j in 0..n {
            for k in 0..n {
                    two_pi * self.wavenumbers[k] * cfd_core::numeric::from_usize(j)?
                sum += exp * u_hat[k];
            u[j] = sum.re;
        Ok(u)
    /// Get wavenumbers
    #[must_use]
    pub fn wavenumbers(&self) -> &[T] {
        &self.wavenumbers
/// Spectral derivative computation
pub struct SpectralDerivative<T: RealField + Copy> {
    transform: FourierTransform<T>,
impl<T: RealField + FromPrimitive + Copy> SpectralDerivative<T> {
        Ok(Self {
            transform: FourierTransform::new(n)?,
        })
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
        // Transform back to physical space
        self.transform.inverse(&du_hat)
