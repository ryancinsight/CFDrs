//! Temporal diagnostics: probe-signal spectra and autocorrelation.
//!
//! Split out of `diagnostics.rs` so the 1D time-series analysis keeps its own
//! file, separate from the 3D shell-integrated spectra.

use crate::atlas_array::{fft_1d_array, value, values};
use cfd_core::error::{Error, Result};
use leto::Array1;
use serde::{Deserialize, Serialize};

pub struct TemporalAutocorrelation {
    /// Number of samples in the original time series.
    pub sample_count: usize,
    /// Sample period used to convert lags into time.
    pub sample_period: f64,
    /// Integer lags included in the report.
    pub lags: Vec<usize>,
    /// Lag times corresponding to `lags`.
    pub lag_times: Vec<f64>,
    /// Normalized autocorrelation values.
    pub values: Vec<f64>,
    /// Mean of the original signal.
    pub mean: f64,
    /// Variance of the demeaned signal.
    pub variance: f64,
}

/// One-dimensional probe-signal power spectrum.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ProbeSignalSpectrum {
    /// Number of samples in the original time series.
    pub sample_count: usize,
    /// Sample period used to convert DFT bins into physical frequencies.
    pub sample_period: f64,
    /// Signed frequency of each bin in hertz.
    pub frequencies_hz: Vec<f64>,
    /// Per-bin fluctuation energy / power.
    pub spectral_energy: Vec<f64>,
    /// Mean of the original signal.
    pub mean: f64,
    /// Variance of the demeaned signal.
    pub variance: f64,
    /// Total fluctuation energy represented by the spectrum.
    pub total_fluctuation_energy: f64,
}

/// Compute the isotropic kinetic energy spectrum of a 3D velocity field.
///
/// # Precision boundary
///
/// The spectrum is a diagnostic, not a simulation state: the FFT backend
/// (`apollo_fft`) and every field of [`KineticEnergySpectrum`] are `f64`. The
/// field samples are therefore lowered to `f64` once, at the top of this
/// routine, and the transform runs entirely in double
/// precision regardless of `T`. Nothing is narrowed back into `T`.
///
/// Fourier coefficients are grouped into integer shells using
/// `floor(sqrt(kx^2 + ky^2 + kz^2))` after mapping FFT indices to signed
/// wavenumbers. The returned shell energies sum to the total kinetic energy
/// implied by Parseval's identity.

pub fn temporal_autocorrelation(
    samples: &[f64],
    sample_period: f64,
    max_lag: usize,
) -> Result<TemporalAutocorrelation> {
    if samples.is_empty() {
        return Err(Error::InvalidConfiguration(
            "temporal_autocorrelation: signal must contain at least one sample".into(),
        ));
    }
    if !sample_period.is_finite() || sample_period <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "temporal_autocorrelation: sample period must be finite and positive".into(),
        ));
    }

    let mean = samples.iter().sum::<f64>() / samples.len() as f64;
    let demeaned: Vec<f64> = samples.iter().map(|value| value - mean).collect();
    let variance = demeaned.iter().map(|value| value * value).sum::<f64>() / samples.len() as f64;
    if variance <= f64::EPSILON {
        return Err(Error::InvalidConfiguration(
            "temporal_autocorrelation: signal variance must be positive".into(),
        ));
    }

    let max_lag = max_lag.min(samples.len().saturating_sub(1));
    let mut lags = Vec::with_capacity(max_lag + 1);
    let mut lag_times = Vec::with_capacity(max_lag + 1);
    let mut values = Vec::with_capacity(max_lag + 1);

    for lag in 0..=max_lag {
        let mut correlation = 0.0;
        let sample_count = samples.len() - lag;
        for index in 0..sample_count {
            correlation += demeaned[index] * demeaned[index + lag];
        }
        correlation /= sample_count as f64;
        lags.push(lag);
        lag_times.push(lag as f64 * sample_period);
        values.push(correlation / variance);
    }

    Ok(TemporalAutocorrelation {
        sample_count: samples.len(),
        sample_period,
        lags,
        lag_times,
        values,
        mean,
        variance,
    })
}

/// Compute the FFT-backed power spectrum of a probe signal.
pub fn probe_signal_spectrum(samples: &[f64], sample_period: f64) -> Result<ProbeSignalSpectrum> {
    if samples.is_empty() {
        return Err(Error::InvalidConfiguration(
            "probe_signal_spectrum: signal must contain at least one sample".into(),
        ));
    }
    if !sample_period.is_finite() || sample_period <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "probe_signal_spectrum: sample period must be finite and positive".into(),
        ));
    }

    let mean = samples.iter().sum::<f64>() / samples.len() as f64;
    let demeaned = Array1::from_shape_fn([samples.len()], |[index]| samples[index] - mean);
    let spectrum = fft_1d_array(&demeaned)?;
    let sample_count = samples.len();
    let normalization = sample_count as f64;
    let frequencies_hz = (0..sample_count)
        .map(|index| {
            signed_frequency_index(index, sample_count) as f64
                / (sample_count as f64 * sample_period)
        })
        .collect();
    let spectral_energy: Vec<f64> = values(&spectrum)
        .into_iter()
        .map(|value| value.norm_sqr() / normalization)
        .collect();
    let total_fluctuation_energy = spectral_energy.iter().sum();
    let variance = total_fluctuation_energy / sample_count as f64;

    Ok(ProbeSignalSpectrum {
        sample_count,
        sample_period,
        frequencies_hz,
        spectral_energy,
        mean,
        variance,
        total_fluctuation_energy,
    })
}

/// The single precision boundary of this module: read one velocity component
/// out of the field and lower it to `f64` for the double-precision FFT.
///
/// `T` is consumed here and never reconstructed, so nothing downstream is a
/// generic computation pretending to run in `T`.

fn signed_frequency_index(index: usize, size: usize) -> isize {
    let half = size / 2;
    if index <= half {
        index as isize
    } else {
        index as isize - size as isize
    }
}

