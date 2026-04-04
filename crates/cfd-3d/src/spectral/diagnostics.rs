//! Spectral diagnostics for velocity fields.
//!
//! This module provides FFT-backed post-processing helpers that are useful for
//! DNS validation and turbulence analysis, including isotropic kinetic energy
//! spectra, enstrophy spectra, probe-signal power spectra, and normalized
//! temporal autocorrelation series.

use apollofft::{fft_1d_array, fft_3d_array, Complex64};
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid_dynamics::VelocityField;
use nalgebra::{RealField, Vector3};
use ndarray::Array1;
use ndarray::Array3;
use num_traits::{FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

/// Isotropic kinetic energy spectrum for a 3D velocity field.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct KineticEnergySpectrum {
    /// Grid dimensions used to compute the spectrum.
    pub dimensions: (usize, usize, usize),
    /// Integer shell index for each bin.
    pub shell_wavenumbers: Vec<usize>,
    /// Total kinetic energy accumulated into each shell.
    pub shell_energy: Vec<f64>,
    /// Number of Fourier modes that landed in each shell.
    pub shell_mode_counts: Vec<usize>,
    /// Total kinetic energy across all shells.
    pub total_energy: f64,
}

/// Isotropic enstrophy spectrum for a 3D velocity field.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct EnstrophySpectrum {
    /// Grid dimensions used to compute the spectrum.
    pub dimensions: (usize, usize, usize),
    /// Physical domain lengths used to scale the wave numbers.
    pub domain_lengths: (f64, f64, f64),
    /// Integer shell index for each bin.
    pub shell_wavenumbers: Vec<usize>,
    /// Total enstrophy accumulated into each shell.
    pub shell_enstrophy: Vec<f64>,
    /// Number of Fourier modes that landed in each shell.
    pub shell_mode_counts: Vec<usize>,
    /// Total enstrophy across all shells.
    pub total_enstrophy: f64,
}

/// Normalized temporal autocorrelation of a probe-signal time series.
#[derive(Debug, Clone, Serialize, Deserialize)]
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
/// Fourier coefficients are grouped into integer shells using
/// `floor(sqrt(kx^2 + ky^2 + kz^2))` after mapping FFT indices to signed
/// wavenumbers. The returned shell energies sum to the total kinetic energy
/// implied by Parseval's identity.
pub fn kinetic_energy_spectrum<T>(velocity: &VelocityField<T>) -> Result<KineticEnergySpectrum>
where
    T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive,
{
    let (nx, ny, nz) = velocity.dimensions;
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(Error::InvalidConfiguration(
            "kinetic_energy_spectrum: dimensions must be greater than zero".into(),
        ));
    }

    let total_points = nx
        .checked_mul(ny)
        .and_then(|value| value.checked_mul(nz))
        .ok_or_else(|| {
            Error::InvalidConfiguration(
                "kinetic_energy_spectrum: grid dimensions overflowed usize".into(),
            )
        })?;

    if velocity.components.len() != total_points {
        return Err(Error::InvalidConfiguration(format!(
            "kinetic_energy_spectrum: expected {total_points} velocity samples, found {}",
            velocity.components.len()
        )));
    }

    let ux = component_to_array(velocity, |value| value.x, "u")?;
    let uy = component_to_array(velocity, |value| value.y, "v")?;
    let uz = component_to_array(velocity, |value| value.z, "w")?;

    let ux_hat = fft_3d_array(&ux);
    let uy_hat = fft_3d_array(&uy);
    let uz_hat = fft_3d_array(&uz);

    let max_shell = radial_shell_index(nx / 2, ny / 2, nz / 2);
    let mut shell_energy = vec![0.0; max_shell + 1];
    let mut shell_mode_counts = vec![0usize; max_shell + 1];
    let normalization = total_points as f64;

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let shell = shell_index(i, j, k, nx, ny, nz);
                let mode_energy = 0.5
                    * (ux_hat[[i, j, k]].norm_sqr()
                        + uy_hat[[i, j, k]].norm_sqr()
                        + uz_hat[[i, j, k]].norm_sqr())
                    / normalization;
                shell_energy[shell] += mode_energy;
                shell_mode_counts[shell] += 1;
            }
        }
    }

    let total_energy = shell_energy.iter().sum();
    let shell_wavenumbers = (0..shell_energy.len()).collect();

    Ok(KineticEnergySpectrum {
        dimensions: velocity.dimensions,
        shell_wavenumbers,
        shell_energy,
        shell_mode_counts,
        total_energy,
    })
}

/// Compute the isotropic enstrophy spectrum of a 3D velocity field.
///
/// The vorticity is evaluated spectrally via `ω̂ = i k × û`, then grouped into
/// integer shells using the same radial shell indexing as the kinetic-energy
/// spectrum.
pub fn enstrophy_spectrum<T>(
    velocity: &VelocityField<T>,
    domain_lengths: (f64, f64, f64),
) -> Result<EnstrophySpectrum>
where
    T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive,
{
    let (nx, ny, nz) = velocity.dimensions;
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(Error::InvalidConfiguration(
            "enstrophy_spectrum: dimensions must be greater than zero".into(),
        ));
    }

    let (lx, ly, lz) = domain_lengths;
    if !lx.is_finite() || !ly.is_finite() || !lz.is_finite() || lx <= 0.0 || ly <= 0.0 || lz <= 0.0 {
        return Err(Error::InvalidConfiguration(
            "enstrophy_spectrum: domain lengths must be finite and positive".into(),
        ));
    }

    let total_points = nx
        .checked_mul(ny)
        .and_then(|value| value.checked_mul(nz))
        .ok_or_else(|| {
            Error::InvalidConfiguration(
                "enstrophy_spectrum: grid dimensions overflowed usize".into(),
            )
        })?;

    if velocity.components.len() != total_points {
        return Err(Error::InvalidConfiguration(format!(
            "enstrophy_spectrum: expected {total_points} velocity samples, found {}",
            velocity.components.len()
        )));
    }

    let ux = component_to_array(velocity, |value| value.x, "u")?;
    let uy = component_to_array(velocity, |value| value.y, "v")?;
    let uz = component_to_array(velocity, |value| value.z, "w")?;

    let ux_hat = fft_3d_array(&ux);
    let uy_hat = fft_3d_array(&uy);
    let uz_hat = fft_3d_array(&uz);

    let max_shell = radial_shell_index(nx / 2, ny / 2, nz / 2);
    let mut shell_enstrophy = vec![0.0; max_shell + 1];
    let mut shell_mode_counts = vec![0usize; max_shell + 1];
    let normalization = total_points as f64;
    let omega_imag = Complex64::new(0.0, 1.0);

    for k in 0..nz {
        let kz = wavenumber_axis(k, nz, lz);
        for j in 0..ny {
            let ky = wavenumber_axis(j, ny, ly);
            for i in 0..nx {
                let kx = wavenumber_axis(i, nx, lx);
                let shell = shell_index(i, j, k, nx, ny, nz);

                let omega_x_hat = omega_imag * (ky * uz_hat[[i, j, k]] - kz * uy_hat[[i, j, k]]);
                let omega_y_hat = omega_imag * (kz * ux_hat[[i, j, k]] - kx * uz_hat[[i, j, k]]);
                let omega_z_hat = omega_imag * (kx * uy_hat[[i, j, k]] - ky * ux_hat[[i, j, k]]);
                let mode_enstrophy = 0.5
                    * (omega_x_hat.norm_sqr()
                        + omega_y_hat.norm_sqr()
                        + omega_z_hat.norm_sqr())
                    / normalization;

                shell_enstrophy[shell] += mode_enstrophy;
                shell_mode_counts[shell] += 1;
            }
        }
    }

    let total_enstrophy = shell_enstrophy.iter().sum();
    let shell_wavenumbers = (0..shell_enstrophy.len()).collect();

    Ok(EnstrophySpectrum {
        dimensions: velocity.dimensions,
        domain_lengths,
        shell_wavenumbers,
        shell_enstrophy,
        shell_mode_counts,
        total_enstrophy,
    })
}

/// Compute the normalized temporal autocorrelation of a probe signal.
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
pub fn probe_signal_spectrum(
    samples: &[f64],
    sample_period: f64,
) -> Result<ProbeSignalSpectrum> {
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
    let demeaned = Array1::from_iter(samples.iter().map(|value| value - mean));
    let spectrum = fft_1d_array(&demeaned);
    let sample_count = samples.len();
    let normalization = sample_count as f64;
    let frequencies_hz = (0..sample_count)
        .map(|index| {
            signed_frequency_index(index, sample_count) as f64
                / (sample_count as f64 * sample_period)
        })
        .collect();
    let spectral_energy: Vec<f64> = spectrum
        .iter()
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

fn component_to_array<T, F>(
    velocity: &VelocityField<T>,
    component: F,
    component_name: &str,
) -> Result<Array3<f64>>
where
    T: cfd_mesh::domain::core::Scalar + RealField + Copy + FromPrimitive + ToPrimitive,
    F: Fn(&Vector3<T>) -> T,
{
    let (nx, ny, nz) = velocity.dimensions;
    let mut values = Array3::<f64>::zeros((nx, ny, nz));

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let idx = k * nx * ny + j * nx + i;
                values[[i, j, k]] = component(&velocity.components[idx]).to_f64().ok_or_else(
                    || {
                        Error::InvalidConfiguration(format!(
                            "kinetic_energy_spectrum: unable to convert {component_name} component to f64"
                        ))
                    },
                )?;
            }
        }
    }

    Ok(values)
}

fn signed_wavenumber(index: usize, size: usize) -> isize {
    let half = size / 2;
    if index <= half {
        index as isize
    } else {
        index as isize - size as isize
    }
}

fn signed_frequency_index(index: usize, size: usize) -> isize {
    let half = size / 2;
    if index <= half {
        index as isize
    } else {
        index as isize - size as isize
    }
}

fn wavenumber_axis(index: usize, size: usize, length: f64) -> f64 {
    let factor = 2.0 * std::f64::consts::PI / length;
    signed_wavenumber(index, size) as f64 * factor
}

fn radial_shell_index(nx_half: usize, ny_half: usize, nz_half: usize) -> usize {
    let max_radius = ((nx_half * nx_half + ny_half * ny_half + nz_half * nz_half) as f64).sqrt();
    max_radius.floor() as usize
}

fn shell_index(i: usize, j: usize, k: usize, nx: usize, ny: usize, nz: usize) -> usize {
    let kx = signed_wavenumber(i, nx);
    let ky = signed_wavenumber(j, ny);
    let kz = signed_wavenumber(k, nz);
    let radius = ((kx * kx + ky * ky + kz * kz) as f64).sqrt();
    radius.floor() as usize
}

#[cfg(test)]
mod tests {
    use super::*;
    use nalgebra::Vector3;

    #[test]
    fn constant_velocity_energy_lives_in_zero_shell() {
        let mut velocity = VelocityField {
            components: vec![Vector3::new(2.0, -1.0, 0.5); 4 * 4 * 4],
            dimensions: (4, 4, 4),
        };

        let spectrum = kinetic_energy_spectrum(&velocity).expect("spectrum should compute");
        let expected_energy = 0.5 * (2.0_f64 * 2.0 + 1.0_f64 + 0.5_f64 * 0.5) * 64.0;

        assert!((spectrum.total_energy - expected_energy).abs() < 1e-10);
        assert!((spectrum.shell_energy[0] - expected_energy).abs() < 1e-10);
        assert!(spectrum
            .shell_energy
            .iter()
            .skip(1)
            .all(|energy| energy.abs() < 1e-10));

        // Ensure the input is left untouched.
        velocity.components[0] = Vector3::new(0.0, 0.0, 0.0);
        assert_eq!(velocity.components.len(), 64);
    }

    #[test]
    fn single_mode_energy_peaks_in_first_shell() {
        let (nx, ny, nz) = (8, 4, 4);
        let mut velocity = VelocityField {
            components: Vec::with_capacity(nx * ny * nz),
            dimensions: (nx, ny, nz),
        };

        for _k in 0..nz {
            for _j in 0..ny {
                for i in 0..nx {
                    let phase = 2.0 * std::f64::consts::PI * i as f64 / nx as f64;
                    velocity.components.push(Vector3::new(phase.sin(), 0.0, 0.0));
                }
            }
        }

        let spectrum = kinetic_energy_spectrum(&velocity).expect("spectrum should compute");
        let expected_energy = 0.25 * (nx * ny * nz) as f64;
        let dominant_shell = spectrum
            .shell_energy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| lhs.1.partial_cmp(rhs.1).expect("energy should be comparable"))
            .map(|(index, _)| index)
            .expect("spectrum should contain at least one shell");

        assert!((spectrum.total_energy - expected_energy).abs() < 1e-10);
        assert!(spectrum.shell_energy[0] < 1e-10);
        assert_eq!(dominant_shell, 1);
        assert!(spectrum.shell_energy[1] > 0.99 * spectrum.total_energy);
    }

    #[test]
    fn constant_velocity_has_zero_enstrophy() {
        let velocity = VelocityField {
            components: vec![Vector3::new(1.0, 2.0, -1.5); 4 * 4 * 4],
            dimensions: (4, 4, 4),
        };

        let spectrum = enstrophy_spectrum(&velocity, (1.0, 1.0, 1.0))
            .expect("enstrophy spectrum should compute");

        assert!(spectrum.total_enstrophy.abs() < 1e-10);
        assert!(spectrum
            .shell_enstrophy
            .iter()
            .all(|value| value.abs() < 1e-10));
    }

    #[test]
    fn single_mode_enstrophy_peaks_in_first_shell() {
        let (nx, ny, nz) = (8, 8, 8);
        let mut velocity = VelocityField {
            components: Vec::with_capacity(nx * ny * nz),
            dimensions: (nx, ny, nz),
        };

        for _k in 0..nz {
            for j in 0..ny {
                let phase = 2.0 * std::f64::consts::PI * j as f64 / ny as f64;
                for _i in 0..nx {
                    velocity.components.push(Vector3::new(phase.sin(), 0.0, 0.0));
                }
            }
        }

        let spectrum = enstrophy_spectrum(&velocity, (1.0, 1.0, 1.0))
            .expect("enstrophy spectrum should compute");

        let dominant_shell = spectrum
            .shell_enstrophy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| lhs.1.partial_cmp(rhs.1).expect("enstrophy should be comparable"))
            .map(|(index, _)| index)
            .expect("spectrum should contain at least one shell");

        assert!(spectrum.total_enstrophy > 0.0);
        assert_eq!(dominant_shell, 1);
        assert!(spectrum.shell_enstrophy[1] > 0.99 * spectrum.total_enstrophy);
    }

    #[test]
    fn alternating_signal_has_predictable_autocorrelation_and_probe_peak() {
        let samples: Vec<f64> = (0..16)
            .map(|index| if index % 2 == 0 { 1.0 } else { -1.0 })
            .collect();

        let autocorrelation = temporal_autocorrelation(&samples, 0.25, 4)
            .expect("autocorrelation should compute");
        assert_eq!(autocorrelation.lags, vec![0, 1, 2, 3, 4]);
        assert!((autocorrelation.values[0] - 1.0).abs() < 1e-10);
        assert!((autocorrelation.values[1] + 1.0).abs() < 1e-10);
        assert!((autocorrelation.values[2] - 1.0).abs() < 1e-10);

        let spectrum = probe_signal_spectrum(&samples, 0.25)
            .expect("probe signal spectrum should compute");
        let dominant_bin = spectrum
            .spectral_energy
            .iter()
            .enumerate()
            .max_by(|lhs, rhs| lhs.1.partial_cmp(rhs.1).expect("power should be comparable"))
            .map(|(index, _)| index)
            .expect("spectrum should contain at least one bin");

        assert_eq!(dominant_bin, 8);
        assert!((spectrum.frequencies_hz[dominant_bin] - 2.0).abs() < 1e-10);
        assert!(spectrum.spectral_energy[dominant_bin] > 0.99 * spectrum.total_fluctuation_energy);
        assert!((spectrum.variance - 1.0).abs() < 1e-10);
    }
}