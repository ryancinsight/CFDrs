//! Spatial spectra: isotropic kinetic energy and enstrophy.
//!
//! Split out of `diagnostics.rs` so the FFT-backed 3D shell integration keeps
//! its own file, separate from the 1D probe-signal diagnostics.

use crate::atlas_array::{fft_3d_array, set, value, values};
use apollo_fft::Complex64;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid_dynamics::VelocityField;
use eunomia::NumericElement;
use leto::Array3;
use leto::geometry::Vector3;
use serde::{Deserialize, Serialize};

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

pub fn kinetic_energy_spectrum<T: NumericElement>(
    velocity: &VelocityField<T>,
) -> Result<KineticEnergySpectrum> {
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

    let ux = component_samples_as_f64(velocity, |value| value.x)?;
    let uy = component_samples_as_f64(velocity, |value| value.y)?;
    let uz = component_samples_as_f64(velocity, |value| value.z)?;

    let ux_hat = fft_3d_array(&ux)?;
    let uy_hat = fft_3d_array(&uy)?;
    let uz_hat = fft_3d_array(&uz)?;

    let max_shell = radial_shell_index(nx / 2, ny / 2, nz / 2);
    let mut shell_energy = vec![0.0; max_shell + 1];
    let mut shell_mode_counts = vec![0usize; max_shell + 1];
    let normalization = total_points as f64;

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let shell = shell_index(i, j, k, nx, ny, nz);
                let index = [i, j, k];
                let mode_energy = 0.5
                    * (value(&ux_hat, index).norm_sqr()
                        + value(&uy_hat, index).norm_sqr()
                        + value(&uz_hat, index).norm_sqr())
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
/// # Precision boundary
///
/// Same as [`kinetic_energy_spectrum`]: the samples are lowered to `f64` once
/// at the top of this routine and the whole transform is `f64`.
///
/// The vorticity is evaluated spectrally via `ω̂ = i k × û`, then grouped into
/// integer shells using the same radial shell indexing as the kinetic-energy
/// spectrum.
pub fn enstrophy_spectrum<T: NumericElement>(
    velocity: &VelocityField<T>,
    domain_lengths: (f64, f64, f64),
) -> Result<EnstrophySpectrum> {
    let (nx, ny, nz) = velocity.dimensions;
    if nx == 0 || ny == 0 || nz == 0 {
        return Err(Error::InvalidConfiguration(
            "enstrophy_spectrum: dimensions must be greater than zero".into(),
        ));
    }

    let (lx, ly, lz) = domain_lengths;
    if !lx.is_finite() || !ly.is_finite() || !lz.is_finite() || lx <= 0.0 || ly <= 0.0 || lz <= 0.0
    {
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

    let ux = component_samples_as_f64(velocity, |value| value.x)?;
    let uy = component_samples_as_f64(velocity, |value| value.y)?;
    let uz = component_samples_as_f64(velocity, |value| value.z)?;

    let ux_hat = fft_3d_array(&ux)?;
    let uy_hat = fft_3d_array(&uy)?;
    let uz_hat = fft_3d_array(&uz)?;

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

                let index = [i, j, k];
                let omega_x_hat =
                    omega_imag * (value(&uz_hat, index) * ky - value(&uy_hat, index) * kz);
                let omega_y_hat =
                    omega_imag * (value(&ux_hat, index) * kz - value(&uz_hat, index) * kx);
                let omega_z_hat =
                    omega_imag * (value(&uy_hat, index) * kx - value(&ux_hat, index) * ky);
                let mode_enstrophy = 0.5
                    * (omega_x_hat.norm_sqr() + omega_y_hat.norm_sqr() + omega_z_hat.norm_sqr())
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

fn component_samples_as_f64<T, F>(velocity: &VelocityField<T>, component: F) -> Result<Array3<f64>>
where
    T: NumericElement,
    F: Fn(&Vector3<T>) -> T,
{
    let (nx, ny, nz) = velocity.dimensions;
    let mut values = Array3::<f64>::zeros([nx, ny, nz]);

    for k in 0..nz {
        for j in 0..ny {
            for i in 0..nx {
                let idx = k * nx * ny + j * nx + i;
                set(
                    &mut values,
                    [i, j, k],
                    <T as NumericElement>::to_f64(component(&velocity.components[idx])),
                );
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

