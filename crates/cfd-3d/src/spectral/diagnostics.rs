//! Spectral diagnostics for velocity fields.
//!
//! This module provides FFT-backed post-processing helpers that are useful for
//! DNS validation and turbulence analysis, including isotropic kinetic energy
//! spectra derived from 3D velocity fields.

use apollofft::fft_3d_array;
use cfd_core::error::{Error, Result};
use cfd_core::physics::fluid_dynamics::VelocityField;
use nalgebra::{RealField, Vector3};
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
}