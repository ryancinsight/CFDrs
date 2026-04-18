//! NUFFT-backed coupling helpers for immersed-boundary marker and probe transfer.

use apollofft::Complex64;
use apollofft::{fft_3d_array, ifft_3d_array, nufft_type1_3d, nufft_type2_3d, UniformGrid3D};
use cfd_core::{
    error::{Error, Result},
    physics::fluid_dynamics::VelocityField,
};
use nalgebra::Vector3;
use ndarray::Array3;

use super::lagrangian::LagrangianPoint;

/// NUFFT-based coupling for irregular probes and Lagrangian markers.
#[derive(Debug, Clone)]
pub struct NufftMarkerCoupler3D {
    grid: UniformGrid3D,
}

impl NufftMarkerCoupler3D {
    /// Create a new coupler from a validated periodic grid.
    #[must_use]
    pub fn new(grid: UniformGrid3D) -> Self {
        Self { grid }
    }

    /// Build a coupler from grid dimensions and physical lengths.
    pub fn from_dimensions(
        dimensions: (usize, usize, usize),
        lengths: (f64, f64, f64),
    ) -> Result<Self> {
        let (nx, ny, nz) = dimensions;
        let (lx, ly, lz) = lengths;
        if nx == 0 || ny == 0 || nz == 0 {
            return Err(Error::InvalidConfiguration(
                "grid dimensions must be strictly positive".to_string(),
            ));
        }
        if !lx.is_finite()
            || !ly.is_finite()
            || !lz.is_finite()
            || lx <= 0.0
            || ly <= 0.0
            || lz <= 0.0
        {
            return Err(Error::InvalidConfiguration(
                "domain lengths must be finite and strictly positive".to_string(),
            ));
        }

        let grid = UniformGrid3D::new(nx, ny, nz, lx / nx as f64, ly / ny as f64, lz / nz as f64)
            .map_err(|error| Error::InvalidConfiguration(error.to_string()))?;
        Ok(Self::new(grid))
    }

    /// Return the periodic grid used by the coupling operators.
    #[must_use]
    pub fn grid(&self) -> UniformGrid3D {
        self.grid
    }

    /// Sample a scalar field at irregular probe locations.
    pub fn sample_scalar_field(
        &self,
        field: &Array3<f64>,
        positions: &[Vector3<f64>],
    ) -> Result<Vec<f64>> {
        self.validate_scalar_field(field)?;

        let field_hat = fft_3d_array(field);
        let normalized_hat = self.normalize_spectrum(&field_hat);
        let probe_positions = self.positions_to_tuples(positions);
        let samples = nufft_type2_3d(&normalized_hat, &probe_positions, self.grid);

        Ok(samples.into_iter().map(|sample| sample.re).collect())
    }

    /// Sample a velocity field at irregular probe locations.
    pub fn sample_velocity_at_positions(
        &self,
        field: &VelocityField<f64>,
        positions: &[Vector3<f64>],
    ) -> Result<Vec<Vector3<f64>>> {
        self.validate_velocity_field(field)?;

        let ux = self.component_to_array(field, 0)?;
        let uy = self.component_to_array(field, 1)?;
        let uz = self.component_to_array(field, 2)?;

        let ux_hat = self.normalize_spectrum(&fft_3d_array(&ux));
        let uy_hat = self.normalize_spectrum(&fft_3d_array(&uy));
        let uz_hat = self.normalize_spectrum(&fft_3d_array(&uz));

        let probe_positions = self.positions_to_tuples(positions);
        let sampled_x = nufft_type2_3d(&ux_hat, &probe_positions, self.grid);
        let sampled_y = nufft_type2_3d(&uy_hat, &probe_positions, self.grid);
        let sampled_z = nufft_type2_3d(&uz_hat, &probe_positions, self.grid);

        Ok(sampled_x
            .into_iter()
            .zip(sampled_y)
            .zip(sampled_z)
            .map(|((x, y), z)| Vector3::new(x.re, y.re, z.re))
            .collect())
    }

    /// Sample a velocity field at Lagrangian marker locations.
    pub fn sample_velocity_at_markers(
        &self,
        field: &VelocityField<f64>,
        markers: &[LagrangianPoint<f64>],
    ) -> Result<Vec<Vector3<f64>>> {
        let positions: Vec<_> = markers.iter().map(|marker| marker.position).collect();
        self.sample_velocity_at_positions(field, &positions)
    }

    /// Spread Lagrangian marker forces onto the Eulerian grid.
    pub fn spread_marker_forces(
        &self,
        markers: &[LagrangianPoint<f64>],
    ) -> Result<VelocityField<f64>> {
        let positions = self.marker_positions(markers);
        let forces_x: Vec<_> = markers
            .iter()
            .map(|marker| Complex64::new(marker.force.x * marker.weight, 0.0))
            .collect();
        let forces_y: Vec<_> = markers
            .iter()
            .map(|marker| Complex64::new(marker.force.y * marker.weight, 0.0))
            .collect();
        let forces_z: Vec<_> = markers
            .iter()
            .map(|marker| Complex64::new(marker.force.z * marker.weight, 0.0))
            .collect();

        let force_hat_x = nufft_type1_3d(&positions, &forces_x, self.grid);
        let force_hat_y = nufft_type1_3d(&positions, &forces_y, self.grid);
        let force_hat_z = nufft_type1_3d(&positions, &forces_z, self.grid);

        let force_x = ifft_3d_array(&force_hat_x);
        let force_y = ifft_3d_array(&force_hat_y);
        let force_z = ifft_3d_array(&force_hat_z);

        Ok(self.velocity_field_from_components(force_x, force_y, force_z))
    }

    fn validate_scalar_field(&self, field: &Array3<f64>) -> Result<()> {
        if field.dim() != (self.grid.nx, self.grid.ny, self.grid.nz) {
            return Err(Error::DimensionMismatch {
                expected: self.grid.nx * self.grid.ny * self.grid.nz,
                actual: field.len(),
            });
        }
        Ok(())
    }

    fn validate_velocity_field(&self, field: &VelocityField<f64>) -> Result<()> {
        let expected = self.grid.nx * self.grid.ny * self.grid.nz;
        if field.dimensions != (self.grid.nx, self.grid.ny, self.grid.nz) {
            return Err(Error::DimensionMismatch {
                expected,
                actual: field.components.len(),
            });
        }
        if field.components.len() != expected {
            return Err(Error::DimensionMismatch {
                expected,
                actual: field.components.len(),
            });
        }
        Ok(())
    }

    fn component_to_array(
        &self,
        field: &VelocityField<f64>,
        component: usize,
    ) -> Result<Array3<f64>> {
        if component > 2 {
            return Err(Error::InvalidInput(format!(
                "invalid velocity component index {component}"
            )));
        }

        Ok(Array3::from_shape_fn(
            (self.grid.nx, self.grid.ny, self.grid.nz),
            |(i, j, k)| {
                let index = k * self.grid.nx * self.grid.ny + j * self.grid.nx + i;
                field.components[index][component]
            },
        ))
    }

    fn velocity_field_from_components(
        &self,
        x: Array3<f64>,
        y: Array3<f64>,
        z: Array3<f64>,
    ) -> VelocityField<f64> {
        let mut components = Vec::with_capacity(self.grid.nx * self.grid.ny * self.grid.nz);

        for k in 0..self.grid.nz {
            for j in 0..self.grid.ny {
                for i in 0..self.grid.nx {
                    components.push(Vector3::new(x[[i, j, k]], y[[i, j, k]], z[[i, j, k]]));
                }
            }
        }

        VelocityField {
            components,
            dimensions: (self.grid.nx, self.grid.ny, self.grid.nz),
        }
    }

    fn positions_to_tuples(&self, positions: &[Vector3<f64>]) -> Vec<(f64, f64, f64)> {
        let (lx, ly, lz) = self.grid.lengths();
        positions
            .iter()
            .map(|position| {
                (
                    position.x.rem_euclid(lx),
                    position.y.rem_euclid(ly),
                    position.z.rem_euclid(lz),
                )
            })
            .collect()
    }

    fn marker_positions(&self, markers: &[LagrangianPoint<f64>]) -> Vec<(f64, f64, f64)> {
        self.positions_to_tuples(
            &markers
                .iter()
                .map(|marker| marker.position)
                .collect::<Vec<_>>(),
        )
    }

    fn normalize_spectrum(&self, spectrum: &Array3<Complex64>) -> Array3<Complex64> {
        let scale = (self.grid.nx * self.grid.ny * self.grid.nz) as f64;
        spectrum.mapv(|value| value / scale)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn build_field(grid: UniformGrid3D) -> VelocityField<f64> {
        let (lx, ly, lz) = grid.lengths();
        let mut components = Vec::with_capacity(grid.nx * grid.ny * grid.nz);

        for k in 0..grid.nz {
            let z = k as f64 * grid.dz;
            for j in 0..grid.ny {
                let y = j as f64 * grid.dy;
                for i in 0..grid.nx {
                    let x = i as f64 * grid.dx;
                    let u = 1.0 + 0.25 * (2.0 * std::f64::consts::PI * x / lx).cos();
                    let v = -0.5 * (2.0 * std::f64::consts::PI * y / ly).sin();
                    let w = 0.125 * (2.0 * std::f64::consts::PI * z / lz).cos();
                    components.push(Vector3::new(u, v, w));
                }
            }
        }

        VelocityField {
            components,
            dimensions: (grid.nx, grid.ny, grid.nz),
        }
    }

    #[test]
    fn samples_irregular_probe_locations() {
        let grid = UniformGrid3D::new(8, 8, 8, 0.25, 0.25, 0.25).unwrap();
        let coupler = NufftMarkerCoupler3D::new(grid);
        let field = build_field(grid);
        let probes = vec![
            Vector3::new(0.13, 0.27, 0.41),
            Vector3::new(0.61, 0.19, 0.03),
        ];

        let samples = coupler
            .sample_velocity_at_positions(&field, &probes)
            .expect("sampling should succeed");

        for (probe, sample) in probes.iter().zip(samples.iter()) {
            let expected = Vector3::new(
                1.0 + 0.25 * (2.0 * std::f64::consts::PI * probe.x / 2.0).cos(),
                -0.5 * (2.0 * std::f64::consts::PI * probe.y / 2.0).sin(),
                0.125 * (2.0 * std::f64::consts::PI * probe.z / 2.0).cos(),
            );

            assert!((sample.x - expected.x).abs() < 1e-9);
            assert!((sample.y - expected.y).abs() < 1e-9);
            assert!((sample.z - expected.z).abs() < 1e-9);
        }
    }

    #[test]
    fn spreads_marker_forces_to_grid() {
        let grid = UniformGrid3D::new(6, 6, 6, 0.2, 0.2, 0.2).unwrap();
        let coupler = NufftMarkerCoupler3D::new(grid);
        let mut marker = LagrangianPoint::new(Vector3::new(0.37, 0.41, 0.43), 1.0);
        marker.force = Vector3::new(1.0, -0.5, 0.25);

        let field = coupler
            .spread_marker_forces(&[marker])
            .expect("spreading should succeed");

        assert_eq!(field.dimensions, (6, 6, 6));
        assert_eq!(field.components.len(), 216);
        assert!(field.components.iter().any(|value| value.x.is_finite()));
        assert!(field.components.iter().any(|value| value.y.is_finite()));
        assert!(field.components.iter().any(|value| value.z.is_finite()));
    }
}
