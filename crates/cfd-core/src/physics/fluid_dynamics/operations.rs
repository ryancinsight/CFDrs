//! Flow field operations
//!
//! Provides operations on flow fields including vorticity, divergence,
//! and other fluid mechanical quantities.

use eunomia::{NumericElement, RealField};
use leto::geometry::Vector3;
use moirai::prelude::{ParallelSlice, ParallelSliceMut};

use super::fields::VelocityField;

/// Operations on flow fields
pub struct FlowOperations;

impl FlowOperations {
    /// Calculate vorticity field (curl of velocity)
    #[must_use]
    pub fn vorticity<T: RealField + NumericElement + Copy + Send + Sync>(
        velocity: &VelocityField<T>,
    ) -> Vec<Vector3<T>> {
        let (nx, ny, nz) = velocity.dimensions;
        let mut vorticity = vec![Vector3::zeros(); nx * ny * nz];

        // Use parallel iteration for better performance
        vorticity.par_mut().enumerate(|idx, vort| {
            let k = idx / (nx * ny);
            let j = (idx % (nx * ny)) / nx;
            let i = idx % nx;

            *vort = vorticity_at_point(velocity, i, j, k, nx, ny, nz);
        });

        vorticity
    }

    /// Calculate divergence field
    #[must_use]
    pub fn divergence<T: RealField + NumericElement + Copy + Send + Sync>(
        velocity: &VelocityField<T>,
    ) -> Vec<T> {
        let (nx, ny, nz) = velocity.dimensions;
        let mut divergence = vec![<T as NumericElement>::ZERO; nx * ny * nz];
        let two = two::<T>();

        // Use parallel iteration for better performance
        divergence.par_mut().enumerate(|idx, div| {
            let k = idx / (nx * ny);
            let j = (idx % (nx * ny)) / nx;
            let i = idx % nx;

            // Central differences with boundary handling
            let dudx = if i > 0 && i < nx - 1 {
                let idx_plus = k * nx * ny + j * nx + (i + 1);
                let idx_minus = k * nx * ny + j * nx + (i - 1);
                (velocity.components[idx_plus].x - velocity.components[idx_minus].x) / two
            } else {
                <T as NumericElement>::ZERO
            };

            let dvdy = if j > 0 && j < ny - 1 {
                let idx_plus = k * nx * ny + (j + 1) * nx + i;
                let idx_minus = k * nx * ny + (j - 1) * nx + i;
                (velocity.components[idx_plus].y - velocity.components[idx_minus].y) / two
            } else {
                <T as NumericElement>::ZERO
            };

            let dwdz = if k > 0 && k < nz - 1 {
                let idx_plus = (k + 1) * nx * ny + j * nx + i;
                let idx_minus = (k - 1) * nx * ny + j * nx + i;
                (velocity.components[idx_plus].z - velocity.components[idx_minus].z) / two
            } else {
                <T as NumericElement>::ZERO
            };

            *div = dudx + dvdy + dwdz;
        });

        divergence
    }

    /// Calculate kinetic energy field
    #[must_use]
    pub fn kinetic_energy<T: RealField + NumericElement + Copy + Send + Sync>(
        velocity: &VelocityField<T>,
    ) -> Vec<T> {
        let half = half::<T>();

        // Use parallel iteration for better performance
        velocity
            .components
            .par()
            .map_collect(|v| half * v.norm_squared())
    }

    /// Calculate enstrophy field
    #[must_use]
    pub fn enstrophy<T: RealField + NumericElement + Copy + Send + Sync>(
        velocity: &VelocityField<T>,
    ) -> Vec<T> {
        let vorticity = Self::vorticity(velocity);
        let half = half::<T>();

        // Use parallel iteration for better performance
        vorticity.par().map_collect(|w| half * w.norm_squared())
    }
}

#[inline]
fn two<T: NumericElement>() -> T {
    <T as NumericElement>::ONE + <T as NumericElement>::ONE
}

#[inline]
fn half<T: NumericElement>() -> T {
    <T as NumericElement>::ONE / two::<T>()
}

// Helper function for vorticity calculation at a point
#[allow(clippy::similar_names)] // CFD derivatives use standard notation: dudx, dvdy, etc.
fn vorticity_at_point<T: RealField + NumericElement + Copy>(
    velocity: &VelocityField<T>,
    i: usize,
    j: usize,
    k: usize,
    nx: usize,
    ny: usize,
    nz: usize,
) -> Vector3<T> {
    let two = two::<T>();

    // Calculate velocity gradients using central differences
    let dvdz = if k > 0 && k < nz - 1 {
        let idx_plus = (k + 1) * nx * ny + j * nx + i;
        let idx_minus = (k - 1) * nx * ny + j * nx + i;
        (velocity.components[idx_plus].y - velocity.components[idx_minus].y) / two
    } else {
        <T as NumericElement>::ZERO
    };

    let dwdy = if j > 0 && j < ny - 1 {
        let idx_plus = k * nx * ny + (j + 1) * nx + i;
        let idx_minus = k * nx * ny + (j - 1) * nx + i;
        (velocity.components[idx_plus].z - velocity.components[idx_minus].z) / two
    } else {
        <T as NumericElement>::ZERO
    };

    let dudz = if k > 0 && k < nz - 1 {
        let idx_plus = (k + 1) * nx * ny + j * nx + i;
        let idx_minus = (k - 1) * nx * ny + j * nx + i;
        (velocity.components[idx_plus].x - velocity.components[idx_minus].x) / two
    } else {
        <T as NumericElement>::ZERO
    };

    let dwdx = if i > 0 && i < nx - 1 {
        let idx_plus = k * nx * ny + j * nx + (i + 1);
        let idx_minus = k * nx * ny + j * nx + (i - 1);
        (velocity.components[idx_plus].z - velocity.components[idx_minus].z) / two
    } else {
        <T as NumericElement>::ZERO
    };

    let dvdx = if i > 0 && i < nx - 1 {
        let idx_plus = k * nx * ny + j * nx + (i + 1);
        let idx_minus = k * nx * ny + j * nx + (i - 1);
        (velocity.components[idx_plus].y - velocity.components[idx_minus].y) / two
    } else {
        <T as NumericElement>::ZERO
    };

    let dudy = if j > 0 && j < ny - 1 {
        let idx_plus = k * nx * ny + (j + 1) * nx + i;
        let idx_minus = k * nx * ny + (j - 1) * nx + i;
        (velocity.components[idx_plus].x - velocity.components[idx_minus].x) / two
    } else {
        <T as NumericElement>::ZERO
    };

    // Vorticity = curl(velocity)
    Vector3::new(dwdy - dvdz, dudz - dwdx, dvdx - dudy)
}

#[cfg(test)]
mod tests {
    use super::{FlowOperations, VelocityField};
    use leto::geometry::Vector3;

    fn velocity_field(
        dimensions: (usize, usize, usize),
        component: impl Fn(usize, usize, usize) -> Vector3<f64>,
    ) -> VelocityField<f64> {
        let (nx, ny, nz) = dimensions;
        let mut components = Vec::with_capacity(nx * ny * nz);
        for k in 0..nz {
            for j in 0..ny {
                for i in 0..nx {
                    components.push(component(i, j, k));
                }
            }
        }
        VelocityField {
            components,
            dimensions,
        }
    }

    #[test]
    fn constant_velocity_has_zero_divergence_and_vorticity() {
        let velocity = velocity_field((3, 3, 3), |_, _, _| Vector3::new(2.0, -1.0, 0.5));

        let divergence = FlowOperations::divergence(&velocity);
        let vorticity = FlowOperations::vorticity(&velocity);

        assert!(divergence.iter().all(|value| value.abs() < 1e-12));
        assert!(vorticity.iter().all(|value| value.norm() < 1e-12));
    }

    #[test]
    fn affine_velocity_divergence_matches_central_difference() {
        let velocity = velocity_field((3, 3, 3), |i, j, k| {
            Vector3::new(i as f64, j as f64, k as f64)
        });

        let divergence = FlowOperations::divergence(&velocity);
        let center_idx = 13;

        assert!((divergence[center_idx] - 3.0).abs() < 1e-12);
        assert_eq!(divergence[0], 0.0);
    }

    #[test]
    fn kinetic_energy_uses_half_velocity_norm_squared() {
        let velocity = velocity_field((1, 1, 1), |_, _, _| Vector3::new(2.0, 3.0, 6.0));

        let energy = FlowOperations::kinetic_energy(&velocity);

        assert_eq!(energy, vec![24.5]);
    }
}
