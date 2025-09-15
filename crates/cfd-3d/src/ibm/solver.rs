//! IBM solver implementation

use super::{config::IbmConfig, forcing::*, interpolation::*, lagrangian::LagrangianPoint};
use cfd_core::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::{FromPrimitive, ToPrimitive};

// Feedback control constants
const DEFAULT_PROPORTIONAL_GAIN: f64 = 10.0;
const DEFAULT_INTEGRAL_GAIN: f64 = 1.0;

/// IBM solver for 3D flow around immersed boundaries
pub struct IbmSolver<T: RealField + FromPrimitive + ToPrimitive + Copy> {
    /// Configuration
    config: IbmConfig,
    /// Lagrangian points representing the immersed boundary
    lagrangian_points: Vec<LagrangianPoint<T>>,
    /// Interpolation kernel
    kernel: InterpolationKernel<T>,
    /// Forcing method
    forcing: Box<dyn ForcingMethod<T>>,
    /// Grid spacing
    dx: Vector3<T>,
    /// Grid dimensions
    grid_size: (usize, usize, usize),
}

impl<T: RealField + FromPrimitive + ToPrimitive + Copy> IbmSolver<T> {
    /// Create a new IBM solver
    pub fn new(config: IbmConfig, dx: Vector3<T>, grid_size: (usize, usize, usize)) -> Self {
        let kernel = InterpolationKernel::new(
            DeltaFunction::RomaPeskin4,
            T::from_f64(config.smoothing_width).unwrap_or_else(T::one),
        );

        let forcing: Box<dyn ForcingMethod<T>> = if config.use_direct_forcing {
            Box::new(DirectForcing::new(
                T::from_f64(config.force_scale).unwrap_or_else(T::one),
            ))
        } else {
            Box::new(FeedbackForcing::new(
                T::from_f64(DEFAULT_PROPORTIONAL_GAIN).unwrap_or_else(T::one),
                T::from_f64(DEFAULT_INTEGRAL_GAIN).unwrap_or_else(T::one),
            ))
        };

        Self {
            config,
            lagrangian_points: Vec::new(),
            kernel,
            forcing,
            dx,
            grid_size,
        }
    }

    /// Add a Lagrangian point
    pub fn add_lagrangian_point(&mut self, point: LagrangianPoint<T>) {
        self.lagrangian_points.push(point);
    }

    /// Interpolate velocity from Eulerian to Lagrangian grid
    pub fn interpolate_velocity(&mut self, eulerian_velocity: &[Vector3<T>]) -> Result<()> {
        // Collect velocities first to avoid borrow conflict
        let velocities: Vec<_> = self
            .lagrangian_points
            .iter()
            .map(|point| self.interpolate_at_point(&point.position, eulerian_velocity))
            .collect::<Result<Vec<_>>>()?;

        // Then update the points
        for (point, velocity) in self.lagrangian_points.iter_mut().zip(velocities) {
            point.velocity = velocity;
        }
        Ok(())
    }

    /// Calculate forcing terms
    pub fn calculate_forcing(&mut self, desired_velocity: &[Vector3<T>], dt: T) -> Result<()> {
        for (i, point) in self.lagrangian_points.iter_mut().enumerate() {
            point.force = self
                .forcing
                .calculate_force(&desired_velocity[i], &point.velocity, dt);
        }
        Ok(())
    }

    /// Spread forces from Lagrangian to Eulerian grid
    pub fn spread_forces(&self) -> Result<Vec<Vector3<T>>> {
        let total_points = self.grid_size.0 * self.grid_size.1 * self.grid_size.2;
        let mut eulerian_forces = vec![Vector3::zeros(); total_points];

        for point in &self.lagrangian_points {
            self.spread_from_point(&point, &mut eulerian_forces)?;
        }

        Ok(eulerian_forces)
    }

    /// Interpolate at a specific point
    fn interpolate_at_point(
        &self,
        position: &Vector3<T>,
        field: &[Vector3<T>],
    ) -> Result<Vector3<T>> {
        let mut result = Vector3::zeros();
        let stencil = self.kernel.stencil_size();

        // Find grid indices and convert to integers
        let i_int = ((position.x / self.dx.x).floor()).to_isize().unwrap_or(0);
        let j_int = ((position.y / self.dx.y).floor()).to_isize().unwrap_or(0);
        let k_int = ((position.z / self.dx.z).floor()).to_isize().unwrap_or(0);

        let i_start = i_int - (stencil as isize / 2);
        let j_start = j_int - (stencil as isize / 2);
        let k_start = k_int - (stencil as isize / 2);

        for di in 0..stencil {
            for dj in 0..stencil {
                for dk in 0..stencil {
                    let ii = (i_start + di as isize) as usize;
                    let jj = (j_start + dj as isize) as usize;
                    let kk = (k_start + dk as isize) as usize;

                    if ii < self.grid_size.0 && jj < self.grid_size.1 && kk < self.grid_size.2 {
                        let idx =
                            ii + jj * self.grid_size.0 + kk * self.grid_size.0 * self.grid_size.1;

                        let rx = (position.x / self.dx.x
                            - T::from_usize(ii).unwrap_or_else(T::zero))
                        .abs();
                        let ry = (position.y / self.dx.y
                            - T::from_usize(jj).unwrap_or_else(T::zero))
                        .abs();
                        let rz = (position.z / self.dx.z
                            - T::from_usize(kk).unwrap_or_else(T::zero))
                        .abs();

                        let weight =
                            self.kernel.delta(rx) * self.kernel.delta(ry) * self.kernel.delta(rz);
                        result += field[idx] * weight;
                    }
                }
            }
        }

        Ok(result)
    }

    /// Spread from a Lagrangian point to Eulerian grid
    fn spread_from_point(
        &self,
        point: &LagrangianPoint<T>,
        field: &mut [Vector3<T>],
    ) -> Result<()> {
        let stencil = self.kernel.stencil_size();

        // Find grid indices and convert to integers
        let i_int = ((point.position.x / self.dx.x).floor())
            .to_isize()
            .unwrap_or(0);
        let j_int = ((point.position.y / self.dx.y).floor())
            .to_isize()
            .unwrap_or(0);
        let k_int = ((point.position.z / self.dx.z).floor())
            .to_isize()
            .unwrap_or(0);

        let i_start = i_int - (stencil as isize / 2);
        let j_start = j_int - (stencil as isize / 2);
        let k_start = k_int - (stencil as isize / 2);

        for di in 0..stencil {
            for dj in 0..stencil {
                for dk in 0..stencil {
                    let ii = (i_start + di as isize) as usize;
                    let jj = (j_start + dj as isize) as usize;
                    let kk = (k_start + dk as isize) as usize;

                    if ii < self.grid_size.0 && jj < self.grid_size.1 && kk < self.grid_size.2 {
                        let idx =
                            ii + jj * self.grid_size.0 + kk * self.grid_size.0 * self.grid_size.1;

                        let rx = (point.position.x / self.dx.x
                            - T::from_usize(ii).unwrap_or_else(T::zero))
                        .abs();
                        let ry = (point.position.y / self.dx.y
                            - T::from_usize(jj).unwrap_or_else(T::zero))
                        .abs();
                        let rz = (point.position.z / self.dx.z
                            - T::from_usize(kk).unwrap_or_else(T::zero))
                        .abs();

                        let weight =
                            self.kernel.delta(rx) * self.kernel.delta(ry) * self.kernel.delta(rz);
                        field[idx] += point.force * weight * point.weight;
                    }
                }
            }
        }

        Ok(())
    }

    /// Update Lagrangian point positions
    pub fn update_positions(&mut self, dt: T) {
        for point in &mut self.lagrangian_points {
            point.update_position(dt);
        }
    }

    /// Get number of Lagrangian points
    pub fn num_points(&self) -> usize {
        self.lagrangian_points.len()
    }

    /// Get reference to Lagrangian points
    pub fn points(&self) -> &[LagrangianPoint<T>] {
        &self.lagrangian_points
    }

    /// Get mutable reference to Lagrangian points
    pub fn points_mut(&mut self) -> &mut [LagrangianPoint<T>] {
        &mut self.lagrangian_points
    }
}
