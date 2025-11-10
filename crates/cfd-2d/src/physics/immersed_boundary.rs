//! Immersed Boundary Methods for Complex Geometries
//!
//! Immersed boundary methods allow simulation of complex geometries on Cartesian
//! grids without the need for body-conforming meshes, greatly simplifying
//! the mesh generation process for complex geometries.
//!
//! ## Mathematical Foundation
//!
//! The immersed boundary method introduces a force density f(x,t) that enforces
//! the boundary conditions on the immersed boundary Γ:
//!
//! ```math
//! ∂u/∂t + u·∇u = -∇p/ν + ∇²u + f(x,t)
//! ∇·u = 0
//! ```
//!
//! The force f(x,t) is related to the boundary velocity U_b through:
//!
//! ```math
//! f(x,t) = ∫_Γ F(s,t) δ(x - X(s,t)) ds
//! ```
//!
//! where F(s,t) is chosen to enforce U_b = u(X(s,t),t).
//!
//! ## Implementation Strategy
//!
//! 1. **Boundary Representation**: Lagrangian markers on immersed boundary
//! 2. **Force Spreading**: Distribute boundary forces to Eulerian grid
//! 3. **Velocity Interpolation**: Interpolate fluid velocity to boundary points
//! 4. **Force Update**: Compute forces to enforce boundary conditions
//!
//! ## Discrete Formulation
//!
//! **Force Spreading:**
//! ```math
//! f_i = Σ_j D(r_{ij}) F_j Δs_j
//! ```
//!
//! **Velocity Interpolation:**
//! ```math
//! U_j = Σ_i D(r_{ij}) u_i
//! ```
//!
//! where D(r) is a discrete delta function (typically 4-point or 6-point).
//!
//! ## Key Advantages
//!
//! 1. **Geometric Flexibility**: Handle arbitrarily complex geometries
//! 2. **Mesh Simplification**: Use simple Cartesian grids
//! 3. **Time Efficiency**: Avoid costly mesh regeneration
//! 4. **Parallel Friendly**: Natural domain decomposition
//!
//! ## Boundary Conditions
//!
//! The method naturally handles:
//! - **Dirichlet BCs**: Fixed velocity boundaries
//! - **Neumann BCs**: Stress-free boundaries
//! - **Moving Boundaries**: Time-dependent boundary motion
//! - **Elastic Boundaries**: Deformable structures
//!
//! ## Literature Compliance
//!
//! - Peskin, C. S. (1972). Flow patterns around heart valves: A numerical method.
//!   Journal of Computational Physics, 10(2), 252-271.
//! - Peskin, C. S. (2002). The immersed boundary method. Acta Numerica, 11, 479-517.
//! - Mittal, R., & Iaccarino, G. (2005). Immersed boundary methods. Annual Review
//!   of Fluid Mechanics, 37, 239-261.
//! - Griffith, B. E., et al. (2009). An efficient hybrid immersed boundary method.
//!   SIAM Journal on Scientific Computing, 31(6), 3488-3511.
//!
//! ## Numerical Considerations
//!
//! - **Delta Function**: Choice affects accuracy and stability
//! - **Time Stepping**: Must be consistent with fluid solver
//! - **Boundary Resolution**: Sufficient boundary point density required
//! - **Force Regularization**: Prevents numerical instabilities

use nalgebra::{DMatrix, Vector2};
use crate::error::{Error, Result};

/// Immersed boundary configuration
#[derive(Debug, Clone, Copy)]
pub struct ImmersedBoundaryConfig {
    /// Support radius for delta function (typically 2-4 grid cells)
    pub delta_support: f64,
    /// Regularization parameter for force computation
    pub regularization: f64,
    /// Maximum iterations for force computation
    pub max_iterations: usize,
    /// Convergence tolerance for force iteration
    pub tolerance: f64,
}

impl Default for ImmersedBoundaryConfig {
    fn default() -> Self {
        Self {
            delta_support: 2.0,
            regularization: 1e-6,
            max_iterations: 50,
            tolerance: 1e-8,
        }
    }
}

/// Lagrangian boundary point representation
#[derive(Debug, Clone)]
pub struct BoundaryPoint {
    /// Position in physical space
    pub position: Vector2<f64>,
    /// Desired velocity (Dirichlet BC)
    pub desired_velocity: Vector2<f64>,
    /// Force exerted on fluid
    pub force: Vector2<f64>,
    /// Boundary segment length
    pub segment_length: f64,
}

/// Immersed boundary method implementation
#[derive(Debug, Clone)]
pub struct ImmersedBoundaryMethod {
    config: ImmersedBoundaryConfig,
    /// Lagrangian boundary points
    boundary_points: Vec<BoundaryPoint>,
    /// Grid dimensions
    grid_size: (usize, usize),
    /// Grid spacing
    dx: f64,
    dy: f64,
}

impl ImmersedBoundaryMethod {
    /// Create immersed boundary method
    ///
    /// # Arguments
    ///
    /// * `grid_size` - (nx, ny) grid dimensions
    /// * `domain_size` - (Lx, Ly) physical domain size
    ///
    /// # Returns
    ///
    /// Immersed boundary method instance
    pub fn new(grid_size: (usize, usize), domain_size: (f64, f64)) -> Self {
        let dx = domain_size.0 / grid_size.0 as f64;
        let dy = domain_size.1 / grid_size.1 as f64;

        Self {
            config: ImmersedBoundaryConfig::default(),
            boundary_points: Vec::new(),
            grid_size,
            dx,
            dy,
        }
    }

    /// Create with custom configuration
    pub fn with_config(config: ImmersedBoundaryConfig, grid_size: (usize, usize), domain_size: (f64, f64)) -> Self {
        let dx = domain_size.0 / grid_size.0 as f64;
        let dy = domain_size.1 / grid_size.1 as f64;

        Self {
            config,
            boundary_points: Vec::new(),
            grid_size,
            dx,
            dy,
        }
    }

    /// Add boundary point to immersed boundary
    pub fn add_boundary_point(&mut self, point: BoundaryPoint) {
        self.boundary_points.push(point);
    }

    /// Add circular boundary
    ///
    /// # Arguments
    ///
    /// * `center` - Circle center
    /// * `radius` - Circle radius
    /// * `num_points` - Number of boundary points
    /// * `velocity` - Desired boundary velocity
    pub fn add_circle(&mut self, center: Vector2<f64>, radius: f64, num_points: usize, velocity: Vector2<f64>) {
        let segment_length = 2.0 * std::f64::consts::PI * radius / num_points as f64;

        for i in 0..num_points {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / num_points as f64;

            let cos_angle = angle.cos();
            let sin_angle = angle.sin();

            let position = Vector2::new(
                center.x + radius * cos_angle,
                center.y + radius * sin_angle,
            );

            let point = BoundaryPoint {
                position,
                desired_velocity: velocity,
                force: Vector2::zeros(),
                segment_length,
            };

            self.boundary_points.push(point);
        }
    }

    /// Spread forces from boundary to fluid grid
    ///
    /// # Arguments
    ///
    /// * `force_matrix` - Output force field (nx x ny x 2)
    pub fn spread_forces(&self, force_matrix: &mut DMatrix<f64>) -> Result<()> {
        // Initialize force matrix to zero
        force_matrix.fill(0.0);

        for boundary_point in &self.boundary_points {
            self.spread_single_force(boundary_point, force_matrix)?;
        }

        Ok(())
    }

    /// Spread force from single boundary point to nearby grid points
    fn spread_single_force(&self, point: &BoundaryPoint, force_matrix: &mut DMatrix<f64>) -> Result<()> {
        let pos_x = point.position.x / self.dx;
        let pos_y = point.position.y / self.dy;

        // Find grid cell containing the boundary point
        let i_center = (pos_x.floor() as usize).min(self.grid_size.0 - 1);
        let j_center = (pos_y.floor() as usize).min(self.grid_size.1 - 1);

        // Spread force using discrete delta function
        let support = self.config.delta_support.ceil() as usize;

        for di in -(support as i32)..=(support as i32) {
            for dj in -(support as i32)..=(support as i32) {
                let i = (i_center as i32 + di) as usize;
                let j = (j_center as i32 + dj) as usize;

                if i >= self.grid_size.0 || j >= self.grid_size.1 {
                    continue;
                }

                // Distance from boundary point to grid point
                let grid_x = i as f64 * self.dx;
                let grid_y = j as f64 * self.dy;

                let dx_dist = (grid_x - point.position.x) / self.dx;
                let dy_dist = (grid_y - point.position.y) / self.dy;

                // Discrete delta function (4-point for 2D)
                let delta_x = self.discrete_delta(dx_dist);
                let delta_y = self.discrete_delta(dy_dist);
                let delta = delta_x * delta_y;

                // Spread forces to u and v components
                let matrix_idx_u = 2 * (j * self.grid_size.0 + i);
                let matrix_idx_v = 2 * (j * self.grid_size.0 + i) + 1;

                force_matrix[matrix_idx_u] += delta * point.force.x * point.segment_length;
                force_matrix[matrix_idx_v] += delta * point.force.y * point.segment_length;
            }
        }

        Ok(())
    }

    /// Interpolate velocities from fluid grid to boundary points
    ///
    /// # Arguments
    ///
    /// * `velocity_matrix` - Fluid velocity field (nx x ny x 2)
    pub fn interpolate_velocities(&self, velocity_matrix: &DMatrix<f64>) -> Result<Vec<Vector2<f64>>> {
        let mut interpolated_velocities = Vec::with_capacity(self.boundary_points.len());

        for boundary_point in &self.boundary_points {
            let velocity = self.interpolate_single_velocity(boundary_point, velocity_matrix)?;
            interpolated_velocities.push(velocity);
        }

        Ok(interpolated_velocities)
    }

    /// Interpolate velocity at single boundary point
    fn interpolate_single_velocity(&self, point: &BoundaryPoint, velocity_matrix: &DMatrix<f64>) -> Result<Vector2<f64>> {
        let pos_x = point.position.x / self.dx;
        let pos_y = point.position.y / self.dy;

        let i_center = (pos_x.floor() as usize).min(self.grid_size.0 - 1);
        let j_center = (pos_y.floor() as usize).min(self.grid_size.1 - 1);

        let mut u_interp = 0.0;
        let mut v_interp = 0.0;

        let support = self.config.delta_support.ceil() as usize;

        for di in -(support as i32)..=(support as i32) {
            for dj in -(support as i32)..=(support as i32) {
                let i = (i_center as i32 + di) as usize;
                let j = (j_center as i32 + dj) as usize;

                if i >= self.grid_size.0 || j >= self.grid_size.1 {
                    continue;
                }

                let grid_x = i as f64 * self.dx;
                let grid_y = j as f64 * self.dy;

                let dx_dist = (grid_x - point.position.x) / self.dx;
                let dy_dist = (grid_y - point.position.y) / self.dy;

                let delta_x = self.discrete_delta(dx_dist);
                let delta_y = self.discrete_delta(dy_dist);
                let delta = delta_x * delta_y;

                let matrix_idx_u = 2 * (j * self.grid_size.0 + i);
                let matrix_idx_v = 2 * (j * self.grid_size.0 + i) + 1;

                u_interp += delta * velocity_matrix[matrix_idx_u];
                v_interp += delta * velocity_matrix[matrix_idx_v];
            }
        }

        Ok(Vector2::new(u_interp, v_interp))
    }

    /// Update boundary forces to enforce desired velocities
    ///
    /// Uses a simple iterative approach to compute forces
    ///
    /// # Arguments
    ///
    /// * `current_velocities` - Interpolated velocities at boundary points
    pub fn update_forces(&mut self, current_velocities: &[Vector2<f64>]) -> Result<()> {
        if current_velocities.len() != self.boundary_points.len() {
            return Err(Error::InvalidConfiguration(
                "Velocity array size mismatch with boundary points".to_string()
            ));
        }

        // Simple force update: F = (U_desired - U_current) / regularization
        for (i, boundary_point) in self.boundary_points.iter_mut().enumerate() {
            let velocity_error = boundary_point.desired_velocity - current_velocities[i];
            boundary_point.force = velocity_error / self.config.regularization;
        }

        Ok(())
    }

    /// Discrete delta function (4-point in 2D)
    ///
    /// Based on Peskin's discrete delta function for immersed boundaries
    fn discrete_delta(&self, r: f64) -> f64 {
        let abs_r = r.abs();

        if abs_r <= 1.0 {
            // 4-point delta function
            0.25 * (1.0 + (std::f64::consts::PI * r).cos()) *
            (1.0 + (std::f64::consts::PI * r * 0.5).cos())
        } else if abs_r <= 2.0 {
            // Smooth transition
            let pi = std::f64::consts::PI;
            0.25 * (-4.0 * r * r + 8.0 * abs_r - 3.0) *
            (1.0 + (pi * r * 0.5).cos())
        } else {
            0.0
        }
    }

    /// Get boundary points (read-only)
    pub fn boundary_points(&self) -> &[BoundaryPoint] {
        &self.boundary_points
    }

    /// Get configuration
    pub fn config(&self) -> &ImmersedBoundaryConfig {
        &self.config
    }

    /// Set configuration
    pub fn set_config(&mut self, config: ImmersedBoundaryConfig) {
        self.config = config;
    }

    /// Get grid information
    pub fn grid_info(&self) -> ((usize, usize), (f64, f64)) {
        (self.grid_size, (self.dx, self.dy))
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;
    use num_traits::ToPrimitive;

    #[test]
    fn test_immersed_boundary_creation() {
        let ibm = ImmersedBoundaryMethod::new((64, 64), (1.0, 1.0));
        assert_eq!(ibm.grid_size, (64, 64));
        assert_eq!(ibm.boundary_points.len(), 0);
    }

    #[test]
    fn test_circle_boundary_creation() {
        let mut ibm = ImmersedBoundaryMethod::new((64, 64), (1.0, 1.0));

        let center = Vector2::new(0.5, 0.5);
        let velocity = Vector2::new(1.0, 0.0);

        ibm.add_circle(center, 0.2, 32, velocity);

        assert_eq!(ibm.boundary_points.len(), 32);

        // Check that points are on circle
        for point in &ibm.boundary_points {
            let dx = point.position.x - center.x;
            let dy = point.position.y - center.y;
            let radius = (dx * dx + dy * dy).sqrt();
            assert_relative_eq!(radius, 0.2, epsilon = 1e-10);
            assert_eq!(point.desired_velocity, velocity);
        }
    }

    #[test]
    fn test_discrete_delta_function() {
        let ibm = ImmersedBoundaryMethod::new((64, 64), (1.0, 1.0));

        // Test delta function properties
        assert_eq!(ibm.discrete_delta(0.0), 1.0); // δ(0) = 1

        // Test at integer points (should be zero)
        assert_eq!(ibm.discrete_delta(1.0), 0.0);
        assert_eq!(ibm.discrete_delta(-1.0), 0.0);
        assert_eq!(ibm.discrete_delta(2.0), 0.0);
        assert_eq!(ibm.discrete_delta(-2.0), 0.0);

        // Test outside support
        assert_eq!(ibm.discrete_delta(3.0), 0.0);
        assert_eq!(ibm.discrete_delta(-3.0), 0.0);
    }

    #[test]
    fn test_force_spreading() {
        let mut ibm = ImmersedBoundaryMethod::new((16, 16), (1.0, 1.0));

        // Add a single boundary point
        let point = BoundaryPoint {
            position: Vector2::new(0.5, 0.5),
            desired_velocity: Vector2::zeros(),
            force: Vector2::new(1.0, 0.5),
            segment_length: 0.1,
        };
        ibm.add_boundary_point(point);

        let mut force_matrix = DMatrix::zeros(16 * 16 * 2, 1);

        ibm.spread_forces(&mut force_matrix).unwrap();

        // Check that forces were spread (should have non-zero entries)
        let has_nonzero = force_matrix.iter().any(|&x| x != 0.0);
        assert!(has_nonzero, "Forces should be spread to grid points");
    }

    #[test]
    fn test_velocity_interpolation() {
        let mut ibm = ImmersedBoundaryMethod::new((16, 16), (1.0, 1.0));

        // Add boundary point at center
        let point = BoundaryPoint {
            position: Vector2::new(0.5, 0.5),
            desired_velocity: Vector2::zeros(),
            force: Vector2::zeros(),
            segment_length: 0.1,
        };
        ibm.add_boundary_point(point);

        // Create velocity field with constant velocity
        let mut velocity_matrix = DMatrix::zeros(16 * 16 * 2, 1);
        let constant_u = 1.5;
        let constant_v = -0.8;

        for i in 0..16 {
            for j in 0..16 {
                let idx_u = 2 * (j * 16 + i);
                let idx_v = 2 * (j * 16 + i) + 1;
                velocity_matrix[idx_u] = constant_u;
                velocity_matrix[idx_v] = constant_v;
            }
        }

        let velocities = ibm.interpolate_velocities(&velocity_matrix).unwrap();

        assert_eq!(velocities.len(), 1);
        assert_relative_eq!(velocities[0].x, constant_u, epsilon = 1e-10);
        assert_relative_eq!(velocities[0].y, constant_v, epsilon = 1e-10);
    }

    #[test]
    fn test_force_update() {
        let mut ibm = ImmersedBoundaryMethod::new((16, 16), (1.0, 1.0));

        let point = BoundaryPoint {
            position: Vector2::new(0.5, 0.5),
            desired_velocity: Vector2::new(1.0, 0.0),
            force: Vector2::zeros(),
            segment_length: 0.1,
        };
        ibm.add_boundary_point(point);

        // Current velocity is zero, desired is (1, 0)
        let current_velocities = vec![Vector2::new(0.0, 0.0)];

        ibm.update_forces(&current_velocities).unwrap();

        // Force should be (1/regularization, 0)
        let expected_force_x = 1.0 / ibm.config().regularization.to_f64().unwrap();
        assert_relative_eq!(ibm.boundary_points[0].force.x, expected_force_x, epsilon = 1e-10);
        assert_eq!(ibm.boundary_points[0].force.y, 0.0);
    }
}
