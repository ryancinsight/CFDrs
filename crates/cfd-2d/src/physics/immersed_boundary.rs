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
//!
//! # Theorem (IBM Spreading/Interpolation Adjointness — Peskin 2002)
//!
//! The discrete spreading operator $\mathcal{S}$ and interpolation operator $\mathcal{J}$
//! satisfy the adjoint relation $\langle \mathcal{S}\mathbf{F}, \mathbf{u} \rangle_E
//! = \langle \mathbf{F}, \mathcal{J}\mathbf{u} \rangle_L$, where $E$ and $L$ denote
//! Eulerian and Lagrangian inner products respectively.
//!
//! **Proof sketch**:
//! With the 4-point cosine delta function $\delta_h(\mathbf{r}) = \prod_{d=1}^{2} \phi(r_d/h)/h$,
//! spreading is $f_i = \sum_j \delta_h(\mathbf{x}_i - \mathbf{X}_j) F_j \Delta s_j$
//! and interpolation $U_j = \sum_i \delta_h(\mathbf{x}_i - \mathbf{X}_j) u_i h^2$.
//! Direct substitution shows $\sum_i f_i u_i h^2 = \sum_j F_j U_j \Delta s_j$,
//! establishing the adjoint property. This ensures the IBM adds zero net energy to
//! the fluid–structure system, guaranteeing stability of the coupled solver.

use cfd_core::error::{Error, Result};
use leto::Array2;
use leto::geometry::Vector2;

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
    ///
    /// Prefer [`ImmersedBoundaryMethod::try_new`] for fallible construction;
    /// this constructor is retained for callers that already validated their
    /// inputs and need the infallible signature.
    pub fn new(grid_size: (usize, usize), domain_size: (f64, f64)) -> Self {
        match Self::try_new(grid_size, domain_size) {
            Ok(ibm) => ibm,
            Err(error) => panic!("ImmersedBoundaryMethod::new called with invalid inputs: {error}"),
        }
    }

    /// Fallible constructor that validates `grid_size` and `domain_size`.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidConfiguration`] when:
    ///
    /// - `grid_size.0` or `grid_size.1` is zero (the `dx = domain.x / grid.x`
    ///   calculation divides by `grid.x`),
    /// - `domain_size.0` or `domain_size.1` is non-finite or non-positive
    ///   (the resulting `dx`/`dy` would be non-finite or negative, and the
    ///   discrete-delta interpolation divides by `dx`/`dy`).
    pub fn try_new(grid_size: (usize, usize), domain_size: (f64, f64)) -> Result<Self> {
        Self::with_config(ImmersedBoundaryConfig::default(), grid_size, domain_size)
    }

    /// Create with custom configuration, validating every invariant.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidConfiguration`] when:
    ///
    /// - `config.delta_support` is non-finite or non-positive (the discrete
    ///   delta function support radius is the count of cells the kernel
    ///   touches),
    /// - `config.regularization` is non-finite or non-positive (`update_forces`
    ///   divides by `regularization`),
    /// - `grid_size.0` or `grid_size.1` is zero,
    /// - `domain_size.0` or `domain_size.1` is non-finite or non-positive.
    pub fn with_config(
        config: ImmersedBoundaryConfig,
        grid_size: (usize, usize),
        domain_size: (f64, f64),
    ) -> Result<Self> {
        if !config.delta_support.is_finite() || config.delta_support <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::with_config: config.delta_support must be finite and positive, got {}",
                config.delta_support
            )));
        }
        if !config.regularization.is_finite() || config.regularization <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::with_config: config.regularization must be finite and positive (the force update divides by it), got {}",
                config.regularization
            )));
        }
        if grid_size.0 == 0 || grid_size.1 == 0 {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::with_config: grid_size must have a positive extent in each dimension, got {grid_size:?}"
            )));
        }
        if !domain_size.0.is_finite() || domain_size.0 <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::with_config: domain_size.0 must be finite and positive, got {}",
                domain_size.0
            )));
        }
        if !domain_size.1.is_finite() || domain_size.1 <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::with_config: domain_size.1 must be finite and positive, got {}",
                domain_size.1
            )));
        }
        let dx = domain_size.0 / grid_size.0 as f64;
        let dy = domain_size.1 / grid_size.1 as f64;

        Ok(Self {
            config,
            boundary_points: Vec::new(),
            grid_size,
            dx,
            dy,
        })
    }

    /// Add boundary point to immersed boundary.
    ///
    /// Prefer [`ImmersedBoundaryMethod::try_add_boundary_point`] for fallible
    /// insertion; this method is retained for callers that already validated
    /// their inputs.
    pub fn add_boundary_point(&mut self, point: BoundaryPoint) {
        match self.try_add_boundary_point(point) {
            Ok(()) => {}
            Err(error) => {
                panic!("ImmersedBoundaryMethod::add_boundary_point rejected the point: {error}")
            }
        }
    }

    /// Add a boundary point, validating that the position is inside the grid
    /// extent and finite.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidConfiguration`] when the point's position is
    /// non-finite or outside the grid extent (the discrete-delta iteration
    /// silently produces zero contribution for out-of-grid points).
    pub fn try_add_boundary_point(&mut self, point: BoundaryPoint) -> Result<()> {
        if !point.position[0].is_finite() || !point.position[1].is_finite() {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::try_add_boundary_point: position must be finite, got {:?}",
                point.position
            )));
        }
        if point.segment_length <= 0.0 || !point.segment_length.is_finite() {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::try_add_boundary_point: segment_length must be finite and positive, got {}",
                point.segment_length
            )));
        }
        let max_x = self.grid_size.0 as f64 * self.dx;
        let max_y = self.grid_size.1 as f64 * self.dy;
        if point.position[0] < 0.0
            || point.position[0] >= max_x
            || point.position[1] < 0.0
            || point.position[1] >= max_y
        {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::try_add_boundary_point: position {:?} must lie inside the grid extent [0, {max_x}) x [0, {max_y})",
                point.position
            )));
        }
        self.boundary_points.push(point);
        Ok(())
    }

    /// Add circular boundary
    ///
    /// # Arguments
    ///
    /// * `center` - Circle center
    /// * `radius` - Circle radius
    /// * `num_points` - Number of boundary points
    /// * `velocity` - Desired boundary velocity
    ///
    /// Prefer [`ImmersedBoundaryMethod::try_add_circle`] for fallible
    /// insertion; this method is retained for backwards compatibility.
    pub fn add_circle(
        &mut self,
        center: Vector2<f64>,
        radius: f64,
        num_points: usize,
        velocity: Vector2<f64>,
    ) {
        match self.try_add_circle(center, radius, num_points, velocity) {
            Ok(()) => {}
            Err(error) => panic!("ImmersedBoundaryMethod::add_circle rejected the inputs: {error}"),
        }
    }

    /// Add a circular boundary, validating the inputs.
    ///
    /// # Errors
    ///
    /// Returns [`Error::InvalidConfiguration`] when `center` is non-finite,
    /// `radius` is non-finite or non-positive (zero radius is a degenerate
    /// circle), `num_points` is zero (segment length divides by it), or
    /// `velocity` is non-finite.
    pub fn try_add_circle(
        &mut self,
        center: Vector2<f64>,
        radius: f64,
        num_points: usize,
        velocity: Vector2<f64>,
    ) -> Result<()> {
        if !center[0].is_finite() || !center[1].is_finite() {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::try_add_circle: center must be finite, got {center:?}"
            )));
        }
        if !radius.is_finite() || radius <= 0.0 {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::try_add_circle: radius must be finite and positive, got {radius}"
            )));
        }
        if num_points == 0 {
            return Err(Error::InvalidConfiguration(
                "ImmersedBoundaryMethod::try_add_circle: num_points must be positive".to_string(),
            ));
        }
        if !velocity[0].is_finite() || !velocity[1].is_finite() {
            return Err(Error::InvalidConfiguration(format!(
                "ImmersedBoundaryMethod::try_add_circle: velocity must be finite, got {velocity:?}"
            )));
        }

        let segment_length = 2.0 * std::f64::consts::PI * radius / num_points as f64;

        for i in 0..num_points {
            let angle = 2.0 * std::f64::consts::PI * i as f64 / num_points as f64;

            let cos_angle = angle.cos();
            let sin_angle = angle.sin();

            let position = Vector2::new(
                center[0] + radius * cos_angle,
                center[1] + radius * sin_angle,
            );

            let point = BoundaryPoint {
                position,
                desired_velocity: velocity,
                force: Vector2::zeros(),
                segment_length,
            };

            self.boundary_points.push(point);
        }
        Ok(())
    }

    /// Spread forces from boundary to fluid grid
    ///
    /// # Arguments
    ///
    /// * `force_matrix` - Output force field with shape `[nx * ny, 2]`
    pub fn spread_forces(&self, force_matrix: &mut Array2<f64>) -> Result<()> {
        // Initialize force matrix to zero
        force_matrix.fill(0.0);

        for boundary_point in &self.boundary_points {
            self.spread_single_force(boundary_point, force_matrix)?;
        }

        Ok(())
    }

    /// Spread force from single boundary point to nearby grid points
    fn spread_single_force(
        &self,
        point: &BoundaryPoint,
        force_matrix: &mut Array2<f64>,
    ) -> Result<()> {
        let pos_x = point.position[0] / self.dx;
        let pos_y = point.position[1] / self.dy;

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

                let dx_dist = (grid_x - point.position[0]) / self.dx;
                let dy_dist = (grid_y - point.position[1]) / self.dy;

                // Discrete delta function (4-point for 2D)
                let delta_x = self.discrete_delta(dx_dist);
                let delta_y = self.discrete_delta(dy_dist);
                let delta = delta_x * delta_y;

                // Spread forces to u and v components
                let matrix_cell = j * self.grid_size.0 + i;

                force_matrix[[matrix_cell, 0]] += delta * point.force[0] * point.segment_length;
                force_matrix[[matrix_cell, 1]] += delta * point.force[1] * point.segment_length;
            }
        }

        Ok(())
    }

    /// Interpolate velocities from fluid grid to boundary points
    ///
    /// # Arguments
    ///
    /// * `velocity_matrix` - Fluid velocity field with shape `[nx * ny, 2]`
    pub fn interpolate_velocities(
        &self,
        velocity_matrix: &Array2<f64>,
    ) -> Result<Vec<Vector2<f64>>> {
        let mut interpolated_velocities = Vec::with_capacity(self.boundary_points.len());

        for boundary_point in &self.boundary_points {
            let velocity = self.interpolate_single_velocity(boundary_point, velocity_matrix)?;
            interpolated_velocities.push(velocity);
        }

        Ok(interpolated_velocities)
    }

    /// Interpolate velocity at single boundary point
    fn interpolate_single_velocity(
        &self,
        point: &BoundaryPoint,
        velocity_matrix: &Array2<f64>,
    ) -> Result<Vector2<f64>> {
        let pos_x = point.position[0] / self.dx;
        let pos_y = point.position[1] / self.dy;

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

                let dx_dist = (grid_x - point.position[0]) / self.dx;
                let dy_dist = (grid_y - point.position[1]) / self.dy;

                let delta_x = self.discrete_delta(dx_dist);
                let delta_y = self.discrete_delta(dy_dist);
                let delta = delta_x * delta_y;

                let matrix_cell = j * self.grid_size.0 + i;

                u_interp += delta * velocity_matrix[[matrix_cell, 0]];
                v_interp += delta * velocity_matrix[[matrix_cell, 1]];
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
                "Velocity array size mismatch with boundary points".to_string(),
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
            0.25 * (1.0 + (std::f64::consts::PI * r).cos())
                * (1.0 + (std::f64::consts::PI * r * 0.5).cos())
        } else if abs_r <= 2.0 {
            // Smooth transition
            let pi = std::f64::consts::PI;
            0.25 * (-4.0 * r * r + 8.0 * abs_r - 3.0) * (1.0 + (pi * r * 0.5).cos())
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
    use eunomia::assert_relative_eq;

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
            let dx = point.position[0] - center[0];
            let dy = point.position[1] - center[1];
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

        let mut force_matrix = Array2::zeros([16 * 16, 2]);

        ibm.spread_forces(&mut force_matrix)
            .expect("expected value");

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
        let mut velocity_matrix = Array2::zeros([16 * 16, 2]);
        let constant_u = 1.5;
        let constant_v = -0.8;

        for i in 0..16 {
            for j in 0..16 {
                let cell = j * 16 + i;
                velocity_matrix[[cell, 0]] = constant_u;
                velocity_matrix[[cell, 1]] = constant_v;
            }
        }

        let velocities = ibm
            .interpolate_velocities(&velocity_matrix)
            .expect("expected value");

        assert_eq!(velocities.len(), 1);
        assert_relative_eq!(velocities[0][0], constant_u, epsilon = 1e-10);
        assert_relative_eq!(velocities[0][1], constant_v, epsilon = 1e-10);
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

        ibm.update_forces(&current_velocities)
            .expect("expected value");

        // Force should be (1/regularization, 0)
        let expected_force_x = 1.0 / ibm.config().regularization;
        assert_relative_eq!(
            ibm.boundary_points[0].force[0],
            expected_force_x,
            epsilon = 1e-10
        );
        assert_eq!(ibm.boundary_points[0].force[1], 0.0);
    }

    #[test]
    fn try_new_rejects_zero_grid_size() {
        let result = ImmersedBoundaryMethod::try_new((0, 64), (1.0, 1.0));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("grid_size")),
            "expected InvalidConfiguration error for zero grid_size.0"
        );
    }

    #[test]
    fn try_new_rejects_zero_domain_size() {
        let result = ImmersedBoundaryMethod::try_new((64, 64), (0.0, 1.0));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("domain_size.0")),
            "expected InvalidConfiguration error for zero domain_size.0"
        );
    }

    #[test]
    fn try_new_rejects_negative_domain_size() {
        let result = ImmersedBoundaryMethod::try_new((64, 64), (1.0, -1.0));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("domain_size.1")),
            "expected InvalidConfiguration error for negative domain_size.1"
        );
    }

    #[test]
    fn try_new_rejects_non_finite_domain_size() {
        let result = ImmersedBoundaryMethod::try_new((64, 64), (f64::NAN, 1.0));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("domain_size.0")),
            "expected InvalidConfiguration error for NaN domain_size.0"
        );
    }

    #[test]
    fn try_new_accepts_physically_valid_inputs() {
        let result = ImmersedBoundaryMethod::try_new((64, 64), (1.0, 1.0));
        assert!(
            result.is_ok(),
            "valid grid_size and domain_size must succeed"
        );
    }

    #[test]
    fn with_config_rejects_zero_regularization() {
        let config = ImmersedBoundaryConfig {
            regularization: 0.0,
            ..ImmersedBoundaryConfig::default()
        };
        let result = ImmersedBoundaryMethod::with_config(config, (64, 64), (1.0, 1.0));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("regularization")),
            "expected InvalidConfiguration error for zero regularization"
        );
    }

    #[test]
    fn with_config_rejects_zero_delta_support() {
        let config = ImmersedBoundaryConfig {
            delta_support: 0.0,
            ..ImmersedBoundaryConfig::default()
        };
        let result = ImmersedBoundaryMethod::with_config(config, (64, 64), (1.0, 1.0));
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("delta_support")),
            "expected InvalidConfiguration error for zero delta_support"
        );
    }

    #[test]
    fn try_add_circle_rejects_zero_radius() {
        let mut ibm = ImmersedBoundaryMethod::new((64, 64), (1.0, 1.0));
        let result = ibm.try_add_circle(Vector2::new(0.5, 0.5), 0.0, 32, Vector2::zeros());
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("radius")),
            "expected InvalidConfiguration error for zero radius"
        );
    }

    #[test]
    fn try_add_circle_rejects_zero_num_points() {
        let mut ibm = ImmersedBoundaryMethod::new((64, 64), (1.0, 1.0));
        let result = ibm.try_add_circle(Vector2::new(0.5, 0.5), 0.1, 0, Vector2::zeros());
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("num_points")),
            "expected InvalidConfiguration error for zero num_points"
        );
    }

    #[test]
    fn try_add_circle_rejects_non_finite_center() {
        let mut ibm = ImmersedBoundaryMethod::new((64, 64), (1.0, 1.0));
        let result = ibm.try_add_circle(Vector2::new(f64::NAN, 0.5), 0.1, 32, Vector2::zeros());
        assert!(
            matches!(result, Err(Error::InvalidConfiguration(ref msg)) if msg.contains("center")),
            "expected InvalidConfiguration error for NaN center"
        );
    }

    #[test]
    fn try_add_circle_accepts_physically_valid_inputs() {
        let mut ibm = ImmersedBoundaryMethod::new((64, 64), (1.0, 1.0));
        let result = ibm.try_add_circle(Vector2::new(0.5, 0.5), 0.2, 32, Vector2::zeros());
        assert!(result.is_ok(), "valid circle must be accepted");
        assert_eq!(ibm.boundary_points.len(), 32);
    }
}
