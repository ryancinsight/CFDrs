//! Convergence criteria for PISO algorithm
//!
//! # Invariant (Residual Monotonicity)
//!
//! Each PISO corrector step reduces the continuity residual
//! $\|\nabla \cdot \mathbf{u}\|_\infty$. Convergence is declared when
//! the residual drops below the user-specified tolerance $\epsilon > 0$.

use crate::fields::SimulationFields;
use crate::grid::StructuredGrid2D;
use crate::solvers::continuity::max_forward_continuity_residual;
use nalgebra::RealField;
use num_traits::FromPrimitive;

/// Convergence criteria for PISO iterations
#[derive(Debug, Clone)]
pub struct ConvergenceCriteria<T: RealField + Copy> {
    /// Tolerance for velocity residual
    pub velocity_tolerance: T,
    /// Tolerance for pressure residual
    pub pressure_tolerance: T,
    /// Tolerance for continuity residual
    pub continuity_tolerance: T,
    /// Maximum iterations
    pub max_iterations: usize,
}

impl<T: RealField + Copy + FromPrimitive> Default for ConvergenceCriteria<T> {
    fn default() -> Self {
        Self {
            max_iterations: 1000,
            velocity_tolerance: T::from_f64(1e-5)
                .unwrap_or_else(|| T::from_f64(1e-6).expect("analytical constant conversion")),
            pressure_tolerance: T::from_f64(1e-4)
                .unwrap_or_else(|| T::from_f64(1e-5).expect("analytical constant conversion")),
            continuity_tolerance: T::from_f64(1e-5)
                .unwrap_or_else(|| T::from_f64(1e-6).expect("analytical constant conversion")),
        }
    }
}

/// Convergence monitor for tracking residuals
pub struct ConvergenceMonitor<T: RealField + Copy> {
    /// Velocity residual history
    pub velocity_residuals: Vec<T>,
    /// Pressure residual history
    pub pressure_residuals: Vec<T>,
    /// Continuity residual history
    pub continuity_residuals: Vec<T>,
    /// Current iteration
    pub iteration: usize,
}

impl<T: RealField + Copy + FromPrimitive + Copy> Default for ConvergenceMonitor<T> {
    fn default() -> Self {
        Self::new()
    }
}

impl<T: RealField + Copy + FromPrimitive + Copy> ConvergenceMonitor<T> {
    /// Create new convergence monitor
    #[must_use]
    pub fn new() -> Self {
        Self {
            velocity_residuals: Vec::new(),
            pressure_residuals: Vec::new(),
            continuity_residuals: Vec::new(),
            iteration: 0,
        }
    }

    /// Check if converged
    pub fn is_converged(&self, criteria: &ConvergenceCriteria<T>) -> bool {
        if self.velocity_residuals.is_empty()
            || self.pressure_residuals.is_empty()
            || self.continuity_residuals.is_empty()
        {
            return false;
        }

        let vel_res = self.velocity_residuals.last().copied().unwrap_or(T::one());
        let pres_res = self.pressure_residuals.last().copied().unwrap_or(T::one());
        let cont_res = self
            .continuity_residuals
            .last()
            .copied()
            .unwrap_or(T::one());

        vel_res < criteria.velocity_tolerance
            && pres_res < criteria.pressure_tolerance
            && cont_res < criteria.continuity_tolerance
    }

    /// Update residuals
    pub fn update(
        &mut self,
        fields_old: &SimulationFields<T>,
        fields_new: &SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) {
        let vel_res = self.calculate_velocity_residual(fields_old, fields_new, grid.nx, grid.ny);
        let pres_res = self.calculate_pressure_residual(fields_old, fields_new, grid.nx, grid.ny);
        let cont_res = self.calculate_continuity_residual(fields_new, grid);

        self.velocity_residuals.push(vel_res);
        self.pressure_residuals.push(pres_res);
        self.continuity_residuals.push(cont_res);
        self.iteration += 1;
    }

    /// Calculate velocity residual (L2 norm)
    fn calculate_velocity_residual(
        &self,
        fields_old: &SimulationFields<T>,
        fields_new: &SimulationFields<T>,
        nx: usize,
        ny: usize,
    ) -> T {
        let mut sum = T::zero();
        let mut count = 0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let du = fields_new.u.at(i, j) - fields_old.u.at(i, j);
                let dv = fields_new.v.at(i, j) - fields_old.v.at(i, j);
                sum = sum + du * du + dv * dv;
                count += 2;
            }
        }

        let count_t = T::from_usize(count).expect("analytical constant conversion");
        (sum / count_t).sqrt()
    }

    /// Calculate pressure residual (L2 norm)
    fn calculate_pressure_residual(
        &self,
        fields_old: &SimulationFields<T>,
        fields_new: &SimulationFields<T>,
        nx: usize,
        ny: usize,
    ) -> T {
        let mut sum = T::zero();
        let mut count = 0;

        for i in 1..nx - 1 {
            for j in 1..ny - 1 {
                let dp = fields_new.p.at(i, j) - fields_old.p.at(i, j);
                sum += dp * dp;
                count += 1;
            }
        }

        let count_t = T::from_usize(count).expect("analytical constant conversion");
        (sum / count_t).sqrt()
    }

    /// Calculate continuity residual (mass imbalance)
    fn calculate_continuity_residual(
        &self,
        fields: &SimulationFields<T>,
        grid: &StructuredGrid2D<T>,
    ) -> T {
        max_forward_continuity_residual(
            grid.nx,
            grid.ny,
            grid.dx,
            grid.dy,
            |i, j| fields.u.at(i, j),
            |i, j| fields.v.at(i, j),
        )
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn continuity_residual_uses_grid_spacing() {
        let mut monitor = ConvergenceMonitor::<f64>::new();
        let fields_old = SimulationFields::new(4, 4);
        let mut fields_new = SimulationFields::new(4, 4);

        let dx = 0.5_f64;
        let dy = 2.0_f64;

        for j in 0..4 {
            for i in 0..4 {
                fields_new.u.set(i, j, 2.0 * (i as f64) * dx);
                fields_new.v.set(i, j, 3.0 * (j as f64) * dy);
            }
        }

        let grid = StructuredGrid2D::new(4, 4, 0.0, 1.5, 0.0, 6.0).unwrap();

        monitor.update(&fields_old, &fields_new, &grid);

        assert!((monitor.continuity_residuals[0] - 5.0).abs() < 1e-12);
        assert_eq!(monitor.iteration, 1);
    }
}
