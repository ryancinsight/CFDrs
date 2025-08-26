//! Boundary condition handling for LBM.
//!
//! This module implements various boundary conditions including
//! bounce-back, velocity, and pressure boundaries.

use crate::solvers::lbm::lattice::{equilibrium, D2Q9};
use cfd_core::BoundaryCondition;
use nalgebra::{RealField, Vector2};
use num_traits::FromPrimitive;
use std::collections::HashMap;
/// Types of LBM boundary conditions
#[derive(Debug, Clone, Copy, PartialEq)]
pub enum BoundaryType {
    /// Bounce-back (no-slip wall)
    BounceBack,
    /// Velocity boundary (Dirichlet)
    Velocity,
    /// Pressure boundary (Dirichlet)
    Pressure,
    /// Open boundary (Neumann)
    Open,
    /// Periodic boundary
    Periodic,
}
/// Boundary handler for applying boundary conditions
pub struct BoundaryHandler<T: RealField + Copy> {
    /// Boundary type for each edge
    boundary_types: HashMap<String, BoundaryType>,
    /// Phantom data for type parameter
    _phantom: std::marker::PhantomData<T>,
impl<T: RealField + Copy + FromPrimitive> Default for BoundaryHandler<T> {
    fn default() -> Self {
        Self::new()
    }
impl<T: RealField + Copy + FromPrimitive> BoundaryHandler<T> {
    /// Create new boundary handler
    #[must_use]
    pub fn new() -> Self {
        Self {
            boundary_types: HashMap::new(),
            _phantom: std::marker::PhantomData,
        }
    /// Set boundary type for an edge
    pub fn set_boundary(&mut self, edge: String, boundary_type: BoundaryType) {
        self.boundary_types.insert(edge, boundary_type);
    /// Apply bounce-back boundary condition
    pub fn apply_bounce_back(f: &mut Vec<Vec<[T; 9]>>, i: usize, j: usize) {
        let f_current = f[j][i];
        for q in 0..9 {
            let q_opp = D2Q9::OPPOSITE[q];
            f[j][i][q] = f_current[q_opp];
    /// Apply velocity boundary condition (Zou-He)
    pub fn apply_velocity_boundary(
        f: &mut Vec<Vec<[T; 9]>>,
        density: &mut Vec<Vec<T>>,
        velocity: &mut Vec<Vec<[T; 2]>>,
        i: usize,
        j: usize,
        u_boundary: Vector2<T>,
    ) {
        // Set velocity
        velocity[j][i][0] = u_boundary.x;
        velocity[j][i][1] = u_boundary.y;
        // Compute density from known distributions
        let rho = Self::compute_boundary_density(f, i, j);
        density[j][i] = rho;
        // Compute equilibrium distributions
            let weight = T::from_f64(D2Q9::WEIGHTS[q]).unwrap_or_else(T::zero);
            let lattice_vel = &D2Q9::VELOCITIES[q];
            let u_arr = [u_boundary.x, u_boundary.y];
            f[j][i][q] = equilibrium(rho, &u_arr, q, weight, lattice_vel);
    /// Apply pressure boundary condition
    pub fn apply_pressure_boundary(
        p_boundary: T,
        // Convert pressure to density (assuming cs^2 = 1/3)
        let cs2 = T::from_f64(1.0 / 3.0).unwrap_or_else(T::zero);
        let rho = p_boundary / cs2;
        // Extrapolate velocity from interior
        let u = Self::extrapolate_velocity(velocity, i, j);
        velocity[j][i] = u;
            f[j][i][q] = equilibrium(rho, &u, q, weight, lattice_vel);
    /// Apply all boundary conditions
    pub fn apply_boundaries(
        &self,
        boundaries: &HashMap<(usize, usize), BoundaryCondition<T>>,
        for ((i, j), bc) in boundaries {
            match bc {
                BoundaryCondition::Wall { .. } => {
                    Self::apply_bounce_back(f, *i, *j);
                }
                BoundaryCondition::VelocityInlet { velocity: vel } => {
                    let u_boundary = Vector2::new(vel[0], vel[1]);
                    Self::apply_velocity_boundary(f, density, velocity, *i, *j, u_boundary);
                BoundaryCondition::PressureInlet { pressure }
                | BoundaryCondition::PressureOutlet { pressure } => {
                    Self::apply_pressure_boundary(f, density, velocity, *i, *j, *pressure);
                _ => {}
            }
    /// Compute density at boundary from known distributions
    fn compute_boundary_density(f: &Vec<Vec<[T; 9]>>, i: usize, j: usize) -> T {
        let mut rho = T::zero();
            rho += f[j][i][q];
        rho
    /// Extrapolate velocity from interior points
    fn extrapolate_velocity(velocity: &Vec<Vec<[T; 2]>>, i: usize, j: usize) -> [T; 2] {
        let ny = velocity.len();
        let nx = if ny > 0 { velocity[0].len() } else { 0 };
        // First-order extrapolation boundary treatment
        if i == 0 && i + 1 < nx {
            velocity[j][i + 1]
        } else if i == nx - 1 && i > 0 {
            velocity[j][i - 1]
        } else if j == 0 && j + 1 < ny {
            velocity[j + 1][i]
        } else if j == ny - 1 && j > 0 {
            velocity[j - 1][i]
        } else {
            [T::zero(), T::zero()]
#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn test_bounce_back() {
        let mut f = vec![vec![[0.1_f64; 9]; 10]; 10];
        // Set specific values
            f[5][5][q] = q as f64;
        BoundaryHandler::<f64>::apply_bounce_back(&mut f, 5, 5);
        // Check that distributions are swapped
            assert_eq!(f[5][5][q], q_opp as f64);
