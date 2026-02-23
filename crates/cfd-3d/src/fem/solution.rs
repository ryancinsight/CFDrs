//! Solution structure for FEM

use crate::fem::constants;
use nalgebra::{DVector, RealField, Vector3};

/// Solution for 3D incompressible flow
#[derive(Debug, Clone)]
pub struct StokesFlowSolution<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Velocity field stored in block order: [u_nodes, v_nodes, w_nodes]
    pub velocity: DVector<T>,
    /// Pressure field (1 per node)
    pub pressure: DVector<T>,
    /// Number of total nodes
    pub n_nodes: usize,
    /// Number of corner nodes (which have pressure DOFs in Taylor-Hood)
    pub n_corner_nodes: usize,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> StokesFlowSolution<T> {
    /// Create a new solution
    #[must_use]
    pub fn new(velocity: DVector<T>, pressure: DVector<T>, n_nodes: usize) -> Self {
        let n_corner = pressure.len();
        Self {
            velocity,
            pressure,
            n_nodes,
            n_corner_nodes: n_corner,
        }
    }

    /// Explicit Taylor-Hood constructor
    pub fn new_with_corners(velocity: DVector<T>, pressure: DVector<T>, n_nodes: usize, n_corners: usize) -> Self {
        Self {
            velocity,
            pressure,
            n_nodes,
            n_corner_nodes: n_corners,
        }
    }

    /// Get velocity at node
    #[must_use]
    pub fn get_velocity(&self, node_idx: usize) -> Vector3<T> {
        let base = node_idx;
        let v_offset = self.n_nodes;
        Vector3::new(
            self.velocity[base],
            self.velocity[base + v_offset],
            self.velocity[base + 2 * v_offset],
        )
    }

    /// Get pressure at node
    #[must_use]
    pub fn get_pressure(&self, node_idx: usize) -> T {
        self.pressure[node_idx]
    }

    /// Set velocity at node
    pub fn set_velocity(&mut self, node_idx: usize, vel: &Vector3<T>) {
        let base = node_idx;
        let v_offset = self.n_nodes;
        self.velocity[base] = vel.x;
        self.velocity[base + v_offset] = vel.y;
        self.velocity[base + 2 * v_offset] = vel.z;
    }

    /// Set pressure at node
    pub fn set_pressure(&mut self, node_idx: usize, p: T) {
        self.pressure[node_idx] = p;
    }

    /// Blend this solution with another one (useful for Picard relaxation)
    /// result = self * omega + other * (1 - omega)
    pub fn blend(&self, other: &Self, omega: T) -> Self {
        let one = T::one();
        let one_minus_omega = one - omega;

        let velocity = &self.velocity * omega + &other.velocity * one_minus_omega;
        let pressure = &self.pressure * omega + &other.pressure * one_minus_omega;

        Self::new(velocity, pressure, self.n_nodes)
    }

    /// Interleave velocity and pressure into a single vector for the linear solver
    /// Format: [u_nodes..., v_nodes..., w_nodes..., p_nodes...]
    pub fn interleave(&self) -> DVector<T> {
        // Blocks: [U...V...W... P...]
        let n_vel = self.n_nodes * constants::VELOCITY_COMPONENTS;
        let n_pres = self.n_corner_nodes;
        let mut data = DVector::zeros(n_vel + n_pres);
        
        for i in 0..n_vel {
            data[i] = self.velocity[i];
        }
        for i in 0..n_pres {
            data[n_vel + i] = self.pressure[i];
        }
        data
    }
}
