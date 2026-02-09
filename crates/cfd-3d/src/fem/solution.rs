//! Solution structure for FEM

use crate::fem::constants;
use nalgebra::{DVector, RealField, Vector3};

/// Solution for 3D incompressible flow
#[derive(Debug, Clone)]
pub struct StokesFlowSolution<T: RealField + Copy> {
    /// Velocity field (3 components per node)
    pub velocity: DVector<T>,
    /// Pressure field (1 per node)
    pub pressure: DVector<T>,
    /// Number of nodes
    pub n_nodes: usize,
}

impl<T: RealField + Copy> StokesFlowSolution<T> {
    /// Create a new solution
    #[must_use]
    pub fn new(velocity: DVector<T>, pressure: DVector<T>, n_nodes: usize) -> Self {
        Self {
            velocity,
            pressure,
            n_nodes,
        }
    }

    /// Get velocity at node
    #[must_use]
    pub fn get_velocity(&self, node_idx: usize) -> Vector3<T> {
        let base = node_idx * constants::VELOCITY_COMPONENTS;
        Vector3::new(
            self.velocity[base],
            self.velocity[base + 1],
            self.velocity[base + 2],
        )
    }

    /// Get pressure at node
    #[must_use]
    pub fn get_pressure(&self, node_idx: usize) -> T {
        self.pressure[node_idx]
    }

    /// Set velocity at node
    pub fn set_velocity(&mut self, node_idx: usize, vel: &Vector3<T>) {
        let base = node_idx * constants::VELOCITY_COMPONENTS;
        self.velocity[base] = vel.x;
        self.velocity[base + 1] = vel.y;
        self.velocity[base + 2] = vel.z;
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
    /// Format: [u0, v0, w0, p0, u1, v1, w1, p1, ...]
    pub fn interleave(&self) -> DVector<T> {
        let dofs_per_node = constants::VELOCITY_COMPONENTS + 1;
        let mut data = DVector::zeros(self.n_nodes * dofs_per_node);
        for i in 0..self.n_nodes {
            let base = i * dofs_per_node;
            let v_base = i * constants::VELOCITY_COMPONENTS;
            data[base] = self.velocity[v_base];
            data[base + 1] = self.velocity[v_base + 1];
            data[base + 2] = self.velocity[v_base + 2];
            data[base + 3] = self.pressure[i];
        }
        data
    }
}
