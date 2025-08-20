//! Solution structure for FEM

use nalgebra::{RealField, Vector3, DVector};
use crate::fem::constants;

/// Solution for 3D incompressible flow
#[derive(Debug, Clone)]
pub struct StokesFlowSolution<T: RealField> {
    /// Velocity field (3 components per node)
    pub velocity: DVector<T>,
    /// Pressure field (1 per node)
    pub pressure: DVector<T>,
    /// Number of nodes
    pub n_nodes: usize,
}

impl<T: RealField + Copy> StokesFlowSolution<T> {
    /// Create a new solution
    pub fn new(velocity: DVector<T>, pressure: DVector<T>, n_nodes: usize) -> Self {
        Self {
            velocity,
            pressure,
            n_nodes,
        }
    }

    /// Get velocity at node
    pub fn get_velocity(&self, node_idx: usize) -> Vector3<T> {
        let base = node_idx * constants::VELOCITY_COMPONENTS;
        Vector3::new(
            self.velocity[base],
            self.velocity[base + 1],
            self.velocity[base + 2],
        )
    }

    /// Get pressure at node
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
}