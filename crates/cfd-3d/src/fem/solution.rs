//! Solution structure for FEM

use crate::fem::{constants, scalar};
use eunomia::{NumericElement, RealField};
use leto::{Array1, Vector3};
use std::ops::{Index, IndexMut};

/// FEM degree-of-freedom vector backed by Leto array storage.
#[derive(Debug, Clone)]
pub struct FemDofVector<T> {
    data: Array1<T>,
}

impl<T> FemDofVector<T> {
    /// Build a FEM DOF vector from dense row-major values.
    pub fn from_vec(values: Vec<T>) -> Self {
        let len = values.len();
        Self {
            data: Array1::from_shape_vec([len], values)
                .expect("invariant: one-dimensional Leto shape matches vector length"),
        }
    }

    /// Number of scalar DOFs.
    #[must_use]
    pub fn len(&self) -> usize {
        self.data.size()
    }

    /// Whether the vector has no scalar DOFs.
    #[must_use]
    pub fn is_empty(&self) -> bool {
        self.len() == 0
    }

    /// Dense slice view of the DOFs.
    #[must_use]
    pub fn as_slice(&self) -> &[T] {
        self.data
            .as_slice()
            .expect("invariant: FEM DOF vectors are dense one-dimensional Leto arrays")
    }

    /// Mutable dense slice view of the DOFs.
    #[must_use]
    pub fn as_mut_slice(&mut self) -> &mut [T] {
        self.data
            .as_slice_mut()
            .expect("invariant: FEM DOF vectors are dense one-dimensional Leto arrays")
    }

    /// Iterator over scalar DOFs.
    pub fn iter(&self) -> std::slice::Iter<'_, T> {
        self.as_slice().iter()
    }

    /// Consume this vector into dense row-major values.
    pub fn into_vec(self) -> Vec<T>
    where
        T: Clone,
    {
        self.data.into_vec()
    }
}

impl<T: Clone> FemDofVector<T> {
    /// Build from an owned Leto dense vector.
    pub fn from_array(array: Array1<T>) -> Self {
        Self { data: array }
    }

    /// Clone into an owned Leto dense vector.
    pub fn to_array(&self) -> Array1<T> {
        Array1::from_shape_vec([self.len()], self.as_slice().to_vec())
            .expect("invariant: one-dimensional Leto shape matches vector length")
    }
}

impl<T: NumericElement> FemDofVector<T> {
    /// Euclidean norm of the vector.
    #[must_use]
    pub fn norm(&self) -> T {
        let squared = self
            .iter()
            .fold(scalar::zero::<T>(), |acc, &value| acc + value * value);
        <T as NumericElement>::sqrt(squared)
    }

    /// Euclidean norm of `self - other`, returning `None` for shape mismatch.
    #[must_use]
    pub fn difference_norm(&self, other: &Self) -> Option<T> {
        if self.len() != other.len() {
            return None;
        }
        let squared =
            self.iter()
                .zip(other.iter())
                .fold(scalar::zero::<T>(), |acc, (&left, &right)| {
                    let diff = left - right;
                    acc + diff * diff
                });
        Some(<T as NumericElement>::sqrt(squared))
    }
}

impl<T> From<Vec<T>> for FemDofVector<T> {
    fn from(value: Vec<T>) -> Self {
        Self::from_vec(value)
    }
}

impl<T: Clone> From<Array1<T>> for FemDofVector<T> {
    fn from(value: Array1<T>) -> Self {
        Self::from_array(value)
    }
}

impl<T> Index<usize> for FemDofVector<T> {
    type Output = T;

    fn index(&self, index: usize) -> &Self::Output {
        &self.data[[index]]
    }
}

impl<T> IndexMut<usize> for FemDofVector<T> {
    fn index_mut(&mut self, index: usize) -> &mut Self::Output {
        &mut self.data[[index]]
    }
}

impl<'a, T> IntoIterator for &'a FemDofVector<T> {
    type Item = &'a T;
    type IntoIter = std::slice::Iter<'a, T>;

    fn into_iter(self) -> Self::IntoIter {
        self.iter()
    }
}

/// Solution for 3D incompressible flow
#[derive(Debug, Clone)]
pub struct StokesFlowSolution<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Velocity field stored in block order: [u_nodes, v_nodes, w_nodes]
    pub velocity: FemDofVector<T>,
    /// Pressure field (1 per node)
    pub pressure: FemDofVector<T>,
    /// Number of total nodes
    pub n_nodes: usize,
    /// Number of corner nodes (which have pressure DOFs in Taylor-Hood)
    pub n_corner_nodes: usize,
}

impl<T: cfd_mesh::domain::core::Scalar + RealField + Copy> StokesFlowSolution<T> {
    /// Create a new solution
    #[must_use]
    pub fn new(
        velocity: impl Into<FemDofVector<T>>,
        pressure: impl Into<FemDofVector<T>>,
        n_nodes: usize,
    ) -> Self {
        let velocity = velocity.into();
        let pressure = pressure.into();
        let n_corner = pressure.len();
        Self {
            velocity,
            pressure,
            n_nodes,
            n_corner_nodes: n_corner,
        }
    }

    /// Explicit Taylor-Hood constructor
    pub fn new_with_corners(
        velocity: impl Into<FemDofVector<T>>,
        pressure: impl Into<FemDofVector<T>>,
        n_nodes: usize,
        n_corners: usize,
    ) -> Self {
        Self {
            velocity: velocity.into(),
            pressure: pressure.into(),
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

    /// Blend this solution with another one.
    ///
    /// Returns `None` when the solution vectors have mismatched shapes.
    pub fn blend(&self, other: &Self, omega: T) -> Option<Self>
    where
        T: NumericElement,
    {
        if self.velocity.len() != other.velocity.len()
            || self.pressure.len() != other.pressure.len()
        {
            return None;
        }
        let one = scalar::one::<T>();
        let one_minus_omega = one - omega;

        let velocity = self
            .velocity
            .iter()
            .zip(other.velocity.iter())
            .map(|(&left, &right)| left * omega + right * one_minus_omega)
            .collect::<Vec<_>>();
        let pressure = self
            .pressure
            .iter()
            .zip(other.pressure.iter())
            .map(|(&left, &right)| left * omega + right * one_minus_omega)
            .collect::<Vec<_>>();

        Some(Self::new(velocity, pressure, self.n_nodes))
    }

    /// Interleave velocity and pressure into a single vector for the linear solver
    /// Format: [u_nodes..., v_nodes..., w_nodes..., p_nodes...]
    pub fn interleave(&self) -> Array1<T> {
        // Blocks: [U...V...W... P...]
        let n_vel = self.n_nodes * constants::VELOCITY_COMPONENTS;
        let n_pres = self.n_corner_nodes;
        let mut data = Array1::zeros([n_vel + n_pres]);

        for i in 0..n_vel {
            data[i] = self.velocity[i];
        }
        for i in 0..n_pres {
            data[n_vel + i] = self.pressure[i];
        }
        data
    }
}
