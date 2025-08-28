//! Mesh operations and transformations

use super::types::Mesh;
use nalgebra::{Matrix3, Point3, RealField, Vector3};
use num_traits::FromPrimitive;

/// Trait for mesh operations
pub trait MeshOperations<T: RealField + Copy> {
    /// Translate the mesh
    fn translate(&mut self, displacement: Vector3<T>);

    /// Rotate the mesh
    fn rotate(&mut self, rotation: Matrix3<T>);

    /// Scale the mesh
    fn scale(&mut self, scale_factor: T);

    /// Compute mesh bounds
    fn bounds(&self) -> (Point3<T>, Point3<T>);

    /// Compute mesh centroid
    fn centroid(&self) -> Point3<T>;
}

/// Mesh transformation operations
pub struct MeshTransform<T: RealField + Copy> {
    _phantom: std::marker::PhantomData<T>,
}

impl<T: RealField + Copy + FromPrimitive> MeshOperations<T> for Mesh<T> {
    fn translate(&mut self, displacement: Vector3<T>) {
        for node in &mut self.nodes {
            *node += displacement;
        }
    }

    fn rotate(&mut self, rotation: Matrix3<T>) {
        let centroid = self.centroid();
        for node in &mut self.nodes {
            let relative = *node - centroid;
            let rotated = rotation * relative;
            *node = centroid + rotated;
        }
    }

    fn scale(&mut self, scale_factor: T) {
        let centroid = self.centroid();
        for node in &mut self.nodes {
            let relative = *node - centroid;
            *node = centroid + relative * scale_factor;
        }
    }

    fn bounds(&self) -> (Point3<T>, Point3<T>) {
        if self.nodes.is_empty() {
            return (Point3::origin(), Point3::origin());
        }

        let mut min = self.nodes[0];
        let mut max = self.nodes[0];

        for node in &self.nodes[1..] {
            for i in 0..3 {
                if node[i] < min[i] {
                    min[i] = node[i];
                }
                if node[i] > max[i] {
                    max[i] = node[i];
                }
            }
        }

        (min, max)
    }

    fn centroid(&self) -> Point3<T> {
        if self.nodes.is_empty() {
            return Point3::origin();
        }

        let sum = self
            .nodes
            .iter()
            .fold(Vector3::zeros(), |acc, p| acc + p.coords);
        let count = T::from_usize(self.nodes.len()).unwrap_or_else(T::one);
        Point3::from(sum / count)
    }
}
