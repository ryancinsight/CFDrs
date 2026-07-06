//! Mesh transformation and geometric operations
//!
//! This module provides tools for manipulating the mesh geometry,
//! such as translation, rotation, and scaling.

use super::Mesh;
use eunomia::FloatElement;
use eunomia::RealField;
use leto::{
    geometry::{Point3, Vector3},
    FixedMatrix, FixedVector,
};

/// Trait for geometric operations on meshes
pub trait MeshOperations<T: RealField + Copy> {
    /// Translate the entire mesh by a displacement vector
    fn translate(&mut self, displacement: Vector3<T>);

    /// Rotate the mesh around its centroid using a rotation matrix
    fn rotate(&mut self, rotation: FixedMatrix<T, 3, 3>);

    /// Scale the mesh relative to its centroid
    fn scale(&mut self, scale_factor: T);

    /// Compute the axis-aligned bounding box (AABB) of the mesh
    fn bounds(&self) -> (Point3<T>, Point3<T>);

    /// Compute the geometric centroid of all nodes in the mesh
    fn centroid(&self) -> Point3<T>;
}

impl<T: RealField + Copy + FloatElement> MeshOperations<T> for Mesh<T> {
    fn translate(&mut self, displacement: Vector3<T>) {
        for node in &mut self.nodes {
            *node += displacement;
        }
    }

    fn rotate(&mut self, rotation: FixedMatrix<T, 3, 3>) {
        let centroid = self.centroid();
        for node in &mut self.nodes {
            let relative = *node - centroid;
            let rotated = Vector3::from_array((rotation * FixedVector::new(relative.data)).into_array());
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

        let mut min = self.nodes[0].coords.data;
        let mut max = self.nodes[0].coords.data;

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

        (Point3::from_array(min), Point3::from_array(max))
    }

    fn centroid(&self) -> Point3<T> {
        if self.nodes.is_empty() {
            return Point3::origin();
        }

        let sum = self
            .nodes
            .iter()
            .fold(Vector3::zeros(), |acc, p| acc + p.coords);
        let count = <T as FloatElement>::from_f64(self.nodes.len() as f64);
        Point3::from(sum / count)
    }
}
