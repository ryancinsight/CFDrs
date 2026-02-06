//! 3D Branching geometry and mesh generation
//!
//! Provides builders for creating complex 3D branching meshes
//! (bifurcations, trifurcations) using multi-block structured grids.

use crate::grid::StructuredGridBuilder;
use crate::mesh::Mesh;
use nalgebra::{Point3, RealField, Rotation3, Vector3};
use num_traits::{Float, FromPrimitive};

/// Builder for 3D branching meshes
#[derive(Debug, Clone)]
pub struct BranchingMeshBuilder<T: RealField + Copy + Float> {
    /// Parent branch diameter [m]
    pub d_parent: T,
    /// Parent branch length [m]
    pub l_parent: T,
    /// Daughter branch diameters [m]
    pub d_daughters: Vec<T>,
    /// Daughter branch lengths [m]
    pub l_daughters: Vec<T>,
    /// Branching angles [radians]
    pub angles: Vec<T>,
    /// Mesh resolution (subdivisions)
    pub resolution: usize,
}

impl<T: RealField + Copy + FromPrimitive + Float> BranchingMeshBuilder<T> {
    /// Create a new symmetric bifurcation builder
    pub fn bifurcation(
        d_parent: T,
        l_parent: T,
        d_daughter: T,
        l_daughter: T,
        angle: T,
        resolution: usize,
    ) -> Self {
        Self {
            d_parent,
            l_parent,
            d_daughters: vec![d_daughter, d_daughter],
            l_daughters: vec![l_daughter, l_daughter],
            angles: vec![angle, -angle],
            resolution,
        }
    }

    /// Create a new symmetric trifurcation builder
    pub fn trifurcation(
        d_parent: T,
        l_parent: T,
        d_daughter: T,
        l_daughter: T,
        angle: T,
        resolution: usize,
    ) -> Self {
        Self {
            d_parent,
            l_parent,
            d_daughters: vec![d_daughter, d_daughter, d_daughter],
            l_daughters: vec![l_daughter, l_daughter, l_daughter],
            angles: vec![angle, T::zero(), -angle],
            resolution,
        }
    }

    /// Build the full branching mesh
    ///
    /// # Algorithm
    /// 1. Create parent block
    /// 2. For each daughter:
    ///    a. Create daughter block
    ///    b. Rotate and translate to junction point
    ///    c. Merge with existing mesh
    pub fn build(&self) -> crate::error::Result<Mesh<T>> {
        let mut mesh = self.build_parent()?;

        for i in 0..self.d_daughters.len() {
            let daughter_mesh = self.build_daughter(i)?;
            mesh.merge(&daughter_mesh, T::from_f64(1e-7).unwrap());
        }

        // Mark remaining external faces as walls
        let boundary_faces = mesh.boundary_faces();
        for f_idx in boundary_faces {
            if mesh.boundary_label(f_idx).is_none() {
                mesh.mark_boundary(f_idx, "wall".to_string());
            }
        }

        Ok(mesh)
    }

    fn build_parent(&self) -> crate::error::Result<Mesh<T>> {
        let half_d = self.d_parent / T::from_f64(2.0).unwrap();
        let mut mesh = StructuredGridBuilder::new(self.resolution * 2, self.resolution, self.resolution)
            .with_bounds((
                (T::zero(), self.l_parent),
                (-half_d, half_d),
                (-half_d, half_d),
            ))
            .build()?;

        // Mark inlet (x = 0)
        for f_idx in 0..mesh.face_count() {
            if let Some(face) = mesh.face(f_idx) {
                let mut all_at_zero = true;
                for &v_idx in &face.vertices {
                    if let Some(v) = mesh.vertex(v_idx) {
                        if Float::abs(v.position.x) > T::from_f64(1e-7).unwrap() {
                            all_at_zero = false;
                            break;
                        }
                    }
                }
                if all_at_zero {
                    mesh.mark_boundary(f_idx, "inlet".to_string());
                }
            }
        }

        Ok(mesh)
    }

    fn build_daughter(&self, idx: usize) -> crate::error::Result<Mesh<T>> {
        let d = self.d_daughters[idx];
        let l = self.l_daughters[idx];
        let angle = self.angles[idx];
        let half_d = d / T::from_f64(2.0).unwrap();

        let mut daughter_mesh = StructuredGridBuilder::new(self.resolution * 2, self.resolution, self.resolution)
            .with_bounds((
                (T::zero(), l),
                (-half_d, half_d),
                (-half_d, half_d),
            ))
            .build()?;

        // Mark outlet (x = l) before transformation
        for f_idx in 0..daughter_mesh.face_count() {
            if let Some(face) = daughter_mesh.face(f_idx) {
                let mut all_at_l = true;
                for &v_idx in &face.vertices {
                    if let Some(v) = daughter_mesh.vertex(v_idx) {
                        if Float::abs(v.position.x - l) > T::from_f64(1e-7).unwrap() {
                            all_at_l = false;
                            break;
                        }
                    }
                }
                if all_at_l {
                    daughter_mesh.mark_boundary(f_idx, format!("outlet_{}", idx));
                }
            }
        }

        // Rotate and translate
        // Junction point is at (l_parent, 0, 0)
        let rotation = Rotation3::from_axis_angle(&Vector3::z_axis(), angle);
        let translation = Vector3::new(self.l_parent, T::zero(), T::zero());

        for v in daughter_mesh.vertices_mut() {
            v.position = rotation * v.position + translation;
        }

        Ok(daughter_mesh)
    }
}
