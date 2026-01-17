//! Mesh refinement utilities

use crate::error::MeshError;
use crate::mesh::Mesh;
use crate::topology::{Face, Vertex};
use nalgebra::RealField;
use std::collections::HashMap;

/// Helper structure to manage mesh refinement context
pub(crate) struct RefinementContext<'a, T: RealField + Copy> {
    pub(crate) old_mesh: &'a Mesh<T>,
    pub(crate) new_mesh: Mesh<T>,
    pub(crate) midpoint_cache: HashMap<(usize, usize), usize>,
    pub(crate) face_cache: HashMap<Vec<usize>, usize>,
}

impl<'a, T: RealField + Copy> RefinementContext<'a, T> {
    pub fn new(old_mesh: &'a Mesh<T>) -> Self {
        let mut new_mesh = Mesh::new();
        // Copy existing vertices
        for v in old_mesh.vertices() {
            new_mesh.add_vertex(*v);
        }
        Self {
            old_mesh,
            new_mesh,
            midpoint_cache: HashMap::new(),
            face_cache: HashMap::new(),
        }
    }

    pub fn get_midpoint(&mut self, v1_idx: usize, v2_idx: usize) -> Result<usize, MeshError> {
        let key = if v1_idx < v2_idx {
            (v1_idx, v2_idx)
        } else {
            (v2_idx, v1_idx)
        };
        if let Some(&idx) = self.midpoint_cache.get(&key) {
            Ok(idx)
        } else {
            let v1 = self
                .old_mesh
                .vertex(v1_idx)
                .ok_or_else(|| MeshError::InvalidMesh(format!("Invalid vertex index {} in midpoint calculation", v1_idx)))?;
            let v2 = self
                .old_mesh
                .vertex(v2_idx)
                .ok_or_else(|| MeshError::InvalidMesh(format!("Invalid vertex index {} in midpoint calculation", v2_idx)))?;
            let center = nalgebra::center(&v1.position, &v2.position);
            let new_v = Vertex::new(center);
            let idx = self.new_mesh.add_vertex(new_v);
            self.midpoint_cache.insert(key, idx);
            Ok(idx)
        }
    }

    pub fn add_face_dedup(&mut self, v_indices: Vec<usize>) -> usize {
        let mut sorted_indices = v_indices.clone();
        sorted_indices.sort_unstable();

        if let Some(&idx) = self.face_cache.get(&sorted_indices) {
            idx
        } else {
            let face = Face::from_vertices(v_indices);
            let idx = self.new_mesh.add_face(face);
            self.face_cache.insert(sorted_indices, idx);
            idx
        }
    }
}
