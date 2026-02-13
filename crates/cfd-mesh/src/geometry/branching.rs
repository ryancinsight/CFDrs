//! 3D Branching geometry and mesh generation
//!
//! Provides builders for creating complex 3D branching meshes
//! (bifurcations, trifurcations) using multi-block structured grids.

use crate::grid::StructuredGridBuilder;
use crate::mesh::Mesh;
use nalgebra::{RealField, Rotation3, Vector3};
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
    ///    c. Stitch to parent
    pub fn build(&self) -> crate::error::Result<Mesh<T>> {
        let mut mesh = self.build_parent()?;

        for (i, &d_daughter) in self.d_daughters.iter().enumerate() {
            let l_daughter = self.l_daughters[i];
            let angle = self.angles[i];
            
            let mut daughter_mesh = self.build_daughter(i, d_daughter, l_daughter, angle)?;
            
            let parent_interface_label = format!("interface_d{}", i);
            let daughter_inlet_label = "inlet"; 
            
            mesh.stitch(&daughter_mesh, &parent_interface_label, daughter_inlet_label, T::from_f64(1e-7).unwrap())
                .map_err(|e| crate::error::MeshError::InvalidMesh(format!("Stitching failed for daughter {}: {}", i, e)))?;
        }
        
        // Note: We no longer aggressively mark all unlabeled topological boundaries as "wall".
        // If the mesh is shattered (e.g. during HexToTet conversion), this would pin internal flow.
        // Instead, we trust the explicit labels from build_parent and build_daughter.
        // Unlabeled external faces will be treated as symmetry/leaky by the solver unless explicitly handled.

        Ok(mesh)
    }

    fn build_parent(&self) -> crate::error::Result<Mesh<T>> {
        let half_d = self.d_parent / T::from_f64(2.0).unwrap();
        let n_daughters = self.d_daughters.len();
        // Ensure y-resolution is divisible by n_daughters for clean N-way splits.
        // Round up resolution to next multiple of n_daughters.
        let res_y = ((self.resolution + n_daughters - 1) / n_daughters) * n_daughters;
        let res_z = self.resolution;
        
        let mut mesh = StructuredGridBuilder::new(self.resolution * 2, res_y, res_z)
            .with_bounds((
                (T::zero(), self.l_parent),
                (-half_d, half_d),
                (-half_d, half_d),
            ))
            .build()?;

        let d = self.d_parent;
        let band_height = d / T::from_usize(n_daughters).unwrap();

        for f_idx in 0..mesh.face_count() {
            if let Some(face) = mesh.face(f_idx) {
                let mut all_at_zero = true;
                for &v_idx in &face.vertices {
                    if let Some(v) = mesh.vertex(v_idx) {
                        if Float::abs(v.position.x) > T::from_f64(1e-7).unwrap() {
                            all_at_zero = false;
                        }
                    }
                }
                if all_at_zero {
                    mesh.mark_boundary(f_idx, "inlet".to_string());
                    continue;
                }
                
                let mut all_at_l = true;
                let mut center_y = T::zero();
                 for &v_idx in &face.vertices {
                    if let Some(v) = mesh.vertex(v_idx) {
                        if Float::abs(v.position.x - self.l_parent) > T::from_f64(1e-7).unwrap() {
                            all_at_l = false;
                        }
                        center_y = center_y + v.position.y;
                    }
                }
                
                if all_at_l {
                    let v_count = T::from_usize(face.vertices.len()).unwrap();
                    center_y = center_y / v_count;
                    
                    // Assign face to daughter based on y-band.
                    // Daughter 0 gets the topmost band, daughter N-1 the bottommost.
                    // Bands: daughter k occupies y ∈ [half_d - (k+1)*band_height, half_d - k*band_height]
                    let y_from_top = half_d - center_y; // distance from top, in [0, d]
                    let band_f = y_from_top / band_height;
                    let band_f64: f64 = nalgebra::try_convert(Float::floor(band_f)).unwrap_or(0.0);
                    let band = (band_f64 as usize).min(n_daughters - 1);
                    mesh.mark_boundary(f_idx, format!("interface_d{}", band));
                } else {
                    // It's a wall face on the parent cylinder
                    mesh.mark_boundary(f_idx, "wall".to_string());
                }
            }
        }

        Ok(mesh)
    }

    fn build_daughter(&self, idx: usize, _d: T, l: T, angle: T) -> crate::error::Result<Mesh<T>> {
        let parent_half_d = self.d_parent / T::from_f64(2.0).unwrap();
        let n_daughters = self.d_daughters.len();
        // Must use the same adjusted res_y as build_parent for vertex alignment.
        let parent_res_y = ((self.resolution + n_daughters - 1) / n_daughters) * n_daughters;
        let ny_daughter = parent_res_y / n_daughters;
        
        let d = self.d_parent;
        let band_height = d / T::from_usize(n_daughters).unwrap();
        
        // Daughter k occupies y ∈ [half_d - (k+1)*band_height, half_d - k*band_height]
        let y_top = parent_half_d - band_height * T::from_usize(idx).unwrap();
        let y_bot = y_top - band_height;
        
        let mut daughter_mesh = StructuredGridBuilder::new(self.resolution * 2, ny_daughter, self.resolution)
            .with_bounds((
                (T::zero(), l),
                (y_bot, y_top), 
                (-parent_half_d, parent_half_d),
            ))
            .build()?;
            
        for f_idx in 0..daughter_mesh.face_count() {
            if let Some(face) = daughter_mesh.face(f_idx) {
                let mut all_at_zero = true;
                let mut all_at_l = true;
                for &v_idx in &face.vertices {
                    if let Some(v) = daughter_mesh.vertex(v_idx) {
                        if Float::abs(v.position.x) > T::from_f64(1e-7).unwrap() { all_at_zero = false; }
                        if Float::abs(v.position.x - l) > T::from_f64(1e-7).unwrap() { all_at_l = false; }
                    }
                }
                if all_at_zero {
                    daughter_mesh.mark_boundary(f_idx, "inlet".to_string());
                } else if all_at_l {
                    daughter_mesh.mark_boundary(f_idx, format!("outlet_{}", idx));
                } else {
                    // It's a wall face on the side of the daughter branch
                    daughter_mesh.mark_boundary(f_idx, "wall".to_string());
                }
            }
        }

        let translation = Vector3::new(self.l_parent, T::zero(), T::zero());

        for v in daughter_mesh.vertices_mut() {
             let x_local = v.position.x;
             let s = x_local / l;
             let theta = angle * s;
             
             let rotation = Rotation3::from_axis_angle(&Vector3::z_axis(), theta);
             
             v.position = rotation.transform_point(&v.position);
             v.position = v.position + translation;
        }

        Ok(daughter_mesh)
    }
}
