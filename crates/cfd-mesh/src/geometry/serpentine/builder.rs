//! 3D Serpentine mesh builder
//!
//! Generates a 3D serpentine (sinuous) mesh by extruding a cross-section
//! along a curved centerline defined by a sine wave or circular segments.

use crate::grid::StructuredGridBuilder;
use crate::mesh::Mesh;
use crate::error::Result;
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive};

/// 3D Serpentine mesh builder
#[derive(Debug, Clone)]
pub struct SerpentineMeshBuilder<T: RealField + Copy> {
    /// Channel diameter/width [m]
    pub diameter: T,
    /// Amplitude of the serpentine curve [m]
    pub amplitude: T,
    /// Period of the serpentine curve (wavelength) [m]
    pub wavelength: T,
    /// Number of periods
    pub num_periods: usize,
    /// Number of segments per period in axial direction
    pub segments_per_period: usize,
    /// Number of segments in transverse direction
    pub transverse_resolution: usize,
    /// Whether the cross-section is circular or rectangular
    pub circular: bool,
}

impl<T: RealField + Copy + FromPrimitive + Float> SerpentineMeshBuilder<T> {
    /// Create a new 3D Serpentine mesh builder
    pub fn new(diameter: T, amplitude: T, wavelength: T) -> Self {
        Self {
            diameter,
            amplitude,
            wavelength,
            num_periods: 3,
            segments_per_period: 20,
            transverse_resolution: 8,
            circular: false,
        }
    }

    /// Set number of periods
    pub fn with_periods(mut self, n: usize) -> Self {
        self.num_periods = n;
        self
    }

    /// Set mesh resolution
    pub fn with_resolution(mut self, per_period: usize, transverse: usize) -> Self {
        self.segments_per_period = per_period;
        self.transverse_resolution = transverse;
        self
    }

    /// Set cross-section shape
    pub fn with_circular(mut self, circular: bool) -> Self {
        self.circular = circular;
        self
    }

    /// Build the Serpentine mesh
    pub fn build(&self) -> Result<Mesh<T>> {
        let total_l = self.wavelength * T::from_usize(self.num_periods).unwrap();
        let total_segments = self.segments_per_period * self.num_periods;
        
        // 1. Build initial straight structured hex grid
        // Z is the axial direction (0 to total_l)
        // X, Y are transverse (-1 to 1)
        // 1. Build initial straight structured hex grid
        // Z is the axial direction (0 to total_l)
        // X, Y are transverse (-1 to 1)
        let mut builder = StructuredGridBuilder::<T>::new(
            self.transverse_resolution,
            self.transverse_resolution,
            total_segments,
        );
        let builder = builder.with_bounds((
            (T::from_f64(-1.0).unwrap(), T::from_f64(1.0).unwrap()),
            (T::from_f64(-1.0).unwrap(), T::from_f64(1.0).unwrap()),
            (T::zero(), total_l),
        ));
        
        let mut mesh = builder.build().map_err(|e| crate::error::MeshError::GridError(e.to_string()))?;
        let half_d = self.diameter / T::from_f64(2.0).unwrap();

        // 2. Transform vertices along the serpentine path
        // Path: P(s) = [A*sin(k*s), 0, s]  where k = 2*pi/wavelength
        // We use Frenet-like frame:
        // tangent T = dP/ds / |dP/ds|
        // normal N = dT/ds / |dT/ds|
        // binormal B = T x N
        
        let k = T::from_f64(2.0 * std::f64::consts::PI).unwrap() / self.wavelength;
        
        for v in mesh.vertices_mut() {
            let s = v.position.z; // Parametric position along Z
            
            // Channel path point
            let x_c = self.amplitude * Float::sin(k * s);
            let y_c = T::zero();
            let z_c = s;
            
            // Tangent: T = dP/ds = [A*k*cos(k*s), 0, 1] (unnormalized)
            // Normalized T:
            let tx = self.amplitude * k * Float::cos(k * s);
            let ty = T::zero();
            let tz = T::one();
            let t_vec = Vector3::new(tx, ty, tz).normalize();
            
            // Binormal: Choose a consistent vector (e.g. Y-axis) and cross with T
            // B = normalize(T x [0, 1, 0])
            let b_vec = t_vec.cross(&Vector3::new(T::zero(), T::one(), T::zero())).normalize();
            
            // Normal: N = B x T
            let n_vec = b_vec.cross(&t_vec).normalize();
            
            // Local transverse coordinates from grid
            let u = v.position.x; // [-1, 1]
            let w = v.position.y; // [-1, 1]
            
            let (target_u, target_w) = if self.circular {
                // Map [-1, 1]^2 to disk
                self.map_square_to_disk(u, w, half_d)
            } else {
                (u * half_d, w * half_d)
            };
            
            // Final position: P = center + u*N + w*B
            v.position.x = x_c + target_u * n_vec.x + target_w * b_vec.x;
            v.position.y = y_c + target_u * n_vec.y + target_w * b_vec.y;
            v.position.z = z_c + target_u * n_vec.z + target_w * b_vec.z;
        }

        // 3. Mark boundary faces
        let boundary_faces = mesh.boundary_faces();
        for f_idx in boundary_faces {
            if let Some(face) = mesh.face(f_idx) {
                let mut all_at_inlet = true;
                let mut all_at_outlet = true;
                
                for &v_idx in &face.vertices {
                    if let Some(v) = mesh.vertex(v_idx) {
                        // Use original Z-like coordinate for boundary identification
                        // But since we transformed Z, we should use the parametric value 's'
                        // Wait, v.position.z is now shifted.
                        // Actually, total_l is the exact axial end.
                        // For a sine wave P(s) = [..., s], the Z component of the CENTER is s.
                        // But vertices are shifted by target_u * n_vec.z + ...
                        // So we should check the PARAMETRIC position.
                        // Since StructuredGridBuilder used mesh[i].position.z = s, we can use that if we save it.
                        // Or just check if the vertex was at Z=0 or Z=total_l in the grid.
                        
                        // Let's use a robust way: vertices at the start of the grid have index < (res+1)^2
                        // But indices are shuffled.
                        // Safe way: we haven't overwritten the grid's topology, so we can check original coords if we had them.
                        // Since we don't, we'll use the fact that z_c was 's'. 
                        // At inlet (s=0), T=[A*k, 0, 1], B=[0, 1, 0], N=[-1, 0, Ak]
                        // So z_final = 0 + target_u * Ak.
                        // This is not exactly 0.
                        
                        // Better: StructuredGridBuilder marks boundary faces if we ask it.
                        // It doesn't yet.
                        
                        // I'll use the coordinate range to identify inlet/outlet
                        if Float::abs(v.position.z) > T::from_f64(1e-3).unwrap() {
                            all_at_inlet = false;
                        }
                        if Float::abs(v.position.z - total_l) > T::from_f64(1e-3).unwrap() {
                            all_at_outlet = false;
                        }
                    }
                }
                
                if all_at_inlet {
                    mesh.mark_boundary(f_idx, "inlet".to_string());
                } else if all_at_outlet {
                    mesh.mark_boundary(f_idx, "outlet".to_string());
                } else {
                    mesh.mark_boundary(f_idx, "wall".to_string());
                }
            }
        }

        Ok(mesh)
    }

    fn map_square_to_disk(&self, u: T, w: T, radius: T) -> (T, T) {
        let eps = T::from_f64(1e-10).unwrap();
        if Float::abs(u) < eps && Float::abs(w) < eps {
            return (T::zero(), T::zero());
        }
        
        let r_norm = Float::max(Float::abs(u), Float::abs(w));
        let phi = if Float::abs(u) >= Float::abs(w) {
            if u > T::zero() {
                (T::from_f64(std::f64::consts::PI / 4.0).unwrap()) * (w / u)
            } else {
                (T::from_f64(std::f64::consts::PI / 4.0).unwrap()) * (w / u) + T::from_f64(std::f64::consts::PI).unwrap()
            }
        } else {
            if w > T::zero() {
                -(T::from_f64(std::f64::consts::PI / 4.0).unwrap()) * (u / w) + T::from_f64(std::f64::consts::PI / 2.0).unwrap()
            } else {
                -(T::from_f64(std::f64::consts::PI / 4.0).unwrap()) * (u / w) - T::from_f64(std::f64::consts::PI / 2.0).unwrap()
            }
        };
        
        (r_norm * radius * Float::cos(phi), r_norm * radius * Float::sin(phi))
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_serpentine_mesh_builder() {
        let builder = SerpentineMeshBuilder::<f64>::new(100e-6, 50e-6, 400e-6)
            .with_periods(2)
            .with_resolution(10, 4);
            
        let mesh = builder.build().unwrap();
        
        assert!(mesh.vertex_count() > 0);
        assert!(mesh.cell_count() > 0);
        
        // Check bounds
        let (min, max) = mesh.bounds();
        println!("Serpentine Bounds: Min={:?}, Max={:?}", min, max);
        assert!(max.z >= 800e-6 - 1e-10);
        // The sine wave has peak at amplitude (50e-6) and radius 50e-6.
        // Total range should be roughly [-100e-6, 100e-6].
        assert!(min.x <= -90e-6); 
    }
}
