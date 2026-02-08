//! 3D Venturi mesh builder implementation
//!
//! This module provides a high-fidelity builder for 3D Venturi tube meshes.
//! It utilizes a mathematical mapping approach where a structured hexahedral
//! grid is transformed into the desired Venturi profile (inlet, convergent,
//! throat, divergent, and outlet segments).

use crate::grid::StructuredGridBuilder;
use crate::mesh::Mesh;
use crate::error::Result;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive};

/// 3D Venturi mesh builder
#[derive(Debug, Clone)]
pub struct VenturiMeshBuilder<T: RealField + Copy + Float> {
    /// Inlet diameter/width [m]
    pub d_inlet: T,
    /// Throat diameter/width [m]
    pub d_throat: T,
    /// Outlet diameter/width [m]
    pub d_outlet: T,
    /// Inlet section length [m]
    pub l_inlet: T,
    /// Convergent section length [m]
    pub l_convergent: T,
    /// Throat section length [m]
    pub l_throat: T,
    /// Divergent section length [m]
    pub l_divergent: T,
    /// Outlet section length [m]
    pub l_outlet: T,
    /// Axial resolution (number of cells along flow)
    pub n_axial: usize,
    /// Transverse resolution (number of cells across section)
    pub n_transverse: usize,
    /// Whether the cross-section is circular (true) or rectangular (false)
    pub circular: bool,
}

impl<T: RealField + Copy + FromPrimitive + Float> VenturiMeshBuilder<T> {
    /// Create a new 3D Venturi builder with symmetric inlet/outlet
    pub fn new(d_in: T, d_th: T, l_in: T, l_conv: T, l_th: T, l_div: T, l_out: T) -> Self {
        Self {
            d_inlet: d_in,
            d_throat: d_th,
            d_outlet: d_in,
            l_inlet: l_in,
            l_convergent: l_conv,
            l_throat: l_th,
            l_divergent: l_div,
            l_outlet: l_out,
            n_axial: 40,
            n_transverse: 10,
            circular: false, // Default to rectangular for better stability in structured grids
        }
    }

    /// Set axial and transverse resolution
    pub fn with_resolution(mut self, axial: usize, transverse: usize) -> Self {
        self.n_axial = axial;
        self.n_transverse = transverse;
        self
    }

    /// Set whether the cross-section is circular
    pub fn with_circular(mut self, circular: bool) -> Self {
        self.circular = circular;
        self
    }

    /// Calculate total length of the Venturi
    pub fn total_length(&self) -> T {
        self.l_inlet + self.l_convergent + self.l_throat + self.l_divergent + self.l_outlet
    }

    /// Build the Venturi mesh
    pub fn build(&self) -> Result<Mesh<T>> {
        let total_l = self.total_length();

        // 1. Build initial structured hex grid in a unit cube-like domain
        // We use Z as the axial direction
        let mut mesh = StructuredGridBuilder::new(self.n_transverse, self.n_transverse, self.n_axial)
            .with_bounds((
                (-T::one(), T::one()), // normalized x
                (-T::one(), T::one()), // normalized y
                (T::zero(), total_l),  // actual z
            ))
            .build()?;

        // 2. Transform the vertices to follow the Venturi profile
        let eps = T::from_f64(1e-10).unwrap();
        
        for v in mesh.vertices_mut() {
            let z = v.position.z;
            
            // Calculate local radius/half-width based on axial position
            let d_local = if z < self.l_inlet {
                self.d_inlet
            } else if z < self.l_inlet + self.l_convergent {
                let s = (z - self.l_inlet) / self.l_convergent;
                self.d_inlet + (self.d_throat - self.d_inlet) * s
            } else if z < self.l_inlet + self.l_convergent + self.l_throat {
                self.d_throat
            } else if z < self.l_inlet + self.l_convergent + self.l_throat + self.l_divergent {
                let s = (z - (self.l_inlet + self.l_convergent + self.l_throat)) / self.l_divergent;
                self.d_throat + (self.d_outlet - self.d_throat) * s
            } else {
                self.d_outlet
            };

            let r_local = d_local / T::from_f64(2.0).unwrap();
            
            if self.circular {
                // Map the square [-1, 1]^2 to a disk of radius r_local
                // Using Shirley-Chiu mapping or simpler radial mapping
                let x_norm = v.position.x;
                let y_norm = v.position.y;
                
                if Float::abs(x_norm) < eps && Float::abs(y_norm) < eps {
                    v.position.x = T::zero();
                    v.position.y = T::zero();
                } else {
                    let r = Float::max(Float::abs(x_norm), Float::abs(y_norm));
                    let phi = if Float::abs(x_norm) >= Float::abs(y_norm) {
                        if x_norm > T::zero() {
                            (T::from_f64(std::f64::consts::PI / 4.0).unwrap()) * (y_norm / x_norm)
                        } else {
                            (T::from_f64(std::f64::consts::PI / 4.0).unwrap()) * (y_norm / x_norm) + T::from_f64(std::f64::consts::PI).unwrap()
                        }
                    } else {
                        if y_norm > T::zero() {
                            -(T::from_f64(std::f64::consts::PI / 4.0).unwrap()) * (x_norm / y_norm) + T::from_f64(std::f64::consts::PI / 2.0).unwrap()
                        } else {
                            -(T::from_f64(std::f64::consts::PI / 4.0).unwrap()) * (x_norm / y_norm) - T::from_f64(std::f64::consts::PI / 2.0).unwrap()
                        }
                    };
                    
                    v.position.x = r * r_local * Float::cos(phi);
                    v.position.y = r * r_local * Float::sin(phi);
                }
            } else {
                // Rectangular case: simple scaling
                v.position.x = v.position.x * r_local;
                v.position.y = v.position.y * r_local;
            }
        }

        // 3. Mark boundary faces
        // We need to re-scan for boundary faces since StructuredGridBuilder doesn't mark them
        let total_l_val = total_l;
        let boundary_faces = mesh.boundary_faces();
        
        for f_idx in boundary_faces {
            if let Some(face) = mesh.face(f_idx) {
                let mut all_at_inlet = true;
                let mut all_at_outlet = true;
                
                for &v_idx in &face.vertices {
                    if let Some(v) = mesh.vertex(v_idx) {
                        if Float::abs(v.position.z) > T::from_f64(1e-7).unwrap() {
                            all_at_inlet = false;
                        }
                        if Float::abs(v.position.z - total_l_val) > T::from_f64(1e-7).unwrap() {
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
}

#[cfg(test)]
mod tests {
    use super::*;
    use approx::assert_relative_eq;

    #[test]
    fn test_venturi_builder_dimensions() {
        let builder = VenturiMeshBuilder::<f64>::new(10.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0);
        let mesh = builder.with_resolution(10, 4).build().unwrap();
        
        assert_eq!(mesh.cell_count(), 10 * 4 * 4);
        
        // Check bounds
        let (min, max) = mesh.bounds();
        assert_relative_eq!(min.z, 0.0, epsilon = 1e-7);
        assert_relative_eq!(max.z, 25.0, epsilon = 1e-7);
    }
}
