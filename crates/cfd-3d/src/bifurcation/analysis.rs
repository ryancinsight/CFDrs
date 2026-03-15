//! Post-processing analysis methods for 3D bifurcation solutions.
//!
//! Provides boundary flow integration, point pressure extraction, and
//! element-level shear-rate computation using the FEM solution fields.

use tracing;

use cfd_core::conversion::SafeFromF64;
use cfd_core::error::{Error, Result};
use cfd_mesh::domain::core::index::{FaceId, VertexId};
use nalgebra::{RealField, Vector3};
use num_traits::{Float, FromPrimitive, ToPrimitive};

use super::solver::BifurcationSolver3D;

impl<
        T: cfd_mesh::domain::core::Scalar
            + RealField
            + Copy
            + FromPrimitive
            + ToPrimitive
            + SafeFromF64
            + Float
            + From<f64>,
    > BifurcationSolver3D<T>
{
    /// Calculate flow rate through a boundary label using u·n integration (f64 precision)
    ///
    /// # Theorem — Boundary Flow Integration
    ///
    /// The volume flow rate through a surface $S$ is
    /// $Q = \int_S \mathbf{u} \cdot \mathbf{n} \, dA$, approximated by
    /// summing face-averaged fluxes over triangulated boundary faces.
    pub(crate) fn calculate_boundary_flow_f64(
        &self,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
        label: &str,
    ) -> Result<f64> {
        let mut total_q = 0.0_f64;
        let mut face_count = 0;

        for f_idx in 0..mesh.face_count() {
            let f_id = FaceId::from_usize(f_idx);
            if mesh.boundary_label(f_id) == Some(label) {
                face_count += 1;
                let face = mesh.faces.get(f_id);
                // FaceData always has exactly 3 vertices
                let v0 = mesh.vertices.position(face.vertices[0]).coords;
                let v1 = mesh.vertices.position(face.vertices[1]).coords;
                let v2 = mesh.vertices.position(face.vertices[2]).coords;

                let n_vec = (v1 - v0).cross(&(v2 - v0));
                let area = n_vec.norm() * 0.5_f64;
                let face_normal = n_vec.normalize();

                let mut u_avg = nalgebra::Vector3::zeros();
                for &v_id in &face.vertices {
                    let u = solution.get_velocity(v_id.as_usize());
                    u_avg += u;
                    if face_count <= 2 {
                        tracing::debug!(vertex = v_id.as_usize(), ?u, "vertex velocity");
                    }
                }
                u_avg /= 3.0_f64;

                let face_flow = u_avg.dot(&face_normal) * area;
                total_q += face_flow;
            }
        }

        tracing::debug!(
            label,
            face_count,
            ?total_q,
            "Flow integration complete"
        );
        Ok(total_q.abs())
    }

    /// Extract pressure at a point in the mesh (closest node, f64 precision)
    pub(crate) fn extract_point_pressure_f64(
        &self,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
        point: nalgebra::Vector3<f64>,
    ) -> Result<f64> {
        let mut best_node = 0;
        let mut min_dist = f64::MAX;

        let n_pressure_nodes = solution.n_corner_nodes.min(mesh.vertex_count());
        for i in 0..n_pressure_nodes {
            let pos = mesh.vertices.position(VertexId::from_usize(i));
            let dist = (pos.coords - point).norm();
            if dist < min_dist {
                min_dist = dist;
                best_node = i;
            }
        }

        Ok(solution.get_pressure(best_node))
    }

    /// Calculate element shear rate (f64 internal precision)
    ///
    /// # Theorem — Strain Rate Tensor
    ///
    /// The rate-of-deformation tensor is $\varepsilon_{ij} = \frac{1}{2}
    /// (\partial u_i / \partial x_j + \partial u_j / \partial x_i)$, and
    /// the scalar shear rate (second invariant) is
    /// $\dot{\gamma} = \sqrt{2 \varepsilon_{ij} \varepsilon_{ij}}$.
    pub(crate) fn calculate_element_shear_rate_f64(
        &self,
        cell: &cfd_mesh::domain::topology::Cell,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
    ) -> Result<f64> {
        let mut idxs: Vec<usize> = Vec::with_capacity(10);
        for &face_idx in &cell.faces {
            if face_idx < mesh.face_count() {
                let face = mesh.faces.get(FaceId::from_usize(face_idx));
                for &v_id in &face.vertices {
                    let v_usize = v_id.as_usize();
                    if !idxs.contains(&v_usize) {
                        idxs.push(v_usize);
                    }
                }
            }
        }

        if idxs.len() < 4 {
            return Err(Error::Solver("Invalid cell topology".to_string()));
        }

        let mut local_verts = Vec::with_capacity(idxs.len());
        for &idx in &idxs {
            local_verts.push(mesh.vertices.position(VertexId::from_usize(idx)).coords);
        }

        let mut l = nalgebra::Matrix3::zeros();

        if idxs.len() == 10_usize {
            // Tet10 (P2): Evaluate gradient at centroid (L = [0.25, 0.25, 0.25, 0.25])
            use crate::fem::shape_functions::LagrangeTet10;

            // 1. Calculate P1 gradients (∇L_i)
            let mut tet4 = crate::fem::element::FluidElement::<f64>::new(idxs[0..4].to_vec());
            let six_v = tet4.calculate_volume(&local_verts);
            if six_v.abs() < 1e-24_f64 {
                return Ok(0.0_f64);
            }
            tet4.calculate_shape_derivatives(&local_verts[0..4]);
            let p1_grads = nalgebra::Matrix3x4::from_columns(&[
                Vector3::new(
                    tet4.shape_derivatives[(0, 0)],
                    tet4.shape_derivatives[(1, 0)],
                    tet4.shape_derivatives[(2, 0)],
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 1)],
                    tet4.shape_derivatives[(1, 1)],
                    tet4.shape_derivatives[(2, 1)],
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 2)],
                    tet4.shape_derivatives[(1, 2)],
                    tet4.shape_derivatives[(2, 2)],
                ),
                Vector3::new(
                    tet4.shape_derivatives[(0, 3)],
                    tet4.shape_derivatives[(1, 3)],
                    tet4.shape_derivatives[(2, 3)],
                ),
            ]);

            // 2. Evaluate P2 gradients at centroid
            let tet10 = LagrangeTet10::new(p1_grads);
            let l_centroid = [0.25_f64; 4];
            let p2_grads = tet10.gradients(&l_centroid);

            // 3. Compute velocity gradient: L = sum(u_i * ∇N_i)
            for i in 0..10 {
                let u = solution.get_velocity(idxs[i]);
                for row in 0..3 {
                    for col in 0..3 {
                        l[(row, col)] += p2_grads[(col, i)] * u[row];
                    }
                }
            }
        } else {
            // Tet4: constant gradient
            let mut element = crate::fem::element::FluidElement::new(idxs.clone());
            element.calculate_shape_derivatives(&local_verts);
            for (i, &idx) in idxs.iter().enumerate().take(4) {
                let u = solution.get_velocity(idx);
                for row in 0..3 {
                    for col in 0..3 {
                        l[(row, col)] += element.shape_derivatives[(col, i)] * u[row];
                    }
                }
            }
        }

        let epsilon = (l + l.transpose()) * 0.5_f64;
        let mut inner_prod = 0.0_f64;
        for i in 0..3 {
            for j in 0..3 {
                inner_prod += epsilon[(i, j)] * epsilon[(i, j)];
            }
        }

        Ok(num_traits::Float::sqrt(2.0_f64 * inner_prod))
    }
}
