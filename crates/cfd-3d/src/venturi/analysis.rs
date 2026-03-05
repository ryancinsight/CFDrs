//! Post-solve analysis helpers for the 3D Venturi solver.
//!
//! Contains shear-rate computation, divergence diagnostics, and boundary
//! flow-rate integration — all factored out of the main `solver` module to
//! satisfy the 500-line Module Size Rule.

use crate::fem::solver::extract_vertex_indices;
use cfd_core::error::{Error, Result};
use cfd_mesh::domain::core::index::VertexId;
use nalgebra::{RealField, Vector3};
use num_traits::Float;

use super::solver::VenturiSolver3D;

impl<
        T: cfd_mesh::domain::core::Scalar
            + RealField
            + Copy
            + num_traits::FromPrimitive
            + num_traits::ToPrimitive
            + cfd_core::conversion::SafeFromF64,
    > VenturiSolver3D<T>
{
    /// Compute the scalar shear rate $\dot\gamma$ for a single mesh cell.
    ///
    /// # Theorem — Symmetric Strain-Rate Tensor
    ///
    /// The strain-rate tensor $\varepsilon_{ij} = \frac{1}{2}(\partial u_i / \partial x_j
    /// + \partial u_j / \partial x_i)$ satisfies:
    ///
    /// ```text
    /// γ̇ = √(2 ε : ε)
    /// ```
    ///
    /// which reduces to $\dot\gamma = \sqrt{2 \sum_{i,j} \varepsilon_{ij}^2}$.
    ///
    /// **Proof sketch**: By definition of the Frobenius inner product and the
    /// symmetry of $\varepsilon$, the scalar magnitude is $|\varepsilon|_F$
    /// scaled by $\sqrt{2}$, giving the second invariant of the strain-rate
    /// tensor used in generalised-Newtonian constitutive models (Carreau,
    /// Casson, etc.).
    pub(crate) fn calculate_cell_shear_rate_f64(
        &self,
        cell: &cfd_mesh::domain::topology::Cell,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
    ) -> Result<f64> {
        let idxs = extract_vertex_indices(cell, mesh, solution.n_corner_nodes)
            .map_err(|e| Error::Solver(e.to_string()))?;
        let vertex_positions: Vec<Vector3<f64>> = mesh
            .vertices
            .iter()
            .map(|(_, v)| v.position.coords)
            .collect();
        let local_verts: Vec<Vector3<f64>> = idxs.iter().map(|&i| vertex_positions[i]).collect();

        let mut l = nalgebra::Matrix3::zeros();

        if idxs.len() == 4 {
            let mut tet4 = crate::fem::element::FluidElement::<f64>::new(idxs.clone());
            let six_v = tet4.calculate_volume(&local_verts);
            if six_v.abs() < 1e-24_f64 {
                return Ok(0.0_f64);
            }
            tet4.calculate_shape_derivatives(&local_verts[0..4]);

            for k in 0..4 {
                let u = solution.get_velocity(idxs[k]);
                for row in 0..3 {
                    for col in 0..3 {
                        l[(row, col)] += u[row] * tet4.shape_derivatives[(col, k)];
                    }
                }
            }
        } else if idxs.len() == 10 {
            use crate::fem::shape_functions::LagrangeTet10;

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

            let tet10 = LagrangeTet10::new(p1_grads);
            let l_centroid = [0.25_f64; 4];
            let p2_grads = tet10.gradients(&l_centroid);

            for i in 0..10 {
                let u = solution.get_velocity(idxs[i]);
                for row in 0..3 {
                    for col in 0..3 {
                        l[(row, col)] += p2_grads[(col, i)] * u[row];
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
        let shear = (2.0_f64 * inner_prod).sqrt();
        Ok(shear)
    }

    /// Print element-wise velocity divergence statistics.
    ///
    /// For an incompressible flow $\nabla \cdot \mathbf{u} = 0$; the
    /// volume-weighted mean of $|\nabla \cdot \mathbf{u}|$ quantifies how
    /// well the discrete solution satisfies incompressibility.
    pub(crate) fn print_divergence_stats(
        &self,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
    ) -> Result<()> {
        let mut min_div = f64::MAX;
        let mut max_div = 0.0_f64;
        let mut sum_div = 0.0_f64;
        let mut count = 0usize;
        let mut total_volume = 0.0_f64;
        let mut signed_div_sum = 0.0_f64;
        let mut abs_div_vol_sum = 0.0_f64;

        for cell in &mesh.cells {
            let idxs = extract_vertex_indices(cell, mesh, solution.n_corner_nodes)?;
            if idxs.len() < 4 {
                continue;
            }

            let mut local_verts = Vec::with_capacity(idxs.len());
            for &idx in &idxs {
                local_verts.push(mesh.vertices.get(VertexId::from_usize(idx)).position.coords);
            }

            let mut div = 0.0_f64;
            let cell_volume: f64;
            if idxs.len() == 10 {
                let mut tet4 = crate::fem::element::FluidElement::new(idxs[0..4].to_vec());
                let six_v = tet4.calculate_volume(&local_verts);
                if Float::abs(six_v) < 1e-24_f64 {
                    continue;
                }
                cell_volume = tet4.volume;
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

                let tet10 = crate::fem::shape_functions::LagrangeTet10::new(p1_grads);
                let l_centroid = [0.25_f64; 4];
                let p2_grads = tet10.gradients(&l_centroid);

                for i in 0..10 {
                    let u = solution.get_velocity(idxs[i]);
                    div += p2_grads[(0, i)] * u.x
                        + p2_grads[(1, i)] * u.y
                        + p2_grads[(2, i)] * u.z;
                }
            } else {
                let mut element = crate::fem::element::FluidElement::new(idxs.clone());
                element.calculate_volume(&local_verts);
                element.calculate_shape_derivatives(&local_verts);
                cell_volume = element.volume;
                for (i, &idx) in idxs.iter().enumerate().take(4) {
                    let u = solution.get_velocity(idx);
                    div += element.shape_derivatives[(0, i)] * u.x
                        + element.shape_derivatives[(1, i)] * u.y
                        + element.shape_derivatives[(2, i)] * u.z;
                }
            }

            let div_abs = Float::abs(div);
            if div_abs < min_div {
                min_div = div_abs;
            }
            if div_abs > max_div {
                max_div = div_abs;
            }
            sum_div += div_abs;
            total_volume += cell_volume;
            signed_div_sum += div * cell_volume;
            abs_div_vol_sum += div_abs * cell_volume;
            count += 1;
        }

        if count > 0 {
            let avg_div = sum_div / count as f64;
            let vol_avg_div = if total_volume > 0.0_f64 {
                abs_div_vol_sum / total_volume
            } else {
                0.0_f64
            };
            println!(
                "Divergence Stats: min={min_div:?}, max={max_div:?}, avg={avg_div:?}, vol_avg={vol_avg_div:?}, net={signed_div_sum:?} (n={count}, vol={total_volume:?})"
            );
        }

        Ok(())
    }

    /// Compute the face-integrated volumetric flow rate across a labeled
    /// boundary patch (e.g. `"inlet"` or `"outlet"`).
    pub(crate) fn calculate_boundary_flow(
        &self,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
        label: &str,
    ) -> Result<f64> {
        let mut total_q = 0.0_f64;
        let mut face_count = 0usize;

        for f_idx in mesh.boundary_faces() {
            if mesh.boundary_label(f_idx) == Some(label) {
                let face = mesh.faces.get(f_idx);
                if face.vertices.len() >= 3 {
                    face_count += 1;
                    let v0 = mesh.vertices.get(face.vertices[0]).position.coords;
                    let v1 = mesh.vertices.get(face.vertices[1]).position.coords;
                    let v2 = mesh.vertices.get(face.vertices[2]).position.coords;

                    let n_vec = (v1 - v0).cross(&(v2 - v0));
                    let area = n_vec.norm() * 0.5_f64;
                    if area <= 0.0_f64 {
                        continue;
                    }
                    let face_normal = n_vec.normalize();

                    let mut u_avg = Vector3::zeros();
                    for &v_idx in &face.vertices {
                        u_avg += solution.get_velocity(v_idx.as_usize());
                    }
                    u_avg /= face.vertices.len() as f64;

                    let mut n_oriented = face_normal;
                    if (label == "inlet" && n_oriented.z > 0.0_f64)
                        || (label == "outlet" && n_oriented.z < 0.0_f64)
                    {
                        n_oriented = -n_oriented;
                    }

                    let signed_flux = u_avg.dot(&n_oriented) * area;
                    let face_flow = if label == "inlet" {
                        -signed_flux
                    } else {
                        signed_flux
                    };
                    total_q += face_flow;
                }
            }
        }

        println!("Venturi Flux: label={label}, faces={face_count}, total_q={total_q:?}");

        Ok(total_q)
    }

    /// Compute the face-integrated volumetric flow rate crossing a constant-$z$
    /// plane over the mesh.
    pub(crate) fn calculate_plane_flux(
        &self,
        mesh: &cfd_mesh::IndexedMesh<f64>,
        solution: &crate::fem::StokesFlowSolution<f64>,
        z_plane: f64,
        tol: f64,
    ) -> Result<(f64, usize)> {
        let mut total_q = 0.0_f64;
        let mut face_count = 0usize;

        for face in mesh.faces.iter() {
            if face.vertices.len() < 3 {
                continue;
            }

            let mut on_plane = true;
            for &v_idx in &face.vertices {
                let v = mesh.vertices.get(v_idx);
                if Float::abs(v.position.z - z_plane) > tol {
                    on_plane = false;
                    break;
                }
            }
            if !on_plane {
                continue;
            }

            let v0 = mesh.vertices.get(face.vertices[0]).position.coords;
            let v1 = mesh.vertices.get(face.vertices[1]).position.coords;
            let v2 = mesh.vertices.get(face.vertices[2]).position.coords;

            let n_vec = (v1 - v0).cross(&(v2 - v0));
            let area = n_vec.norm() * 0.5_f64;
            if area <= 0.0_f64 {
                continue;
            }
            let face_normal = n_vec.normalize();

            let mut u_avg = Vector3::zeros();
            for &v_idx in &face.vertices {
                u_avg += solution.get_velocity(v_idx.as_usize());
            }
            u_avg /= face.vertices.len() as f64;

            let mut face_flow = u_avg.dot(&face_normal) * area;
            if face_normal.z < 0.0_f64 {
                face_flow = -face_flow;
            }

            total_q += face_flow;
            face_count += 1;
        }

        Ok((total_q, face_count))
    }
}
