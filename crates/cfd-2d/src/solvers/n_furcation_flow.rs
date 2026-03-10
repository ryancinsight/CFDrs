//! 2D N-furcation flow solver for branching microfluidic junctions
//!
//! Generalises bifurcation and trifurcation geometries to support
//! arbitrary N-branch junctions (quadfurcation, pentafurcation, etc.).

use super::ns_fvm::{BloodModel, NavierStokesSolver2D, SIMPLEConfig, StaggeredGrid2D};
use cfd_core::error::Result as CfdResult;
use nalgebra::RealField;
use num_traits::{Float, FromPrimitive, ToPrimitive};
use serde::{Deserialize, Serialize};

#[derive(Debug, Clone, Serialize, Deserialize)]
/// One daughter branch in an idealized 2D N-furcation junction.
pub struct BranchGeometry<T: RealField + Copy> {
    /// Branch width measured normal to the centerline.
    pub width: T,
    /// Centerline length from junction to outlet.
    pub length: T,
    /// Branch angle in radians relative to the parent centerline.
    pub angle: T, // radians relative to parent
}

#[derive(Debug, Clone, Serialize, Deserialize)]
/// Parent channel plus daughter-branch geometry for a planar N-furcation.
pub struct NFurcationGeometry<T: RealField + Copy> {
    /// Parent-channel width upstream of the junction.
    pub parent_width: T,
    /// Parent-channel length upstream of the junction.
    pub parent_length: T,
    /// Daughter-branch geometries downstream of the junction.
    pub branches: Vec<BranchGeometry<T>>,
}

impl<T: RealField + Copy + FromPrimitive> NFurcationGeometry<T> {
    /// Build a symmetric N-furcation with evenly spaced daughter-branch angles.
    pub fn new_symmetric(parent_width: T, parent_length: T, daughter_width: T, daughter_length: T, num_branches: usize, spread_angle: T) -> Self {
        let mut branches = Vec::new();
        if num_branches == 1 {
            branches.push(BranchGeometry { width: daughter_width, length: daughter_length, angle: T::zero() });
        } else {
            let half_spread = spread_angle / T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);
            let step = spread_angle / T::from_f64((num_branches - 1) as f64).unwrap_or_else(num_traits::Zero::zero);
            for i in 0..num_branches {
                let angle = half_spread - T::from_f64(i as f64).unwrap_or_else(num_traits::Zero::zero) * step;
                branches.push(BranchGeometry { width: daughter_width, length: daughter_length, angle });
            }
        }
        Self { parent_width, parent_length, branches }
    }

    /// Return whether the point `(x, y)` lies inside the parent channel or any daughter branch.
    pub fn contains(&self, x: T, y: T) -> bool {
        let half_pw = self.parent_width / T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);
        if x >= T::zero() && x <= self.parent_length && y >= -half_pw && y <= half_pw {
            return true;
        }

        for branch in &self.branches {
            if self.in_segment(x, y, self.parent_length, T::zero(), branch.angle, branch.length, branch.width) {
                return true;
            }
        }
        false
    }

    fn in_segment(&self, x: T, y: T, start_x: T, start_y: T, angle: T, length: T, width: T) -> bool {
        let dx = x - start_x;
        let dy = y - start_y;
        let cos_a = angle.cos();
        let sin_a = angle.sin();
        let lx = dx * cos_a + dy * sin_a;
        let ly = -dx * sin_a + dy * cos_a;
        let half_w = width / T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);
        lx >= T::zero() && lx <= length && ly >= -half_w && ly <= half_w
    }

    /// Compute an axis-aligned bounding box `[min_x, max_x, min_y, max_y]`.
    pub fn bounding_box(&self) -> [T; 4] {
        let half_pw = self.parent_width / T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);
        let min_x = T::zero();
        let mut max_x = self.parent_length;
        let mut min_y = -half_pw;
        let mut max_y = half_pw;

        for branch in &self.branches {
            let end_x = self.parent_length + branch.length * branch.angle.cos();
            let end_y = branch.length * branch.angle.sin();
            let half_w = branch.width / T::from_f64(2.0).unwrap_or_else(num_traits::Zero::zero);

            max_x = max_x.max(end_x);
            min_y = min_y.min(end_y - half_w).min(end_y); // Simple approximation
            max_y = max_y.max(end_y + half_w).max(end_y);
        }
        [min_x, max_x, min_y, max_y]
    }
}

/// Structured-grid 2D Navier-Stokes solver specialized to an N-furcation mask.
pub struct NFurcationSolver2D<T: RealField + Copy + Float + FromPrimitive> {
    /// Branching junction geometry used to build the computational mask.
    pub geometry: NFurcationGeometry<T>,
    /// Embedded finite-volume Navier-Stokes solver.
    pub ns_solver: NavierStokesSolver2D<T>,
}

impl<T: RealField + Copy + Float + FromPrimitive + ToPrimitive> NFurcationSolver2D<T> {
    /// Construct a masked 2D solver over the N-furcation bounding box.
    pub fn new(geometry: NFurcationGeometry<T>, blood: BloodModel<T>, density: T, nx: usize, ny: usize, config: SIMPLEConfig<T>) -> Self {
        let bbox = geometry.bounding_box();
        let width = bbox[1] - bbox[0];
        let height = bbox[3] - bbox[2];
        let grid = StaggeredGrid2D::new(nx, ny, width, height);
        let mut ns_solver = NavierStokesSolver2D::new(grid, blood, density, config);

        for i in 0..nx {
            for j in 0..ny {
                let x = ns_solver.grid.x_center(i) + bbox[0];
                let y = ns_solver.grid.y_center(j) + bbox[2];
                ns_solver.field.mask[i][j] = geometry.contains(x, y);
            }
        }
        Self { geometry, ns_solver }
    }

    /// Solve the masked N-furcation flow field for a prescribed inlet velocity.
    pub fn solve(&mut self, u_inlet: T) -> CfdResult<NFurcationSolution<T>> {
        let _solve_res = self.ns_solver.solve(u_inlet).map_err(|e| cfd_core::error::Error::Solver(e.to_string()))?;

        let nx = self.ns_solver.grid.nx;
        let ny = self.ns_solver.grid.ny;
        let dy = self.ns_solver.grid.dy;
        let dx = self.ns_solver.grid.dx;

        let mut q_parent = T::zero();
        for j in 0..ny {
            if self.ns_solver.field.mask[0][j] {
                q_parent += self.ns_solver.field.u[0][j] * dy;
            }
        }

        // To properly attribute flows to N daughters at the outlet boundary,
        // we sum out-flux at the boundaries. If N > 2, it's safer to just sum all right-hand boundary flows, 
        // or upper/lower boundary flows. Since branches can exit right, top, or bottom.
        // let mut q_daughters = vec![T::zero(); self.geometry.branches.len()];
        let mut total_out = T::zero();

        for j in 0..ny {
            if self.ns_solver.field.mask[nx - 1][j] {
                total_out += self.ns_solver.field.u[nx][j] * dy;
            }
        }
        for i in 0..nx {
            if self.ns_solver.field.mask[i][ny - 1] {
                total_out += self.ns_solver.field.v[i][ny] * dx;
            }
            if self.ns_solver.field.mask[i][0] {
                total_out -= self.ns_solver.field.v[i][0] * dx; // v is negative, so -v is out flux
            }
        }

        let mass_balance_error = Float::abs(q_parent - total_out) / q_parent;

        Ok(NFurcationSolution {
            q_parent,
            q_total_out: total_out,
            mass_balance_error,
        })
    }
}

#[derive(Debug, Serialize, Deserialize)]
/// Integral flow-rate summary extracted from an N-furcation solution.
pub struct NFurcationSolution<T: RealField + Copy> {
    /// Total inlet flow through the parent branch.
    pub q_parent: T,
    /// Total outflow summed over all daughter exits.
    pub q_total_out: T,
    /// Relative mass-balance defect `|Q_in - Q_out| / Q_in`.
    pub mass_balance_error: T,
}
