use super::super::solver::MomentumComponent;
use cfd_core::CfdScalar;
use eunomia::FloatElement;
use leto::geometry::Vector3;

/// Apply rotating wall boundary condition: u_wall = ω × r
/// where r is the position vector from center of rotation
pub(super) fn apply_rotating_wall_bc<T: CfdScalar + Copy + FloatElement>(
    component: MomentumComponent,
    omega: &Vector3<T>,
    center: &Vector3<T>,
    grid: &crate::grid::StructuredGrid2D<T>,
    idx: usize,
) -> T {
    let nx = grid.nx;
    let i = idx % nx;
    let j = idx / nx;

    let cell_center = grid.cell_center(i, j).expect("expected value");
    let x = cell_center[0];
    let y = cell_center[1];

    let r_x = x - center.x;
    let r_y = y - center.y;
    let omega_z = omega.z;

    match component {
        MomentumComponent::U => -omega_z * r_y,
        MomentumComponent::V => omega_z * r_x,
    }
}
