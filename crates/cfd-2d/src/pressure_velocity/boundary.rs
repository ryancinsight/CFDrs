use cfd_core::physics::boundary::BoundaryCondition;
use nalgebra::RealField;
use std::collections::HashMap;

pub(crate) fn has_pressure_anchor<T: RealField + Copy>(
    boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
) -> bool {
    boundary_conditions
        .values()
        .any(is_pressure_anchor_boundary)
}

pub(crate) fn is_pressure_anchor_boundary<T: RealField + Copy>(
    boundary_condition: &BoundaryCondition<T>,
) -> bool {
    matches!(
        boundary_condition,
        BoundaryCondition::Dirichlet { .. }
            | BoundaryCondition::PressureInlet { .. }
            | BoundaryCondition::PressureOutlet { .. }
            | BoundaryCondition::CharacteristicOutlet { .. }
            | BoundaryCondition::CharacteristicInlet {
                pressure: Some(_),
                ..
            }
    )
}

pub(crate) fn pressure_neighbor_for_side<T: RealField + Copy>(
    side: &str,
    i: usize,
    j: usize,
    nx: usize,
    ny: usize,
    boundary_conditions: &HashMap<String, BoundaryCondition<T>>,
) -> usize {
    if let Some(BoundaryCondition::Periodic { partner }) = boundary_conditions.get(side) {
        if let Some(index) = periodic_partner_index(side, partner, i, j, nx, ny) {
            return index;
        }
    }

    match side {
        "west" => j * nx + 1,
        "east" => j * nx + nx - 2,
        "south" => nx + i,
        "north" => (ny - 2) * nx + i,
        _ => j * nx + i,
    }
}

fn periodic_partner_index(
    side: &str,
    partner: &str,
    i: usize,
    j: usize,
    nx: usize,
    ny: usize,
) -> Option<usize> {
    match (side, partner) {
        ("west", "east") => Some(j * nx + nx - 1),
        ("east", "west") => Some(j * nx),
        ("south", "north") => Some((ny - 1) * nx + i),
        ("north", "south") => Some(i),
        _ => None,
    }
}
