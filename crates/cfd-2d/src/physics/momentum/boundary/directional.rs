//! Directional boundary condition application for momentum equations.
//!
//! Private helpers called by [`super::apply_momentum_boundaries`] — one function
//! per wall of the rectangular domain.

use super::apply_rotating_wall_bc;
use crate::physics::momentum::solver::MomentumComponent;
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_math::sparse::SparseMatrixBuilder;
use nalgebra::RealField;
use num_traits::FromPrimitive;

pub(super) fn apply_west_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
    grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for j in 0..ny {
        let idx = j * nx;

        match bc {
            BoundaryCondition::Dirichlet {
                value,
                component_values,
            } => {
                let val = if let Some(comps) = component_values {
                    let c_idx = match component {
                        MomentumComponent::U => 0,
                        MomentumComponent::V => 1,
                    };
                    if c_idx < comps.len() {
                        comps[c_idx].unwrap_or(*value)
                    } else {
                        *value
                    }
                } else {
                    *value
                };
                rhs[idx] = val;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                let component_idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                rhs[idx] = velocity[component_idx];
            }
            BoundaryCondition::Wall { wall_type } => {
                match wall_type {
                    cfd_core::physics::boundary::WallType::NoSlip => {
                        rhs[idx] = T::zero();
                    }
                    cfd_core::physics::boundary::WallType::Slip => {}
                    cfd_core::physics::boundary::WallType::Moving { velocity } => {
                        let component_idx = match component {
                            MomentumComponent::U => 0,
                            MomentumComponent::V => 1,
                        };
                        rhs[idx] = velocity[component_idx];
                    }
                    cfd_core::physics::boundary::WallType::Rotating { omega, center } => {
                        rhs[idx] = apply_rotating_wall_bc(component, omega, center, grid, idx);
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + 1, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::Periodic { partner: _ } => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + nx - 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Symmetry => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::PressureInlet { .. } | BoundaryCondition::PressureOutlet { .. } => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Outflow => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::CharacteristicInlet { velocity, .. } => {
                if let Some(vel) = velocity {
                    let component_idx = match component {
                        MomentumComponent::U => 0,
                        MomentumComponent::V => 1,
                    };
                    rhs[idx] = vel[component_idx];
                } else {
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::CharacteristicOutlet {
                extrapolate_velocity,
                ..
            } => {
                if *extrapolate_velocity {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + 1, -T::one())?;
                    rhs[idx] = T::zero();
                } else {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + 1, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            _ => {}
        }
    }

    Ok(())
}

pub(super) fn apply_east_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
    grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for j in 0..ny {
        let idx = j * nx + nx - 1;

        match bc {
            BoundaryCondition::Dirichlet {
                value,
                component_values,
            } => {
                let val = if let Some(comps) = component_values {
                    let c_idx = match component {
                        MomentumComponent::U => 0,
                        MomentumComponent::V => 1,
                    };
                    if c_idx < comps.len() {
                        comps[c_idx].unwrap_or(*value)
                    } else {
                        *value
                    }
                } else {
                    *value
                };
                rhs[idx] = val;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                let component_idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                rhs[idx] = velocity[component_idx];
            }
            BoundaryCondition::Wall { wall_type } => {
                match wall_type {
                    cfd_core::physics::boundary::WallType::NoSlip => {
                        rhs[idx] = T::zero();
                    }
                    cfd_core::physics::boundary::WallType::Slip => {}
                    cfd_core::physics::boundary::WallType::Moving { velocity } => {
                        let component_idx = match component {
                            MomentumComponent::U => 0,
                            MomentumComponent::V => 1,
                        };
                        rhs[idx] = velocity[component_idx];
                    }
                    cfd_core::physics::boundary::WallType::Rotating { omega, center } => {
                        rhs[idx] = apply_rotating_wall_bc(component, omega, center, grid, idx);
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx - 1, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::Periodic { partner: _ } => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, j * nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Symmetry => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::PressureInlet { .. } | BoundaryCondition::PressureOutlet { .. } => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Outflow => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - 1, -T::one())?;
                rhs[idx] = T::zero();
            }
            _ => {}
        }
    }

    Ok(())
}

pub(super) fn apply_north_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
    grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for i in 0..nx {
        let idx = (ny - 1) * nx + i;

        match bc {
            BoundaryCondition::Dirichlet {
                value,
                component_values,
            } => {
                let val = if let Some(comps) = component_values {
                    let c_idx = match component {
                        MomentumComponent::U => 0,
                        MomentumComponent::V => 1,
                    };
                    if c_idx < comps.len() {
                        comps[c_idx].unwrap_or(*value)
                    } else {
                        *value
                    }
                } else {
                    *value
                };
                rhs[idx] = val;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                let component_idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                rhs[idx] = velocity[component_idx];
            }
            BoundaryCondition::Wall { wall_type } => {
                match wall_type {
                    cfd_core::physics::boundary::WallType::NoSlip => {
                        rhs[idx] = T::zero();
                    }
                    cfd_core::physics::boundary::WallType::Slip => {}
                    cfd_core::physics::boundary::WallType::Moving { velocity } => {
                        let component_idx = match component {
                            MomentumComponent::U => 0,
                            MomentumComponent::V => 1,
                        };
                        rhs[idx] = velocity[component_idx];
                    }
                    cfd_core::physics::boundary::WallType::Rotating { omega, center } => {
                        rhs[idx] = apply_rotating_wall_bc(component, omega, center, grid, idx);
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx - nx, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::Periodic { partner: _ } => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, i, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Symmetry => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::PressureInlet { .. } | BoundaryCondition::PressureOutlet { .. } => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Outflow => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx - nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            _ => {}
        }
    }

    Ok(())
}

pub(super) fn apply_south_boundary<T: RealField + Copy + FromPrimitive>(
    matrix: &mut SparseMatrixBuilder<T>,
    rhs: &mut nalgebra::DVector<T>,
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
    grid: &crate::grid::StructuredGrid2D<T>,
    nx: usize,
    ny: usize,
) -> cfd_core::error::Result<()> {
    for i in 0..nx {
        let idx = i;

        match bc {
            BoundaryCondition::Dirichlet {
                value,
                component_values,
            } => {
                let val = if let Some(comps) = component_values {
                    let c_idx = match component {
                        MomentumComponent::U => 0,
                        MomentumComponent::V => 1,
                    };
                    if c_idx < comps.len() {
                        comps[c_idx].unwrap_or(*value)
                    } else {
                        *value
                    }
                } else {
                    *value
                };
                rhs[idx] = val;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                let component_idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                rhs[idx] = velocity[component_idx];
            }
            BoundaryCondition::Wall { wall_type } => {
                match wall_type {
                    cfd_core::physics::boundary::WallType::NoSlip => {
                        rhs[idx] = T::zero();
                    }
                    cfd_core::physics::boundary::WallType::Slip => {}
                    cfd_core::physics::boundary::WallType::Moving { velocity } => {
                        let component_idx = match component {
                            MomentumComponent::U => 0,
                            MomentumComponent::V => 1,
                        };
                        rhs[idx] = velocity[component_idx];
                    }
                    cfd_core::physics::boundary::WallType::Rotating { omega, center } => {
                        rhs[idx] = apply_rotating_wall_bc(component, omega, center, grid, idx);
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                if *gradient == T::zero() {
                    matrix.add_entry(idx, idx, T::one())?;
                    matrix.add_entry(idx, idx + nx, -T::one())?;
                    rhs[idx] = T::zero();
                }
            }
            BoundaryCondition::Periodic { partner: _ } => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, (ny - 1) * nx + i, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Symmetry => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::PressureInlet { .. } | BoundaryCondition::PressureOutlet { .. } => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            BoundaryCondition::Outflow => {
                matrix.add_entry(idx, idx, T::one())?;
                matrix.add_entry(idx, idx + nx, -T::one())?;
                rhs[idx] = T::zero();
            }
            _ => {}
        }
    }

    Ok(())
}
