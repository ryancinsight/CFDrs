//! Physical invariant validation for FEM Stokes problems.

use super::problem::StokesFlowProblem;
use cfd_core::error::{Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use nalgebra::{RealField, Vector3};

pub(super) fn validate_physical_invariants<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    problem: &StokesFlowProblem<T>,
) -> Result<()> {
    validate_fluid_properties(problem)?;
    validate_pressure_space(problem)?;
    validate_body_force(problem)?;
    validate_boundary_values(problem)?;
    validate_element_viscosities(problem)
}

fn validate_fluid_properties<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    problem: &StokesFlowProblem<T>,
) -> Result<()> {
    use num_traits::Float;

    if !Float::is_finite(problem.fluid.density) || problem.fluid.density <= T::zero() {
        return Err(Error::InvalidConfiguration(
            "FEM fluid density must be finite and positive".to_string(),
        ));
    }
    if !Float::is_finite(problem.fluid.viscosity) || problem.fluid.viscosity <= T::zero() {
        return Err(Error::InvalidConfiguration(
            "FEM dynamic viscosity must be finite and positive".to_string(),
        ));
    }
    Ok(())
}

fn validate_pressure_space<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    problem: &StokesFlowProblem<T>,
) -> Result<()> {
    if problem.n_corner_nodes == 0 || problem.n_corner_nodes > problem.mesh.vertex_count() {
        return Err(Error::InvalidConfiguration(format!(
            "FEM pressure corner-node count must be in 1..={} but got {}",
            problem.mesh.vertex_count(),
            problem.n_corner_nodes
        )));
    }
    Ok(())
}

fn validate_body_force<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    problem: &StokesFlowProblem<T>,
) -> Result<()> {
    if let Some(force) = problem.body_force {
        validate_vector("FEM body force", force)?;
    }
    Ok(())
}

fn validate_boundary_values<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    problem: &StokesFlowProblem<T>,
) -> Result<()> {
    for (&node, bc) in &problem.boundary_conditions {
        if node >= problem.mesh.vertex_count() {
            return Err(Error::InvalidConfiguration(format!(
                "FEM boundary condition references node {node}, but mesh has {} vertices",
                problem.mesh.vertex_count()
            )));
        }

        match bc {
            BoundaryCondition::Dirichlet {
                value,
                component_values,
            } => {
                validate_scalar("FEM Dirichlet boundary value", *value)?;
                if let Some(values) = component_values {
                    for maybe_value in values.iter().flatten() {
                        validate_scalar("FEM Dirichlet component boundary value", *maybe_value)?;
                    }
                }
            }
            BoundaryCondition::Neumann { gradient } => {
                validate_scalar("FEM Neumann boundary gradient", *gradient)?;
            }
            BoundaryCondition::Robin { alpha, beta, gamma } => {
                validate_scalar("FEM Robin alpha", *alpha)?;
                validate_scalar("FEM Robin beta", *beta)?;
                validate_scalar("FEM Robin gamma", *gamma)?;
            }
            BoundaryCondition::VelocityInlet { velocity } => {
                validate_vector("FEM velocity inlet", *velocity)?;
            }
            BoundaryCondition::PressureInlet {
                pressure,
                velocity_direction,
            } => {
                validate_scalar("FEM pressure inlet", *pressure)?;
                if let Some(direction) = velocity_direction {
                    validate_vector("FEM pressure-inlet direction", *direction)?;
                    if direction.norm() <= T::zero() {
                        return Err(Error::InvalidConfiguration(
                            "FEM pressure-inlet direction must have nonzero magnitude".to_string(),
                        ));
                    }
                }
            }
            BoundaryCondition::PressureOutlet { pressure }
            | BoundaryCondition::CharacteristicOutlet { pressure, .. } => {
                validate_scalar("FEM outlet pressure", *pressure)?;
            }
            BoundaryCondition::MassFlowInlet {
                mass_flow_rate,
                temperature,
            } => {
                validate_scalar("FEM mass-flow inlet", *mass_flow_rate)?;
                if let Some(value) = temperature {
                    validate_scalar("FEM mass-flow inlet temperature", *value)?;
                }
            }
            BoundaryCondition::VolumeFlowInlet { volume_flow_rate } => {
                validate_scalar("FEM volume-flow inlet", *volume_flow_rate)?;
            }
            BoundaryCondition::CharacteristicInlet {
                riemann_invariant_r1,
                riemann_invariant_r2,
                entropy,
                velocity,
                pressure,
            } => {
                for value in [
                    *riemann_invariant_r1,
                    *riemann_invariant_r2,
                    *entropy,
                    *pressure,
                ]
                .into_iter()
                .flatten()
                {
                    validate_scalar("FEM characteristic inlet scalar", value)?;
                }
                if let Some(value) = velocity {
                    validate_vector("FEM characteristic inlet velocity", *value)?;
                }
            }
            BoundaryCondition::Wall { .. }
            | BoundaryCondition::Periodic { .. }
            | BoundaryCondition::Symmetry
            | BoundaryCondition::Outflow => {}
        }
    }

    Ok(())
}

fn validate_element_viscosities<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    problem: &StokesFlowProblem<T>,
) -> Result<()> {
    use num_traits::Float;

    if let Some(viscosities) = &problem.element_viscosities {
        if viscosities.len() != problem.mesh.cells.len() {
            return Err(Error::InvalidConfiguration(format!(
                "FEM element viscosity field length must match cell count: expected {}, got {}",
                problem.mesh.cells.len(),
                viscosities.len()
            )));
        }
        for (idx, &viscosity) in viscosities.iter().enumerate() {
            if !Float::is_finite(viscosity) || viscosity <= T::zero() {
                return Err(Error::InvalidConfiguration(format!(
                    "FEM element viscosity at cell {idx} must be finite and positive"
                )));
            }
        }
    }

    Ok(())
}

fn validate_scalar<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    label: &str,
    value: T,
) -> Result<()> {
    use num_traits::Float;

    if !Float::is_finite(value) {
        return Err(Error::InvalidConfiguration(format!(
            "{label} must be finite"
        )));
    }
    Ok(())
}

fn validate_vector<T: cfd_mesh::domain::core::Scalar + RealField + Copy>(
    label: &str,
    vector: Vector3<T>,
) -> Result<()> {
    validate_scalar(label, vector.x)?;
    validate_scalar(label, vector.y)?;
    validate_scalar(label, vector.z)
}
