use super::super::solver::MomentumComponent;
use crate::scalar;
use cfd_core::error::{BoundaryErrorKind, Error, Result};
use cfd_core::physics::boundary::BoundaryCondition;
use cfd_core::CfdScalar;
use eunomia::FloatElement;
use std::collections::HashMap;
use std::hash::BuildHasher;

fn get_dirichlet_value<T: CfdScalar + Copy + FloatElement>(
    bc: &BoundaryCondition<T>,
    component: MomentumComponent,
) -> Option<T> {
    match bc {
        BoundaryCondition::Dirichlet {
            value,
            component_values,
        } => {
            if let Some(comps) = component_values {
                let idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                if idx < comps.len() {
                    comps[idx].or(Some(*value))
                } else {
                    Some(*value)
                }
            } else {
                Some(*value)
            }
        }
        BoundaryCondition::VelocityInlet { velocity } => {
            let idx = match component {
                MomentumComponent::U => 0,
                MomentumComponent::V => 1,
            };
            Some(velocity[idx])
        }
        BoundaryCondition::Wall { wall_type } => match wall_type {
            cfd_core::physics::boundary::WallType::NoSlip => Some(scalar::zero()),
            cfd_core::physics::boundary::WallType::Moving { velocity } => {
                let idx = match component {
                    MomentumComponent::U => 0,
                    MomentumComponent::V => 1,
                };
                Some(velocity[idx])
            }
            _ => None,
        },
        _ => None,
    }
}

/// Validate boundary condition consistency
pub fn validate_boundary_consistency<T, S>(
    boundaries: &HashMap<String, BoundaryCondition<T>, S>,
    _grid: &crate::grid::StructuredGrid2D<T>,
) -> Result<()>
where
    T: CfdScalar + Copy + FloatElement,
    S: BuildHasher,
{
    for direction in &["north", "south", "east", "west"] {
        if !boundaries.contains_key(*direction) {
            return Err(Error::Boundary(BoundaryErrorKind::InvalidRegion(format!(
                "Missing required boundary: {direction}"
            ))));
        }
    }

    for (name, bc) in boundaries {
        if let BoundaryCondition::Periodic { partner } = bc {
            if !boundaries.contains_key(partner) {
                return Err(Error::Boundary(BoundaryErrorKind::InvalidRegion(format!(
                    "Periodic partner '{partner}' not found for boundary '{name}'"
                ))));
            }
            if let Some(partner_bc) = boundaries.get(partner) {
                if let BoundaryCondition::Periodic { partner: p2 } = partner_bc {
                    if p2 != name {
                        return Err(Error::Boundary(BoundaryErrorKind::InvalidRegion(format!(
                            "Periodic mismatch: '{name}' points to '{partner}', but '{partner}' points to '{p2}'"
                        ))));
                    }
                } else {
                    return Err(Error::Boundary(BoundaryErrorKind::InvalidRegion(format!(
                        "Periodic partner '{partner}' is not periodic"
                    ))));
                }
            }
        }
    }

    let corners = [
        ("west", "south"),
        ("east", "south"),
        ("west", "north"),
        ("east", "north"),
    ];

    for (b1_name, b2_name) in corners {
        let b1 = boundaries
            .get(b1_name)
            .expect("invariant: canonical corner boundary exists after boundary validation");
        let b2 = boundaries
            .get(b2_name)
            .expect("invariant: canonical corner boundary exists after boundary validation");

        for component in [MomentumComponent::U, MomentumComponent::V] {
            if let (Some(v1), Some(v2)) = (
                get_dirichlet_value(b1, component),
                get_dirichlet_value(b2, component),
            ) {
                let diff = scalar::abs(v1 - v2);
                let epsilon = T::default_epsilon() * <T as FloatElement>::from_f64(100.0);
                if diff > epsilon {
                    tracing::debug!(
                        corner = %format!("{b1_name}-{b2_name}"),
                        component = ?component,
                        left = ?v1,
                        right = ?v2,
                        "Corner boundary values differ; retaining both side conditions"
                    );
                }
            }
        }
    }

    Ok(())
}
