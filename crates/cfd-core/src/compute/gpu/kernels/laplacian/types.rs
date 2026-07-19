//! Types for the 2D Laplacian GPU kernel.

use hephaestus_wgpu::BoundaryCondition;

/// Boundary condition type for Laplacian operator
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BoundaryType {
    /// Dirichlet boundary condition (u = 0)
    Dirichlet,
    /// Neumann boundary condition (du/dn = 0)
    Neumann,
    /// Periodic boundary condition
    Periodic,
}

impl BoundaryType {
    #[allow(missing_docs)]
    #[must_use]
    pub fn description(self) -> &'static str {
        match self {
            BoundaryType::Dirichlet => "Dirichlet: u = 0, odd reflection ghosting for endpoints",
            BoundaryType::Neumann => {
                "Neumann: du/dn = 0, one-sided second derivatives at endpoints"
            }
            BoundaryType::Periodic => "Periodic: endpoint-inclusive wrapping to inner indices",
        }
    }

    #[allow(missing_docs)]
    #[must_use]
    pub fn formulation(self) -> &'static str {
        match self {
            BoundaryType::Dirichlet => "d2u/dx2|i=0 = (-2u0)/dx2; d2u/dx2|i=nx-1 = (-2u_{nx-1})/dx2",
            BoundaryType::Neumann => "d2u/dx2|i=0 = (2u0-5u1+4u2-u3)/dx2 (nx>=4) with mirrored fallback",
            BoundaryType::Periodic => "left neighbor wraps to nx-2, right neighbor wraps to 1 for endpoint-inclusive grids",
        }
    }
}

impl From<BoundaryType> for BoundaryCondition {
    fn from(value: BoundaryType) -> Self {
        match value {
            BoundaryType::Dirichlet => BoundaryCondition::Dirichlet,
            BoundaryType::Neumann => BoundaryCondition::Neumann,
            BoundaryType::Periodic => BoundaryCondition::Periodic,
        }
    }
}

