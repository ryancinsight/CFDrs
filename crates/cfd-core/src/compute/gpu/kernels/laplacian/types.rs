//! Types for the 2D Laplacian GPU kernel.

use bytemuck::{Pod, Zeroable};

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
    pub(super) fn as_u32(self) -> u32 {
        match self {
            BoundaryType::Dirichlet => 0,
            BoundaryType::Neumann => 1,
            BoundaryType::Periodic => 2,
        }
    }

    #[allow(missing_docs)]
    #[must_use]
    pub fn description(self) -> &'static str {
        match self {
            BoundaryType::Dirichlet => "Dirichlet: u = 0, odd reflection ghosting for endpoints",
            BoundaryType::Neumann => "Neumann: du/dn = 0, one-sided second derivatives at endpoints",
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

/// Uniform parameters for 2D Laplacian
/// Use 16-byte aligned fields to guarantee consistent WGSL uniform layout.
#[repr(C)]
#[derive(Debug, Clone, Copy, Pod, Zeroable)]
pub(super) struct Laplacian2DUniforms {
    pub dims_bc: [u32; 4], // (nx, ny, bc_type, pad)
    pub inv2: [f32; 4],    // (dx_inv2, dy_inv2, 0.0, 0.0)
}
