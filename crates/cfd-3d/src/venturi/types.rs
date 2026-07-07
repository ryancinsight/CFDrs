//! Configuration and solution types for the 3D Venturi solver.

use crate::scalar;
use eunomia::FloatElement;
use eunomia::RealField;
use serde::{Deserialize, Serialize};

// ============================================================================
// Solver Configuration
// ============================================================================

/// Configuration for 3D Venturi solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct VenturiConfig3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Inlet volumetric flow rate [m³/s]
    pub inlet_flow_rate: T,
    /// Inlet pressure \[Pa]
    pub inlet_pressure: T,
    /// Outlet pressure \[Pa]
    pub outlet_pressure: T,

    /// Maximum iterations for nonlinear (Picard) solver
    pub max_nonlinear_iterations: usize,
    /// Convergence tolerance for nonlinear iterations
    pub nonlinear_tolerance: T,

    /// Mesh resolution (axial, transverse)
    pub resolution: (usize, usize),
    /// Whether the Venturi is circular or rectangular
    pub circular: bool,
    /// Channel height \[m] for rectangular cross-sections.
    ///
    /// When `circular = false`, the cross-section is `width × height` where
    /// `width` varies axially (inlet → throat → outlet) and `height` is
    /// constant.  When `None`, falls back to `width × width` (square).
    #[serde(default, skip_serializing_if = "Option::is_none")]
    pub rect_height: Option<T>,
}

impl<T> Default for VenturiConfig3D<T>
where
    T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy,
{
    fn default() -> Self {
        Self {
            inlet_flow_rate: scalar::from_f64::<T>(1e-7),
            inlet_pressure: scalar::from_f64::<T>(100.0),
            outlet_pressure: scalar::zero::<T>(),
            max_nonlinear_iterations: 15,
            nonlinear_tolerance: scalar::from_f64::<T>(1e-4),
            resolution: (60, 10),
            circular: false,
            rect_height: None,
        }
    }
}

// ============================================================================
// Solution Result
// ============================================================================

/// Complete solution to 3D Venturi problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct VenturiSolution3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Inlet mean velocity \[m/s]
    pub u_inlet: T,
    /// Maximum velocity in the throat \[m/s]
    pub u_throat: T,
    /// Inlet pressure \[Pa]
    pub p_inlet: T,
    /// Average pressure in the throat \[Pa]
    pub p_throat: T,
    /// Outlet pressure \[Pa]
    pub p_outlet: T,
    /// Pressure drop from inlet to throat \[Pa]
    pub dp_throat: T,
    /// Net pressure recovery/loss from inlet to outlet \[Pa]
    pub dp_recovery: T,
    /// Pressure coefficient at the throat, scaled by throat dynamic pressure
    pub cp_throat: T,
    /// Pressure recovery coefficient at the outlet, scaled by throat dynamic pressure
    pub cp_recovery: T,
    /// Mass balance error (relative)
    pub mass_error: T,
    /// Face-integrated inlet volumetric flow rate [m³/s]
    pub q_in_face: T,
}

impl<T> VenturiSolution3D<T>
where
    T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy,
{
    /// Create a zero-initialized Venturi solution
    pub fn new() -> Self {
        Self {
            u_inlet: scalar::zero::<T>(),
            u_throat: scalar::zero::<T>(),
            p_inlet: scalar::zero::<T>(),
            p_throat: scalar::zero::<T>(),
            p_outlet: scalar::zero::<T>(),
            dp_throat: scalar::zero::<T>(),
            dp_recovery: scalar::zero::<T>(),
            cp_throat: scalar::zero::<T>(),
            cp_recovery: scalar::zero::<T>(),
            mass_error: scalar::zero::<T>(),
            q_in_face: scalar::zero::<T>(),
        }
    }
}

impl<T> Default for VenturiSolution3D<T>
where
    T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy,
{
    fn default() -> Self {
        Self::new()
    }
}
