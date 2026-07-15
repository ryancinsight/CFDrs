//! Data types for the 3D bifurcation solver.

use crate::scalar;
use eunomia::FloatElement;
use eunomia::RealField;
use serde::{Deserialize, Serialize};

use super::geometry::BifurcationGeometry3D;

/// Configuration for 3D bifurcation solver
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct BifurcationConfig3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Inlet volumetric flow rate [m³/s]
    pub inlet_flow_rate: T,
    /// Inlet pressure \[Pa]
    pub inlet_pressure: T,
    /// Outlet pressure \[Pa]
    pub outlet_pressure: T,

    /// Time step size \[s] (for transient)
    pub time_step: T,
    /// Number of time steps
    pub num_time_steps: usize,
    /// Use steady-state (num_time_steps = 1)
    pub steady_state: bool,

    /// Maximum iterations for nonlinear solver
    pub max_nonlinear_iterations: usize,
    /// Convergence tolerance for nonlinear iterations
    pub nonlinear_tolerance: T,

    /// Maximum iterations for linear solver
    pub max_linear_iterations: usize,
    /// Convergence tolerance for linear solver
    pub linear_tolerance: T,
    /// Base mesh resolution factor for branching surface generation
    pub mesh_resolution: usize,
}

impl<T> Default for BifurcationConfig3D<T>
where
    T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy,
{
    fn default() -> Self {
        Self {
            inlet_flow_rate: scalar::from_f64::<T>(1e-8),
            inlet_pressure: scalar::from_f64::<T>(100.0),
            outlet_pressure: scalar::zero::<T>(),
            time_step: scalar::from_f64::<T>(0.001),
            num_time_steps: 1,
            steady_state: true,
            max_nonlinear_iterations: 20,
            nonlinear_tolerance: scalar::from_f64::<T>(1e-4),
            max_linear_iterations: 1000,
            linear_tolerance: scalar::from_f64::<T>(1e-6),
            mesh_resolution: 8,
        }
    }
}

/// Complete solution to 3D bifurcation problem
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub struct BifurcationSolution3D<T: cfd_mesh::domain::core::Scalar + RealField + Copy> {
    /// Volume flow rate in the parent branch [m³/s]
    pub q_parent: T,
    /// Volume flow rate in the first daughter branch [m³/s]
    pub q_daughter1: T,
    /// Volume flow rate in the second daughter branch [m³/s]
    pub q_daughter2: T,
    /// Mean velocity in the parent branch \[m/s]
    pub u_parent_mean: T,
    /// Mean velocity in the first daughter branch \[m/s]
    pub u_daughter1_mean: T,
    /// Mean velocity in the second daughter branch \[m/s]
    pub u_daughter2_mean: T,
    /// Pressure at the inlet cross-section \[Pa]
    pub p_inlet: T,
    /// Pressure at the junction midpoint \[Pa]
    pub p_junction_mid: T,
    /// Pressure at the first daughter outlet \[Pa]
    pub p_daughter1_outlet: T,
    /// Pressure at the second daughter outlet \[Pa]
    pub p_daughter2_outlet: T,
    /// Mean pressure at the outlet \[Pa]
    pub p_outlet: T,
    /// Pressure drop across the parent branch \[Pa]
    pub dp_parent: T,
    /// Pressure drop across the first daughter branch \[Pa]
    pub dp_daughter1: T,
    /// Pressure drop across the second daughter branch \[Pa]
    pub dp_daughter2: T,
    /// Volume-averaged wall shear stress in the parent branch \[Pa]
    pub wall_shear_stress_parent: T,
    /// Volume-averaged wall shear stress in the first daughter \[Pa]
    pub wall_shear_stress_daughter1: T,
    /// Volume-averaged wall shear stress in the second daughter \[Pa]
    pub wall_shear_stress_daughter2: T,
    /// Relative mass conservation error: |Q_in - Q_d1 - Q_d2| / Q_in
    pub mass_conservation_error: T,
}

impl<T> BifurcationSolution3D<T>
where
    T: cfd_mesh::domain::core::Scalar + RealField + FloatElement + Copy,
{
    /// Create a zero-initialized bifurcation solution for the given geometry
    pub fn new(_geometry: &BifurcationGeometry3D<T>) -> Self {
        Self {
            q_parent: scalar::zero::<T>(),
            q_daughter1: scalar::zero::<T>(),
            q_daughter2: scalar::zero::<T>(),
            u_parent_mean: scalar::zero::<T>(),
            u_daughter1_mean: scalar::zero::<T>(),
            u_daughter2_mean: scalar::zero::<T>(),
            p_inlet: scalar::zero::<T>(),
            p_junction_mid: scalar::zero::<T>(),
            p_daughter1_outlet: scalar::zero::<T>(),
            p_daughter2_outlet: scalar::zero::<T>(),
            p_outlet: scalar::zero::<T>(),
            dp_parent: scalar::zero::<T>(),
            dp_daughter1: scalar::zero::<T>(),
            dp_daughter2: scalar::zero::<T>(),
            wall_shear_stress_parent: scalar::zero::<T>(),
            wall_shear_stress_daughter1: scalar::zero::<T>(),
            wall_shear_stress_daughter2: scalar::zero::<T>(),
            mass_conservation_error: scalar::zero::<T>(),
        }
    }

    /// Check whether mass is conserved within the given tolerance
    pub fn is_mass_conserved(&self, tolerance: T) -> bool {
        self.mass_conservation_error < tolerance
    }
}
