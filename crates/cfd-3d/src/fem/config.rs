//! FEM solver configuration

use crate::fem::constants;
use nalgebra::RealField;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};

/// FEM configuration for fluid dynamics
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct FemConfig<T: RealField + Copy> {
    /// Base solver configuration
    pub base: cfd_core::solver::SolverConfig<T>,
    /// Use SUPG/PSPG stabilization
    pub use_stabilization: bool,
    /// Stabilization parameter
    pub tau: T,
    /// Time step (for transient problems)
    pub dt: Option<T>,
    /// Reynolds number (for scaling)
    pub reynolds: Option<T>,
    /// Element type to use
    pub element_type: ElementType,
    /// Quadrature order
    pub quadrature_order: usize,
}

/// Element types supported
#[derive(Debug, Clone, Copy, Serialize, Deserialize)]
pub enum ElementType {
    /// Linear tetrahedron (4 nodes)
    Tet4,
    /// Quadratic tetrahedron (10 nodes)
    Tet10,
    /// Linear hexahedron (8 nodes)
    Hex8,
    /// Quadratic hexahedron (20 nodes)
    Hex20,
}

impl<T: RealField + FromPrimitive + Copy> Default for FemConfig<T> {
    fn default() -> Self {
        Self {
            base: cfd_core::solver::SolverConfig::default(),
            use_stabilization: true,
            tau: T::from_f64(constants::DEFAULT_STABILIZATION).unwrap_or_else(|| T::zero()),
            dt: None,
            reynolds: None,
            element_type: ElementType::Tet4,
            quadrature_order: 2,
        }
    }
}