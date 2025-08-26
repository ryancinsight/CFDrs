//! FEM solver configuration

use nalgebra::RealField;
use cfd_core::numeric;
use num_traits::FromPrimitive;
use serde::{Deserialize, Serialize};
use super::constants;
// Use ElementType from cfd-core as the single source of truth
pub use cfd_core::domains::mesh_operations::ElementType;
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
impl<T: RealField + FromPrimitive + Copy> Default for FemConfig<T> {
    fn default() -> Self {
        Self {
            base: cfd_core::solver::SolverConfig::default(),
            use_stabilization: true,
            tau: cfd_core::numeric::from_f64(constants::DEFAULT_STABILIZATION)?,
            dt: Some(cfd_core::numeric::from_f64(constants::DEFAULT_TIME_STEP)?),
            reynolds: Some(cfd_core::numeric::from_f64(constants::DEFAULT_REYNOLDS)?),
            element_type: ElementType::Tetrahedron,
            quadrature_order: constants::DEFAULT_QUADRATURE_ORDER,
        }
    }
