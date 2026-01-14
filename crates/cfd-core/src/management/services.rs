//! Domain services for CFD computations.
//!
//! This module provides domain services that encapsulate complex business logic
//! and coordinate between different domain entities following DDD principles.

// Note: Domain services are being moved to their respective domain modules.
// For example, MeshQualityService is now in geometry::mesh.

pub use crate::geometry::mesh::{MeshOperationsService, MeshQualityService, QualityLevel};
pub use crate::geometry::mesh::{MetricStatistics as QualityStatistics, QualityAssessment};
pub use crate::physics::fluid_dynamics::flow_regimes::FlowRegime;
pub use crate::physics::fluid_dynamics::service::FluidDynamicsService;
