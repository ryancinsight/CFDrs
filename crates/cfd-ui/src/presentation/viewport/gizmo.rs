//! Axis gizmo — legacy configuration re-exported from `axis_indicator`.
//!
//! This module re-exports `AxisIndicatorConfig` under the old name for
//! backward compatibility during the Phase 1+2 migration.

pub use super::axis_indicator::AxisIndicatorConfig;

/// Legacy alias for `AxisIndicatorConfig`.
pub type AxisGizmo = AxisIndicatorConfig;
