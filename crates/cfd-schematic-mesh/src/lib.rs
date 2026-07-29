//! Millifluidic schematic-to-mesh bridge.
//!
//! Converts a `cfd_schematics::NetworkBlueprint` into watertight
//! boundary-surface `IndexedMesh` objects suitable for CFD simulation and
//! STL/manufacturing output.
//!
//! # Why this lives in CFDrs
//!
//! This bridge previously sat in gaia behind a `cfdrs-integration` feature,
//! which made gaia depend back on `cfd-schematics` and produced a
//! repository-level cycle CFDrs -> gaia -> CFDrs. The integrator owns its
//! bridge to a provider, not the reverse, so the pipeline moved here: gaia
//! keeps general geometry, and millifluidic chip design stays in CFDrs.

#![forbid(unsafe_code)]

pub mod blueprint_mesh;
pub mod constraint;
pub mod region_map;
pub mod scheme_io;
pub mod shell_mesh;
pub mod topology;
pub mod well_plate;

pub use blueprint_mesh::{
    BlueprintMeshPipeline, ChannelVolumeTrace, PipelineConfig, PipelineOutput, PipelineVolumeTrace,
    SegmentCenterline,
};
pub use constraint::{DiameterConstraintError, InletOutletConstraint, WallClearanceConstraint};
pub use region_map::{BoundaryLabelMap, RegionMap};
pub use shell_mesh::{ShellMeshPipeline, ShellPipelineConfig, ShellPipelineOutput};
pub use topology::{NetworkTopology, TopologyClass};
pub use well_plate::SbsWellPlate96;
